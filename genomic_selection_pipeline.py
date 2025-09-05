'''
### TODO: Fix section where the AUC is not printing for each fold
### TODO: Fix model saving in correct directory
### TODO: Suppress warnings in script
'''

import logging
import numpy as np
import pandas as pd
import joblib, os
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F

from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import roc_auc_score, brier_score_loss, log_loss
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.frozen import FrozenEstimator
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from xgboost import XGBClassifier

from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical

from joblib import Parallel, delayed

# ------------------- Logging Setup -------------------
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


# ------------------- Model Search Spaces -------------------
def get_search_spaces():
    return {
        "LR": (
            LogisticRegression(max_iter=1000, solver="saga"),
            {
                "C": Real(1e-5, 10, prior="log-uniform"),
                "penalty": Categorical(["l1", "l2"]),
                # "solver": Categorical(["liblinear", "saga", "lbfgs", "newton-cg"]), # invalid combos w/ L1 and L2
            },
        ),
        "RF": (
            RandomForestClassifier(n_jobs=-1),
            {
                "n_estimators": Integer(100, 2000),
                "max_depth": Integer(3, 50),
                "max_features": Categorical(["sqrt", "log2"]),
                "min_samples_split": Integer(2, 20),
                "min_samples_leaf": Integer(1, 10),
            },
        ),
        "GB": (
            XGBClassifier(
                tree_method="hist", device="cuda", eval_metric="logloss"
            ),
            {
                "n_estimators": Integer(100, 2000),
                "max_depth": Integer(3, 15),
                "learning_rate": Real(1e-3, 0.3, prior="log-uniform"),
                "subsample": Real(0.5, 1.0),
                "colsample_bytree": Real(0.5, 1.0),
                "min_child_weight": Integer(1, 10),
                "gamma": Real(0, 5),
            },
        ),
        # "MLP": (
        #     MLPClassifier(max_iter=200),
        #     {
        #         "hidden_layer_sizes": Categorical(
        #             [(64,), (128,), (256,), (128, 64), (256, 128), (128, 128, 64)]
        #         ),
        #         "alpha": Real(1e-6, 1e-1, prior="log-uniform"),
        #         "learning_rate_init": Real(1e-4, 1e-1, prior="log-uniform"),
        #         "batch_size": Integer(32, 512),
        #         "max_iter": Integer(100, 1000),
        #     },
        # ),
    }


# ------------------- Hyperparameter Tuning -------------------
def tune_model(X, y, model_name, n_iter=50):
    base_model, search_space = get_search_spaces()[model_name]
    logger.info(f"Starting Bayesian optimization for {model_name}...")

    opt = BayesSearchCV(
        estimator=base_model,
        search_spaces=search_space,
        n_iter=n_iter,
        cv=5,
        scoring="roc_auc",
        n_jobs=-1,
        verbose=0,
    )
    opt.fit(X, y)
    logger.info(f"Best {model_name} params: {opt.best_params_}")
    return opt.best_estimator_


# ------------------- Outer Fold Run -------------------
def run_outer_fold(fold, train_val_idx, test_idx, X, y, ids, tuned_models, seed=42):
    logger.info(f"Outer Fold {fold+1}")

    X_train_val, X_test = X[train_val_idx], X[test_idx]
    y_train_val, y_test = y[train_val_idx], y[test_idx]
    ids_test = ids[test_idx]

    # Split train_val into train and calibration sets (64/16)
    X_train, X_cal, y_train, y_cal = train_test_split(
        X_train_val, y_train_val, test_size=0.2, stratify=y_train_val, random_state=seed
    )

    fold_results = {"ID": ids_test}

    for name, tuned_model in tuned_models.items():
        # Train on train split
        tuned_model.fit(X_train, y_train)

        # Wrap frozen estimator for calibration
        frozen = FrozenEstimator(tuned_model)
        calibrated = CalibratedClassifierCV(frozen)
        calibrated.fit(X_cal, y_cal)

        # Predict on test set
        probs = calibrated.predict_proba(X_test)[:, 1]
        fold_results[name] = probs

        # Evaluate (log)
        auc = roc_auc_score(y_test, probs)
        brier = brier_score_loss(y_test, probs)
        logloss = log_loss(y_test, probs)
        logger.info(f"{name} | Fold {fold+1}: AUC={auc:.3f}, Brier={brier:.3f}, LogLoss={logloss:.3f}")

    return pd.DataFrame(fold_results)


# ------------------- Nested CV -------------------
def nested_cv(X, y, ids, tuned_models, outer_folds=10, seed=42, n_jobs=-1):
    skf_outer = StratifiedKFold(n_splits=outer_folds, shuffle=True, random_state=seed)

    fold_dfs = Parallel(n_jobs=n_jobs)(
        delayed(run_outer_fold)(fold, train_val_idx, test_idx, X, y, ids, tuned_models, seed)
        for fold, (train_val_idx, test_idx) in enumerate(skf_outer.split(X, y))
    )

    return pd.concat(fold_dfs, axis=0).sort_values("ID")


# ------------------- Main -------------------

def main():
    logger.info("Loading dataset...")
    chunksize = 100
    list_of_dataframes = []
    for df in pd.read_csv('DarpaQCGenoPheno.csv', chunksize=chunksize, index_col=0):
        list_of_dataframes.append(df)
    df = pd.concat(list_of_dataframes)

    ids = df["ID"].values
    ax_columns = [col for col in df.columns if col.startswith('AX')]
    X = df[ax_columns]
    y = df["Status"]
    X = X.to_numpy()
    y = y.to_numpy()

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    # ---- Hyperparameter tuning ----
    tuned_models = {}
    for model_name in get_search_spaces().keys():
        tuned_models[model_name] = tune_model(X, y, model_name)

    # ---- Nested CV with calibration ----
    logger.info("Running nested CV with calibration...")
    prob_df = nested_cv(X, y, ids, tuned_models, n_jobs=-1)

    prob_df.to_csv("GEBVs_10foldCV.csv", index=False)
    logger.info("Predicted calibrated probabilities saved")
    
    # ---- Save trained models ----
    os.makedirs("models", exist_ok=True)
    for name, model in tuned_models.items():
        path = os.path.join("models/sep4", f"{name}_model_10foldCV.joblib")
        joblib.dump(model, path)
        logger.info(f"Saved {name} model to {path}")


if __name__ == "__main__":
    main()
