'''
TODO: Check run logs to make sure logging at each fold is working
'''

import logging
import numpy as np
import pandas as pd
import joblib, os, json

from datetime import datetime
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import roc_auc_score, log_loss
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.frozen import FrozenEstimator
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

# ------------------- Logging Setup -------------------
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO,
    force=True
)
logger = logging.getLogger(__name__)


# ------------------- Build model from params -------------------
def build_model(name, params):
    if name == "LR":
        return LogisticRegression(max_iter=1000, solver="saga", **params)
    elif name == "RF":
        return RandomForestClassifier(n_jobs=-1, **params)
    elif name == "GB":
        return XGBClassifier(tree_method="hist", eval_metric="logloss", **params)
    # elif name == "MLP":
    #     return MLPClassifier(max_iter=200, **params)
    else:
        raise ValueError(f"Unknown model name: {name}")


# ------------------- Outer Fold Run -------------------
def run_outer_fold(fold, train_val_idx, test_idx, X, y, ids, best_params, seed=42):
    logger.info(f"Outer Fold {fold+1}/10")

    X_train_val, X_test = X[train_val_idx], X[test_idx]
    y_train_val, y_test = y[train_val_idx], y[test_idx]
    ids_test = ids[test_idx]

    # Split train_val into train and calibration sets (64/16)
    X_train, X_cal, y_train, y_cal = train_test_split(
        X_train_val, y_train_val, test_size=0.2, stratify=y_train_val, random_state=seed
    )

    fold_results = pd.DataFrame({"ID": ids_test, "Status": y_test})
    fold_metrics = []

    for name, params in best_params.items():
        # Build model with saved hyperparameters
        tuned_model = build_model(name, params)
        tuned_model.fit(X_train, y_train)

        # Wrap frozen estimator for calibration (sigmoid safer for small cal sets)
        frozen = FrozenEstimator(tuned_model)
        calibrated = CalibratedClassifierCV(frozen, method="sigmoid")
        calibrated.fit(X_cal, y_cal)

        # Predict on test set
        probs = calibrated.predict_proba(X_test)[:, 1]
        fold_results[name] = probs

        # Evaluate
        auc = roc_auc_score(y_test, probs)
        logloss_val = log_loss(y_test, probs)
        logger.info(f"{name} | Fold {fold+1}: AUC={auc:.3f}, LogLoss={logloss_val:.3f}")

        fold_metrics.append({
            "fold": fold+1,
            "model": name,
            "AUC": auc,
            "Brier": brier,
            "LogLoss": logloss_val
        })

    return fold_results, pd.DataFrame(fold_metrics)


# ------------------- Main -------------------
def main():
    logger.info("Loading dataset...")
    chunksize = 100
    list_of_dataframes = []
    for df in pd.read_csv("DarpaQCGenoPheno.csv", chunksize=chunksize, index_col=0):
        list_of_dataframes.append(df)
    df = pd.concat(list_of_dataframes)

    ids = df["ID"].values
    ax_columns = [col for col in df.columns if col.startswith("AX")]
    X = df[ax_columns]
    y = df["Status"]
    X = X.to_numpy()
    y = y.to_numpy()

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    today_str = datetime.now().strftime("%b%d").lower()  # e.g., "sep05"

    # ---- Load hyperparameters from Script 1 ----
    with open(f"models/{today_str}/best_hyperparams.json", "r") as f:
        best_params = json.load(f)
        logger.info("Loaded best hyperparameters")

    # ---- Nested CV with calibration ----
    logger.info("Running nested CV with calibration")
    skf_outer = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    all_folds = []
    all_metrics = []

    for fold, (train_val_idx, test_idx) in enumerate(skf_outer.split(X, y)):
        fold_results, fold_metrics = run_outer_fold(fold, train_val_idx, test_idx, X, y, ids, best_params)
        all_folds.append(fold_results)
        all_metrics.append(fold_metrics)

    # Save predictions and fold metrics
    prob_df = pd.concat(all_folds, axis=0).sort_values("ID")
    prob_df.to_csv(f"{today_str}_GEBVs_10foldCV.csv", index=False)
    logger.info("Saved predicted breeding values")

    metrics_df = pd.concat(all_metrics, axis=0)
    metrics_df.to_csv(f"{today_str}_fold_metrics.csv", index=False)
    logger.info("Saved fold metrics")


if __name__ == "__main__":
    main()
