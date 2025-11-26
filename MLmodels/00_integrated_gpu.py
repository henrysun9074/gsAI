import logging
import numpy as np
import pandas as pd
import joblib, os, json
import argparse
from datetime import datetime
from collections import defaultdict
from scipy.stats import pearsonr

from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import make_scorer, roc_auc_score, log_loss, brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical

# setup
parser = argparse.ArgumentParser(description="Train model with selectable inputs")
parser.add_argument("--outdir", "-o", type=str, required=True,
                    help="Directory to save model outputs and fold metrics (will be created if missing inside /gebvs/)")
parser.add_argument("--filename", "-f", type=str, required=True,
                    help="CSV filename to load (path relative to current dir, e.g., 'data/my_file.csv')")
parser.add_argument("--hyperparams-file", "-p", type=str, required=True,
                    help="JSON filename to save all best hyperparameters (will be saved in /models/)")
parser.add_argument("--generation", "-g", type=str, default="all",
                    help='Generation filter: "F0", "F1", "F2", or "all" (default "all")')
parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
args = parser.parse_args()

outdir = args.outdir
filename = args.filename
hyperparams_file = args.hyperparams_file
generation = args.generation

logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                    format="%(asctime)s %(levelname)s:%(message)s")
logger = logging.getLogger(__name__)
logger.info("Running with outdir=%s filename=%s hyperparams_file=%s generation=%s", outdir, filename, hyperparams_file, generation)

# custom pearson correlation scorer
def pearson_corr_func(y_true, y_pred):
    """Return Pearson correlation between predicted probs and true labels"""
    if y_pred.ndim == 2:
        y_pred = y_pred[:, 1]
    if np.allclose(y_pred, y_pred[0]):
        return 0.0
    corr, _ = pearsonr(y_pred, y_true)
    return corr if not np.isnan(corr) else 0.0

pearson_scorer = make_scorer(
    pearson_corr_func,
    response_method="predict_proba"
)

# search spaces
def get_search_spaces():
    return {
        "LR": (
            LogisticRegression(max_iter=1000, solver="saga", random_state=88),
            {
                "C": Real(1e-5, 10, prior="log-uniform"),
                "penalty": Categorical(["l1", "l2"]),
            },
        ),
        "RF": (
            RandomForestClassifier(n_jobs=-1, random_state=88),
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
                tree_method="hist", device="cuda", eval_metric="logloss", random_state=88
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
    }

# --- Hyperparameter Tuning ---

def tune_model(X_train, y_train, model_name, n_iter=100):
    base_model, search_space = get_search_spaces()[model_name]
    logger.info(f"Starting Bayesian optimization for {model_name} on training data...")

    opt = BayesSearchCV(
        estimator=base_model,
        search_spaces=search_space,
        n_iter=n_iter,
        cv=5, # 5-fold CV within the training set for tuning
        scoring=pearson_scorer,
        n_jobs=-1,
        verbose=0,
        random_state=88
    )
    opt.fit(X_train, y_train)
    logger.info(f"Best {model_name} params: {opt.best_params_}")
    
    # Extract the best model's estimator and parameters
    best_estimator = opt.best_estimator_
    best_params = opt.best_params_
    
    # Clean up non-serializable objects from best_params if necessary
    final_params = {k: v for k, v in best_params.items()}

    return final_params, best_estimator

def build_model(name, params):
    """Builds a model instance from name and parameters."""
    if name == "LR":
        return LogisticRegression(max_iter=1000, solver="saga", random_state=88, **params)
    elif name == "RF":
        return RandomForestClassifier(n_jobs=-1, random_state=88, **params)
    elif name == "GB":
        return XGBClassifier(tree_method="hist", eval_metric="logloss", random_state=88, **params)
    else:
        raise ValueError(f"Unknown model name: {name}")

def run_outer_fold(fold, train_idx, test_idx, X, y, ids, best_params_all, seed=88):
    """Runs a single fold of the final evaluation CV (Training + Validation data)."""
    logger.info(f"Outer Fold {fold+1}/5")

    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    ids_test = ids[test_idx]

    fold_results = pd.DataFrame({"ID": ids_test, "Status": y_test})
    fold_metrics = []

    for name, params in best_params_all.items():
        # Build model with saved hyperparameters
        tuned_model = build_model(name, params)
        tuned_model.fit(X_train, y_train)

        # Predict on test set
        probs = tuned_model.predict_proba(X_test)[:, 1]
        fold_results[name] = probs

        # Evaluate
        auc = roc_auc_score(y_test, probs)
        logloss_val = log_loss(y_test, probs)
        brier_val = brier_score_loss(y_test, probs)
        r_val, _ = pearsonr(probs, y_test)
        logger.info(
            f"{name} | Fold {fold+1}: "
            f"AUC={auc:.3f}, LogLoss={logloss_val:.3f}, "
            f"Brier={brier_val:.3f}, PearsonR={r_val:.3f}"
        )
        fold_metrics.append({
            "fold": fold + 1,
            "model": name,
            "AUC": auc,
            "LogLoss": logloss_val,
            "Brier": brier_val,
            "PearsonR": r_val
        })

    return fold_results, pd.DataFrame(fold_metrics)


def main():
    chunksize = 100
    list_of_dataframes = []
    logger.info(f"Loading data from {filename}...")
    
    path_to_file = os.path.join("data", filename)
    if not os.path.exists(path_to_file):
        path_to_file = filename 
    
    for df in pd.read_csv(path_to_file, chunksize=chunksize):
        list_of_dataframes.append(df)
    df = pd.concat(list_of_dataframes)
    
    if generation != "all":
        df = df[df["Generation"] == generation]

    ids_full = df["ID"].to_numpy()
    ax_columns = [col for col in df.columns if col.startswith("AX")]
    X_full = df[ax_columns].to_numpy()
    y_full = df["Status"].to_numpy()
    
    scaler = StandardScaler()
    X_full = scaler.fit_transform(X_full)
    
    X_tune_cv, _, y_tune_cv, _, ids_tune_cv, _ = train_test_split(
        X_full, y_full, ids_full, test_size=0.20, stratify=y_full, random_state=88
    )
    
    X_train, X_val, y_train, y_val = train_test_split(
        X_tune_cv, y_tune_cv, test_size=0.20, stratify=y_tune_cv, random_state=88
    )

    logger.info(f"Data Split: Train={X_train.shape[0]} (64%), Validation={X_val.shape[0]} (16%), Test={X_full.shape[0] - X_tune_cv.shape[0]} (20%)")
    
    X_cv, y_cv, ids_cv = X_tune_cv, y_tune_cv, ids_tune_cv
        
    tuned_params_all = {}
    validation_scores = {}
    
    for model_name in get_search_spaces().keys():
        best_params, _ = tune_model(X_train, y_train, model_name, n_iter=100)
        tuned_params_all[model_name] = best_params
        
        final_model = build_model(model_name, best_params)
        final_model.fit(X_train, y_train)
        
        val_probs = final_model.predict_proba(X_val)[:, 1]
        val_score = pearson_corr_func(y_val, val_probs)
        validation_scores[model_name] = val_score
        
        logger.info(f"{model_name} Validation PearsonR: {val_score:.4f}")

    #save hyperparams    
    hyperparams_dir = os.path.join("MLmodels", "models")
    os.makedirs(hyperparams_dir, exist_ok=True)

    hyperparams_path = os.path.join(hyperparams_dir, hyperparams_file)
    with open(hyperparams_path, "w") as f:
        json.dump(tuned_params_all, f, indent=4)
    logger.info(f"Saved all best hyperparameters to {hyperparams_path}")

    #determine the overall best model based on Validation score
    best_model_name = max(validation_scores, key=validation_scores.get)
    logger.info(f"\nOverall best model from validation set: **{best_model_name}** (PearsonR: {validation_scores[best_model_name]:.4f})\n")
        
    logger.info("Starting final evaluation: 10x repeated 5-fold CV on Training+Validation Data (80% of total data)")
    all_predictions = defaultdict(list)
    all_metrics = []

    for repeat in range(10):
        logger.info(f"=== Repetition {repeat+1}/10 ===")
        skf_outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=88 + repeat)

        for fold, (train_idx, test_idx) in enumerate(skf_outer.split(X_cv, y_cv)):
            fold_results, fold_metrics = run_outer_fold(fold, train_idx, test_idx, X_cv, y_cv, ids_cv, tuned_params_all)
            all_metrics.append(fold_metrics)

            # Store predictions for GEBV calculation
            for _, row in fold_results.iterrows():
                ID = row["ID"]
                for col in fold_results.columns:
                    if col not in ["ID", "Status"]:
                        all_predictions[(ID, col)].append(row[col])

    # save results
    final_records = []
    unique_ids = np.unique(ids_cv)

    for ID in unique_ids:
        record = {"ID": ID}
        record["Status"] = df.loc[df["ID"] == ID, "Status"].iloc[0]
        for model_name in tuned_params_all.keys():
            preds = all_predictions.get((ID, model_name), [])
            if preds: 
                record[model_name] = np.mean(preds)
                record[f"{model_name}_SD"] = np.std(preds, ddof=1)
            else:
                record[model_name] = np.nan
                record[f"{model_name}_SD"] = np.nan
        final_records.append(record)

    prob_df = pd.DataFrame(final_records).sort_values("ID")
    
    # Save predictions and fold metrics
    gebvs_dir = os.path.join("MLmodels", "gebvs", outdir)
    os.makedirs(gebvs_dir, exist_ok=True)
    
    prob_df.to_csv(os.path.join(gebvs_dir, f"{outdir}_GEBVs_10foldCV.csv"), index=False)
    logger.info("Saved predicted breeding values")
    logger.info(f"Prediction file shape: {prob_df.shape}")

    metrics_df = pd.concat(all_metrics, axis=0)
    logger.info(f"Metrics file shape: {metrics_df.shape}")
    metrics_df.to_csv(os.path.join(gebvs_dir, f"{outdir}_fold_metrics.csv"), index=False)
    logger.info("Saved fold metrics")


if __name__ == "__main__":
    main()