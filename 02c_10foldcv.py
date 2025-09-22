import logging
import numpy as np
import pandas as pd
import joblib, os, json

from datetime import datetime
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import roc_auc_score, log_loss, brier_score_loss
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.frozen import FrozenEstimator
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier
from scipy.stats import pearsonr
from collections import defaultdict

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
        brier_val = brier_score_loss(y_test, probs)
        r_val, _ = pearsonr(probs, y_test)
        logger.info(
            f"{name} | Fold {fold+1}: "
            f"AUC={auc:.3f}, LogLoss={logloss_val:.3f}, "
            f"Brier={brier_val:.3f}, PearsonR={r_val:.3f}"
        )
        fold_metrics.append({
            "fold": fold+1,
            "model": name,
            "AUC": auc,
            "LogLoss": logloss_val,
            "Brier": brier_val,
            "PearsonR": r_val
        })

    return fold_results, pd.DataFrame(fold_metrics)


# ------------------- Main -------------------
# filename = "DarpaQCGenoPheno.csv"
filename = "noQC/MeanImputedScaledData.csv"
logger.info(f"Loading data from {filename}...")

def main():
    chunksize = 100
    list_of_dataframes = []
    for df in pd.read_csv(filename, chunksize=chunksize, index_col=0):
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

    today_str = "sep22_allhpt_noqc"

    # ---- Load hyperparameters from Script 1 ----
    with open(f"models/{today_str}/best_hyperparams.json", "r") as f:
        best_params = json.load(f)
        logger.info("Loaded best hyperparameters")

    # ---- Nested CV with calibration ----
    logger.info("Running 5x repeated 10-fold CV with calibration for NOQC + ALL GENERATIONS")
    all_predictions = defaultdict(list) 
    all_metrics = []

    for repeat in range(5):
        logger.info(f"=== Repetition {repeat+1}/5 ===")
        skf_outer = StratifiedKFold(n_splits=10, shuffle=True, random_state=42 + repeat)

        for fold, (train_val_idx, test_idx) in enumerate(skf_outer.split(X, y)):
            fold_results, fold_metrics = run_outer_fold(fold, train_val_idx, test_idx, X, y, ids, best_params)
            all_metrics.append(fold_metrics)

            # Store predictions for each model
            for _, row in fold_results.iterrows():
                ID = row["ID"]
                for col in fold_results.columns:
                    if col not in ["ID", "Status"]:
                        all_predictions[(ID, col)].append(row[col])

    # ---- Aggregate predictions (mean and SD per model per animal) ----
    final_records = []
    unique_ids = np.unique(ids)

    for ID in unique_ids:
        record = {"ID": ID}
        record["Status"] = df.loc[df["ID"] == ID, "Status"].iloc[0]
        for model_name in best_params.keys():
            preds = all_predictions[(ID, model_name)]
            record[model_name] = np.mean(preds)
            record[f"{model_name}_SD"] = np.std(preds, ddof=1)
        final_records.append(record)

    prob_df = pd.DataFrame(final_records).sort_values("ID")

    # Save predictions and fold metrics
    os.makedirs("gebvs", exist_ok=True)
    prob_df.to_csv(os.path.join("gebvs", f"{today_str}_GEBVs_10foldCV.csv"), index=False)
    logger.info("Saved predicted breeding values")
    logger.info(f"Prediction file shape: {prob_df.shape}")
    logger.info("Prediction dataframe preview:")
    logger.info("\n" + str(prob_df.head()))

    metrics_df = pd.concat(all_metrics, axis=0)
    metrics_df.to_csv(os.path.join("gebvs", f"{today_str}_fold_metrics.csv"), index=False)
    logger.info("Saved fold metrics")


if __name__ == "__main__":
    main()
