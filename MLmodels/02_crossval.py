import logging
import numpy as np
import pandas as pd
import joblib, os, json
import argparse

from datetime import datetime
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import roc_auc_score, log_loss, brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier
from scipy.stats import pearsonr
from collections import defaultdict


# CLI Variables 
parser = argparse.ArgumentParser(description="Train model with selectable inputs")
parser.add_argument("--indir", "-i", type=str, required=True,
                    help="Directory to load hyperparameters from (searches inside /models/)")
parser.add_argument("--outdir", "-o", type=str, required=True,
                    help="Directory to save model outputs and fold metrics to inside /gebvs/")
parser.add_argument("--filename", "-f", type=str, required=True,
                    help="CSV filename to load, must be in main/top-level directory")
parser.add_argument("--generation", "-g", type=str, default="all",
                    help='Generation filter: "F0", "F1", "F2", or "all" (default "all")')
parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
args = parser.parse_args()

indir = args.indir
outdir = args.outdir
filename = args.filename
generation = args.generation

logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                    format="%(asctime)s %(levelname)s:%(message)s")
logger = logging.getLogger(__name__)
logger.info("Running with outdir=%s filename=%s generation=%s", outdir, filename, generation)

#  Build model from params 
def build_model(name, params):
    if name == "LR":
        return LogisticRegression(max_iter=1000, solver="saga", **params)
    elif name == "RF":
        return RandomForestClassifier(n_jobs=-1, **params)
    elif name == "GB":
        return XGBClassifier(tree_method="hist", eval_metric="logloss", **params)
    else:
        raise ValueError(f"Unknown model name: {name}")


#  Outer Fold Run 
def run_outer_fold(fold, train_idx, test_idx, X, y, ids, best_params, seed=42):
    logger.info(f"Outer Fold {fold+1}/5")

    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    ids_test = ids[test_idx]

    fold_results = pd.DataFrame({"ID": ids_test, "Status": y_test})
    fold_metrics = []

    for name, params in best_params.items():
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
            "fold": fold+1,
            "model": name,
            "AUC": auc,
            "LogLoss": logloss_val,
            "Brier": brier_val,
            "PearsonR": r_val
        })

    return fold_results, pd.DataFrame(fold_metrics)


#  Main 
def main():
    logger.info(f"Loading data from {filename}...")

    chunksize = 100
    list_of_dataframes = []
    path_to_file = os.path.join("data", filename)
    for df_chunk in pd.read_csv(path_to_file, chunksize=chunksize):
        list_of_dataframes.append(df_chunk)
    df = pd.concat(list_of_dataframes)

    if generation != "all":
        df = df[df["Generation"] == generation]

    ids = df["ID"].values
    ax_columns = [col for col in df.columns if col.startswith("AX")]
    X = df[ax_columns].to_numpy()
    y = df["Status"].to_numpy()

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    #  Load hyperparameters 
    with open(f"MLmodels/models/{indir}/best_hyperparams.json", "r") as f:
        best_params = json.load(f)
        logger.info("Loaded best hyperparameters")

    #  Nested CV 
    logger.info("Running 10x repeated 5-fold CV without calibration")
    all_predictions = defaultdict(list)
    all_metrics = []

    for repeat in range(10):
        logger.info(f"=== Repetition {repeat+1}/10 ===")
        skf_outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42 + repeat)

        for fold, (train_idx, test_idx) in enumerate(skf_outer.split(X, y)):
            fold_results, fold_metrics = run_outer_fold(fold, train_idx, test_idx, X, y, ids, best_params)
            all_metrics.append(fold_metrics)

            # Store predictions for each model
            for _, row in fold_results.iterrows():
                ID = row["ID"]
                for col in fold_results.columns:
                    if col not in ["ID", "Status"]:
                        all_predictions[(ID, col)].append(row[col])

    #  Aggregate predictions (mean and SD per model per animal) 
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
    os.makedirs("MLmodels/gebvs", exist_ok=True)
    prob_df.to_csv(os.path.join("MLmodels/gebvs", f"{outdir}_GEBVs_10foldCV.csv"), index=False)
    logger.info("Saved predicted breeding values")
    logger.info(f"Prediction file shape: {prob_df.shape}")
    logger.info("\n" + str(prob_df.head()))

    metrics_df = pd.concat(all_metrics, axis=0)
    logger.info(f"Metrics file shape: {metrics_df.shape}")
    metrics_df.to_csv(os.path.join("MLmodels/gebvs", f"{outdir}_fold_metrics.csv"), index=False)
    logger.info("Saved fold metrics")


if __name__ == "__main__":
    main()
