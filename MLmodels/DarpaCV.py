import logging
import numpy as np
import pandas as pd
import joblib, os, json
import argparse
from datetime import datetime
from collections import defaultdict
from scipy.stats import pearsonr

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import make_scorer, roc_auc_score, log_loss, brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical

# setup
parser = argparse.ArgumentParser(description="Train model with nested CV")
parser.add_argument("--outdir", "-o", type=str, required=True,
                    help="Directory to save model outputs (will be created inside /gebvs/)")
parser.add_argument("--filename", "-f", type=str, required=True,
                    help="CSV filename to load")
parser.add_argument("--generation", "-g", type=str, default="all",
                    help='Generation filter')
parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
args = parser.parse_args()

outdir = args.outdir
filename = args.filename
generation = args.generation

logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                    format="%(asctime)s %(levelname)s:%(message)s")
logger = logging.getLogger(__name__)

# custom pearson correlation scorer
def pearson_corr_func(y_true, y_pred):
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
        # "LR": (
        #     LogisticRegression(max_iter=1000, solver="saga", random_state=88),
        #     {
        #         "C": Real(1e-5, 10, prior="log-uniform"),
        #         "penalty": Categorical(["l1", "l2"]),
        #     },
        # ),
        # "RF": (
        #     RandomForestClassifier(n_jobs=-1, random_state=88),
        #     {
        #         "n_estimators": Integer(100, 2000),
        #         "max_depth": Integer(3, 50),
        #         "max_features": Categorical(["sqrt", "log2"]),
        #         "min_samples_split": Integer(2, 20),
        #         "min_samples_leaf": Integer(1, 10),
        #     },
        # ),
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

def tune_and_train_model(X_train, y_train, model_name, n_iter=20):
    base_model, search_space = get_search_spaces()[model_name]
    
    opt = BayesSearchCV(
        estimator=base_model,
        search_spaces=search_space,
        n_iter=n_iter, 
        cv=5,          # 5-fold CV for the tuning process
        scoring=pearson_scorer,
        n_jobs=-1,
        verbose=0,
        random_state=88
    )
    opt.fit(X_train, y_train)
    
    # Return best estimator (already refitted on X_train by BayesSearchCV) and its params
    return opt.best_estimator_, opt.best_params_

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
    
    logger.info(f"Total samples: {X_full.shape[0]}")

    all_predictions = defaultdict(list)
    all_metrics = []
    
    all_best_params = [] 

    for repeat in range(10):
        logger.info(f"=== Repetition {repeat+1}/10 ===")
        skf_outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=88 + repeat)

        for fold, (train_idx, test_idx) in enumerate(skf_outer.split(X_full, y_full)):
            logger.info(f"  Fold {fold+1}/5")
            
            X_train, X_test = X_full[train_idx], X_full[test_idx]
            y_train, y_test = y_full[train_idx], y_full[test_idx]
            ids_test = ids_full[test_idx]
            
            fold_results = pd.DataFrame({"ID": ids_test, "Status": y_test})
            
            for model_name in get_search_spaces().keys():
                best_model, best_params = tune_and_train_model(X_train, y_train, model_name, n_iter=50)
                
                # Save params for record
                all_best_params.append({
                    "rep": repeat + 1,
                    "fold": fold + 1,
                    "model": model_name,
                    "params": dict(best_params)
                })

                probs = best_model.predict_proba(X_test)[:, 1]
                fold_results[model_name] = probs

                auc = roc_auc_score(y_test, probs)
                logloss_val = log_loss(y_test, probs)
                brier_val = brier_score_loss(y_test, probs)
                r_val, _ = pearsonr(probs, y_test)
                
                logger.info(f"    {model_name}: PearsonR={r_val:.3f}, AUC={auc:.3f}")
                
                all_metrics.append({
                    "rep": repeat + 1,
                    "fold": fold + 1,
                    "model": model_name,
                    "AUC": auc,
                    "LogLoss": logloss_val,
                    "Brier": brier_val,
                    "PearsonR": r_val
                })

                for _, row in fold_results.iterrows():
                    ID = row["ID"]
                    all_predictions[(ID, model_name)].append(row[model_name])

    # aggregation of all predictions
    final_records = []
    unique_ids = np.unique(ids_full)

    for ID in unique_ids:
        record = {"ID": ID}
        record["Status"] = df.loc[df["ID"] == ID, "Status"].iloc[0]
        for model_name in get_search_spaces().keys():
            preds = all_predictions.get((ID, model_name), [])
            if preds: 
                record[model_name] = np.mean(preds)
                record[f"{model_name}_SD"] = np.std(preds, ddof=1)
            else:
                record[model_name] = np.nan
                record[f"{model_name}_SD"] = np.nan
        final_records.append(record)

    prob_df = pd.DataFrame(final_records).sort_values("ID")
    
    # save outputs
    gebvs_dir = os.path.join("MLmodels", "gebvs", outdir)
    os.makedirs(gebvs_dir, exist_ok=True)
    
    prob_df.to_csv(os.path.join(gebvs_dir, f"{outdir}_GEBVs_NestedCV.csv"), index=False)
    logger.info(f"Saved GEBVs. Shape: {prob_df.shape}")

    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv(os.path.join(gebvs_dir, f"{outdir}_nested_metrics.csv"), index=False)
    logger.info("Saved nested CV metrics")

    params_path = os.path.join(gebvs_dir, f"{outdir}_all_nested_params.json")
    # Convert numpy types to python types for JSON serialization
    def convert(o):
        if isinstance(o, np.int64): return int(o)
        if isinstance(o, np.float64): return float(o)
        return o
        
    with open(params_path, "w") as f:
        json.dump(all_best_params, f, indent=4, default=convert)
    logger.info("Saved all hyperparameter sets from nested loops")

if __name__ == "__main__":
    main()