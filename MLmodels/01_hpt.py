import logging
import numpy as np
import pandas as pd
import joblib, os
import json
import argparse

from sklearn.metrics import make_scorer, brier_score_loss, roc_auc_score
from scipy.stats import pearsonr
from datetime import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler

from skopt import BayesSearchCV
from skopt.space import Real, Integer, Categorical

#  CLI Variables 
parser = argparse.ArgumentParser(description="Train model with selectable inputs")
parser.add_argument("--outdir", "-o", type=str, required=True,
                    help="Directory to save model outputs (will be created if missing)")
parser.add_argument("--filename", "-f", type=str, required=True,
                    help="CSV filename to load (path relative to current dir)")
parser.add_argument("--generation", "-g", type=str, default="all",
                    help='Generation filter: "F0", "F1", "F2", or "all" (default "all")')
parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
args = parser.parse_args()

outdir = args.outdir
filename = args.filename
generation = args.generation

logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                    format="%(asctime)s %(levelname)s:%(message)s")
logger = logging.getLogger(__name__)
logger.info("Running with outdir=%s filename=%s generation=%s", outdir, filename, generation)


#  Logging Setup 
logging.basicConfig(
    format="%(asctime)s  %(levelname)s  %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)

#  Model Search Spaces 
def get_search_spaces():
    return {
        "LR": (
            LogisticRegression(max_iter=1000, solver="saga"),
            {
                "C": Real(1e-5, 10, prior="log-uniform"),
                "penalty": Categorical(["l1", "l2"]),
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
                tree_method="hist", device="cpu", eval_metric="logloss"
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

#  Hyperparameter Tuning 

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

def tune_model(X, y, model_name, n_iter=100):
    base_model, search_space = get_search_spaces()[model_name]
    logger.info(f"Starting Bayesian optimization for {model_name}...")

    opt = BayesSearchCV(
        estimator=base_model,
        search_spaces=search_space,
        n_iter=n_iter,
        cv=5,
        scoring=pearson_scorer,
        n_jobs=-1,
        verbose=0,
    )
    opt.fit(X, y)
    logger.info(f"Best {model_name} params: {opt.best_params_}")
    return opt.best_estimator_, opt.best_params_

#  Main 
def main():
    chunksize = 100
    list_of_dataframes = []

    logger.info(f"Loading data from {filename}...")
    path_to_file = os.path.join("data", filename)
    for df in pd.read_csv(path_to_file, chunksize=chunksize):
        list_of_dataframes.append(df)
    df = pd.concat(list_of_dataframes)
    if generation != "all":
        df = df[df["Generation"] == generation]
    
    ax_columns = [col for col in df.columns if col.startswith("AX")]
    X = df[ax_columns]
    y = df["Status"]
    X = X.to_numpy()
    y = y.to_numpy()

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    tuned_models = {}
    tuned_params = {}

    for model_name in get_search_spaces().keys():
        best_model, best_params = tune_model(X, y, model_name)
        tuned_models[model_name] = best_model
        tuned_params[model_name] = best_params

    run_dir = os.path.join("MLmodels/models", outdir)
    os.makedirs(run_dir, exist_ok=True)

    #  Save hyperparameters 
    with open(os.path.join(run_dir, "best_hyperparams.json"), "w") as f:
        json.dump(tuned_params, f, indent=4)
    logger.info(f"Saved best hyperparameters to {run_dir}/best_hyperparams.json")

    #  Save full tuned models 
    for name, model in tuned_models.items():
        path = os.path.join(run_dir, f"{name}_best_model.joblib")
        joblib.dump(model, path)
        logger.info(f"Saved {name} model to {path}")


if __name__ == "__main__":
    main()
