#!/usr/bin/env python3
"""Nested hyperparameter tuning for genomic prediction CV using RF/GB only.

Execution:
    python3 nested_hpt_crossval.py -f data/genotypes.csv -o my_run -g all

Optionally, disable hyperparameter tuning using --no-tuning, or change the number of inner folds and search iterations
"""

import argparse
import importlib.metadata
import json
import logging
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


LOGGER = logging.getLogger(__name__)
N_REPEATS = 10
N_OUTER_FOLDS = 5
MODEL_NAMES = ("RF", "GB")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Nested hyperparameter tuning for genomic prediction CV using RF/GB only.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--filename", required=True, help="Input CSV path")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output name beneath MLmodels/gebvs and MLmodels/models")
    parser.add_argument("--model-outdir", "--indir", "-i", dest="model_outdir",
                        help="Optional separate model output name (no parameters are loaded)")
    parser.add_argument("-g", "--generation", default="all",
                        choices=["all", "F0", "F1", "F2"])
    parser.add_argument("--no-tuning", action="store_true",
                        help="Disable all hyperparameter searches; use fixed defaults")
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--search-iterations", type=int, default=5,
                        help="Total candidates per model per search, including initialization")
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--device", default="cpu", help="XGBoost device: cpu or cuda")
    parser.add_argument("--n-jobs", type=int, default=1,
                        help="Parallel inner-CV fits; use 1 for a single GPU")
    parser.add_argument("--model-jobs", type=int, default=1,
                        help="Threads per RF/GB fit; avoid oversubscribing CPUs")
    parser.add_argument("--models", nargs="+", choices=MODEL_NAMES, default=list(MODEL_NAMES),
                        help="Models to run; default retains all three original models")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    args.models = list(dict.fromkeys(args.models))
    if args.inner_folds < 2 or args.search_iterations < 1:
        parser.error("--inner-folds must be >= 2 and --search-iterations >= 1")
    if args.n_jobs == 0 or args.model_jobs == 0:
        parser.error("Job counts cannot be zero")
    if not 0 <= args.seed <= 2**32 - 1 - 10000:
        parser.error("--seed must be between 0 and 2**32 - 1 - 10000")
    # Keep both kinds of output inside their original project directories.
    for name in (args.outdir, args.model_outdir):
        if name is not None and (not name.strip() or Path(name).is_absolute()
                                 or ".." in Path(name).parts or name == "."):
            parser.error("Output names must be nonempty relative names without '..'")
    return args


def load_data(filename, generation):
    path = Path(filename).expanduser()
    if not path.is_file():
        path = Path("data") / filename
    if not path.is_file():
        raise ValueError(f"Input CSV not found: {filename} (also tried data/<filename>)")
    # String IDs preserve leading zeroes and avoid numeric ID coercion.
    df = pd.read_csv(path, dtype={"ID": str})
    required = {"ID", "Status"} | ({"Generation"} if generation != "all" else set())
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    if generation != "all":
        df = df.loc[df["Generation"] == generation].copy()
    df = df.reset_index(drop=True)
    if df.empty:
        raise ValueError("No animals remain after filtering")
    if df["ID"].isna().any() or df["ID"].duplicated().any():
        raise ValueError("ID must be nonmissing and unique: duplicate animals can leak across folds")
    features = [col for col in df.columns if col.startswith("AX")]
    if not features:
        raise ValueError("No AX-prefixed feature columns were found")
    X = df[features].apply(pd.to_numeric, errors="raise").astype(float)
    if not np.isfinite(X.to_numpy()).all():
        raise ValueError("AX features must be finite and nonmissing; no imputation is performed")
    y = pd.to_numeric(df["Status"], errors="raise").to_numpy()
    if not np.isin(y, [0, 1]).all() or len(np.unique(y)) != 2:
        raise ValueError("Status must contain both classes, encoded 0 and 1")
    return df, X, y.astype(int), path.resolve()


def build_pipeline(name, seed, args):
    if name == "RF":
        model = RandomForestClassifier(n_jobs=args.model_jobs, random_state=seed)
    elif name == "GB":
        from xgboost import XGBClassifier
        model = XGBClassifier(tree_method="hist", device=args.device,
                              eval_metric="logloss", n_jobs=args.model_jobs,
                              random_state=seed)
    else:
        raise ValueError(f"Unknown model: {name}")
    preprocessing = "passthrough"
    return Pipeline([("scaler", preprocessing), ("model", model)])


def get_search_space(name):
    from skopt.space import Categorical, Integer, Real
    spaces = {
        "RF": {"n_estimators": Integer(100, 2000), "max_depth": Integer(3, 50),
               "max_features": Categorical(["sqrt", "log2"]),
               "min_samples_split": Integer(2, 20), "min_samples_leaf": Integer(1, 10)},
        "GB": {"n_estimators": Integer(100, 2000), "max_depth": Integer(3, 15),
               "learning_rate": Real(1e-3, 0.3, prior="log-uniform"),
               "subsample": Real(0.5, 1.0), "colsample_bytree": Real(0.5, 1.0),
               "min_child_weight": Integer(1, 10), "gamma": Real(0, 5)},
    }
    return {f"model__{key}": value for key, value in spaces[name].items()}


def pearson_value(y, probs):
    if np.ptp(probs) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(pearsonr(probs, y).statistic)


def pearson_scorer(estimator, X, y):
    probs = estimator.predict_proba(X)[:, 1]
    if np.allclose(probs, probs[0]):
        return 0.0
    value = pearson_value(y, probs)
    return value if np.isfinite(value) else 0.0


def fit_model(name, X_train, y_train, seed, args):
    """Only training data enter this function, including during inner CV."""
    pipeline = build_pipeline(name, seed, args)
    if args.no_tuning:
        pipeline.fit(X_train, y_train)
        return pipeline, pipeline.named_steps["model"].get_params(deep=False), None

    from skopt import BayesSearchCV
    inner_cv = StratifiedKFold(n_splits=args.inner_folds, shuffle=True, random_state=seed)
    search = BayesSearchCV(
        estimator=pipeline,
        search_spaces=get_search_space(name),
        n_iter=args.search_iterations,
        # modify skopt's default which is 10 initialization points, which would make a five-candidate run entirely random
        optimizer_kwargs={"n_initial_points": min(2, max(1, args.search_iterations - 1))},
        n_points=1,
        cv=inner_cv,
        scoring=pearson_scorer,
        n_jobs=args.n_jobs,
        random_state=seed,
        refit=True,
        error_score="raise",
        return_train_score=False,
    )
    search.fit(X_train, y_train)
    params = {key.removeprefix("model__"): value for key, value in search.best_params_.items()}
    return search.best_estimator_, params, float(search.best_score_)


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_json(path, value):
    with path.open("w") as handle:
        json.dump(json_safe(value), handle, indent=4, allow_nan=False)
        handle.write("\n")


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s")
    # Fail early for dependency issues
    try:
        if "GB" in args.models:
            import xgboost  # noqa: F401
        if not args.no_tuning:
            import skopt  # noqa: F401
    except ImportError as exc:
        raise SystemExit(f"Missing dependency: {exc}. See installation command in the script.") from exc
    df, X, y, input_path = load_data(args.filename, args.generation)
    if np.bincount(y).min() < N_OUTER_FOLDS:
        raise ValueError("Each Status class needs at least 5 animals for outer CV")

    # Validate every split before fitting any models or creating outputs
    splits = []
    for repeat in range(N_REPEATS):
        outer = StratifiedKFold(n_splits=N_OUTER_FOLDS, shuffle=True,
                                random_state=args.seed + repeat)
        for fold, (train, test) in enumerate(outer.split(X, y)):
            if not args.no_tuning and np.bincount(y[train], minlength=2).min() < args.inner_folds:
                raise ValueError("An outer training fold has too few animals per class for inner CV")
            splits.append((repeat, fold, train, test))

    model_dir = Path("MLmodels/models") / (args.model_outdir or args.outdir)
    fold_dir = model_dir / "fold_models"
    fold_dir.mkdir(parents=True, exist_ok=True)
    gebv_path = Path("MLmodels/gebvs") / f"{args.outdir}_GEBVs.csv"
    metrics_path = Path("MLmodels/gebvs") / f"{args.outdir}_fold_metrics.csv"
    gebv_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "input_file": str(input_path), "arguments": vars(args),
        "outer_repeats": N_REPEATS, "outer_folds": N_OUTER_FOLDS,
        "features": list(X.columns), "positive_class": 1,
        "tuning_enabled": not args.no_tuning,
        "full_data_models_used_for_cv": False,
        "metric_row_order": "repeat, fold, model (in --models order)",
        "sample_ids_in_input_order": df["ID"].tolist(),
        "split_index_convention": "zero-based row positions after generation filtering",
        "outer_splits": [{"repeat": r + 1, "fold": f + 1,
                          "test_indices": test.tolist()} for r, f, _, test in splits],
        "versions": {package: importlib.metadata.version(package) for package in
                     ["numpy", "pandas", "scipy", "scikit-learn", "joblib"]
                     + (["xgboost"] if "GB" in args.models else [])
                     + ([] if args.no_tuning else ["scikit-optimize"])},
        "completed": False,
    }
    save_json(model_dir / "run_metadata.json", metadata)
    predictions = {name: np.full((N_REPEATS, len(y)), np.nan) for name in args.models}
    metrics = []
    fold_parameters = []
    LOGGER.info("Running 10 x 5-fold CV; tuning=%s, inner folds=%d, candidates=%d",
                not args.no_tuning, args.inner_folds, args.search_iterations)
    for repeat, fold, train, test in splits:
        seed = args.seed + repeat * N_OUTER_FOLDS + fold
        for name in args.models:
            LOGGER.info("Repeat %d/10, fold %d/5, model %s", repeat + 1, fold + 1, name)
            fitted, params, inner_score = fit_model(name, X.iloc[train], y[train], seed, args)
            probs = fitted.predict_proba(X.iloc[test])[:, 1]
            predictions[name][repeat, test] = probs
            metrics.append({"fold": fold + 1, "model": name,
                            "AUC": roc_auc_score(y[test], probs),
                            "LogLoss": log_loss(y[test], probs, labels=[0, 1]),
                            "Brier": brier_score_loss(y[test], probs),
                            "PearsonR": pearson_value(y[test], probs)})
            fold_metrics = metrics[-1]
            accuracy = accuracy_score(y[test], (probs > 0.5).astype(int))
            LOGGER.info(
                "Repeat %d/%d | Fold %d/%d | %s | "
                "Accuracy=%.4f | PearsonR=%.4f | AUC=%.4f | "
                "LogLoss=%.4f | Brier=%.4f",
                repeat + 1, N_REPEATS, fold + 1, N_OUTER_FOLDS, name,
                accuracy, fold_metrics["PearsonR"], fold_metrics["AUC"],
                fold_metrics["LogLoss"], fold_metrics["Brier"],
            )
            model_path = fold_dir / f"repeat_{repeat + 1:02d}_fold_{fold + 1:02d}_{name}.joblib"
            joblib.dump(fitted, model_path)
            fold_parameters.append({"repeat": repeat + 1, "fold": fold + 1,
                                    "model": name, "seed": seed, "best_params": params,
                                    "best_inner_PearsonR": inner_score,
                                    "model_file": str(model_path.relative_to(model_dir))})
        save_json(model_dir / "nested_cv_hyperparams.json", fold_parameters)

    result = pd.DataFrame({"ID": df["ID"], "Status": y})
    for name in args.models:
        if not np.isfinite(predictions[name]).all():
            raise RuntimeError(f"Missing/nonfinite held-out predictions for {name}")
        result[name] = predictions[name].mean(axis=0)
        result[f"{name}_SD"] = predictions[name].std(axis=0, ddof=1)
    result.sort_values("ID").to_csv(gebv_path, index=False)
    pd.DataFrame(metrics).to_csv(metrics_path, index=False)
    LOGGER.info("Saved held-out GEBVs to %s and fold metrics to %s", gebv_path, metrics_path)

    # Deployment only
    final_params = {}
    for name in args.models:
        LOGGER.info("Fitting full-data deployment model %s", name)
        fitted, params, _ = fit_model(name, X, y, args.seed + 10000, args)
        joblib.dump(fitted, model_dir / f"{name}_best_model.joblib")
        final_params[name] = params
    save_json(model_dir / "best_hyperparams.json", final_params)
    metadata["completed"] = True
    save_json(model_dir / "run_metadata.json", metadata)
    LOGGER.info("Complete: saved full-data models and parameters to %s", model_dir)


if __name__ == "__main__":
    main()
