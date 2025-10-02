from datetime import datetime
def print_with_time(message):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)

print_with_time("Starting Python script and importing packages...")
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from scipy.stats import pearsonr
from sklearn.frozen import FrozenEstimator
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix
from sklearn.metrics import accuracy_score, roc_auc_score
import seaborn as sns
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import make_scorer
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
print_with_time("Packages loaded")

# Use GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

###############################################################################
print_with_time("Loading and concatenating data chunks...")
chunksize = 100
list_of_dataframes = []

for df in pd.read_csv('MeanImputedScaledData.csv', chunksize=chunksize, index_col=0):
    list_of_dataframes.append(df)

result = pd.concat(list_of_dataframes)
df = result
df['ID'] = df.index
ids = df["ID"].to_numpy()

print_with_time("Preparing features and labels...")
ax_columns = [col for col in df.columns if col.startswith('AX')]
# len(ax_columns)
X = df[ax_columns]
y = df["Status"]
X = X.to_numpy()
y = y.to_numpy()

scaler = StandardScaler()
X = scaler.fit_transform(X)

###############################################################################

fold_accuracies = []
gebv_records = []

# 50 total folds. Ensures representative use of dataset + 10 predictions for each oyster

for loop_idx in range(10):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=loop_idx)

    print(f"\nStarting fold {loop_idx} out of 10")

    for fold_idx, (train_val_idx, test_idx) in enumerate(skf.split(X, y)):
        X_train_val, X_test = X[train_val_idx], X[test_idx]
        y_train_val, y_test = y[train_val_idx], y[test_idx]
        ids_train_val, ids_test = ids[train_val_idx], ids[test_idx]

        X_train, X_calib, y_train, y_calib, ids_train, ids_calib = train_test_split(
            X_train_val, y_train_val, ids_train_val, test_size=0.2,
            stratify=y_train_val, random_state=loop_idx
        )

        gb_model = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)
        gb_model.fit(X_train, y_train)

        gb_frozen = FrozenEstimator(gb_model)
        calibrated_gb = CalibratedClassifierCV(estimator=gb_frozen, method='sigmoid')
        calibrated_gb.fit(X_calib, y_calib)

        y_probs = calibrated_gb.predict_proba(X_test)
        y_preds = np.argmax(y_probs, axis=1)
        acc = accuracy_score(y_test, y_preds)
        fold_accuracies.append(acc)

        # Record GEBVs with ID
        for id_val, true_label, prob_class1 in zip(ids_test, y_test, y_probs[:, 1]):
            gebv_records.append({
                'loop': loop_idx,
                'fold': fold_idx,
                'ID': id_val,
                'true_label': true_label,
                'gebv': prob_class1
            })

        print_with_time(f"    Accuracy for fold {fold_idx + 1}: {acc:.4f}")

print("\nAverage Accuracy over 50 folds:", np.mean(fold_accuracies))
print("Standard Deviation:", np.std(fold_accuracies))

gebv_df = pd.DataFrame(gebv_records)
gebv_df = gebv_df.sort_values('ID')
gebv_df.to_csv("gebv_df_GB.csv", index=False)
print_with_time("Script complete.")

