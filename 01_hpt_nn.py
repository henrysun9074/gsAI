## TODO: Write function for optimizing R-value. Then run 100-iteration BayesOpt for HP tuning for MLP and CNN.
## Draft hidden layer architecture

import logging
import numpy as np
import pandas as pd
import joblib, os
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import seaborn as sns
import matplotlib.pyplot as plt

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
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix
from scipy.stats import pearsonr

from xgboost import XGBClassifier
from bayes_opt import BayesianOptimization

from joblib import Parallel, delayed

logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s", level=logging.INFO)
logger = logging.getLogger(__name__)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)
print("PyTorch path:", torch.__path__ )

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

ids = df["ID"].values  
X = df[ax_columns].values
y = df["Status"].values

X_train, X_test, y_train, y_test, id_train, id_test = train_test_split(
    X, y, ids, test_size=0.2, random_state=42
)

scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# Convert to tensors
X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
y_train_tensor = torch.tensor(y_train, dtype=torch.float32).unsqueeze(1)
X_test_tensor = torch.tensor(X_test, dtype=torch.float32)
y_test_tensor = torch.tensor(y_test, dtype=torch.float32).unsqueeze(1)

train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)

class SimpleMLP(nn.Module):
    def __init__(self, input_dim):
        super(SimpleMLP, self).__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 1),   
            nn.Sigmoid()       
        )
    def forward(self, x):
        return self.network(x)

model = SimpleMLP(input_dim=X_train.shape[1]).to(device)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

n_epochs = 20

for epoch in range(n_epochs):
    model.train()
    epoch_loss = 0
    for X_batch, y_batch in train_loader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)
        optimizer.zero_grad()
        y_pred = model(X_batch)
        loss = criterion(y_pred, y_batch)
        loss.backward()
        optimizer.step()
        epoch_loss += loss.item()
    
    print(f"Epoch {epoch+1}/{n_epochs}, Loss: {epoch_loss/len(train_loader):.4f}")

model.eval()
with torch.no_grad():
    probs = model(X_test_tensor.to(device))
    probs = probs.squeeze().cpu().numpy()   
    preds = (probs > 0.5).astype(int)  

acc = accuracy_score(y_test, preds)
auc = roc_auc_score(y_test, probs)
r_val, _ = pearsonr(probs, y_test)

print(f"Accuracy: {acc:.3f}, AUC: {auc:.3f}, R-value: {r_val:.3f}")
results_df = pd.DataFrame({
    "ID": id_test,
    "TrueLabel": y_test,
    "PredProb": probs,
    "PredClass": preds
})
results_df
