# 🧬 Genomic Selection for dermo Resistance in Eastern Oyster (*Crassostrea virginica*)

**Contact:** Henry Sun, hs325[at]duke.edu

---

## 📖 Overview

We test 9 different genomic selection models to predict survival to dermo upon infection in eastern oysters. Scripts here 

---

## 📂 Repository Structure

*/MLmodels* has code for training, hyperparameter tuning, and cross-validation of 3 different machine learning models.  
*/Rmodels* has code for training and cross-validation of 6 genomic selection models.  
*/analysis* has code for generating all visualizations associated with the project. 

---

## ⚙️ Installation

All Python dependencies are listed in [`requirements.txt`](./requirements.txt).  
We recommend setting up a **Conda environment** for reproducibility.

```bash
# Clone the repository
git clone https://github.com/henrysun9074/gsAI.git
cd gsAI

# Create a new conda environment and install with pip
conda create --name gsAI_env python=3.10
pip install -r requirements.txt

# Activate the environment
conda activate gsAI_env

# Install necessary packages
pip install -r requirements.txt
