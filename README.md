# 🧬 Genomic Selection for Dermo Resistance in Eastern Oyster (*Crassostrea virginica*)

**Author:** Henry Sun, hs325[at]duke.edu

---

## 📖 Overview

This project applies **genomic selection (GS)** and **machine learning (ML)** approaches to predict **genomic estimated breeding values (GEBVs)** for **resistance to Dermo disease** in the Eastern oyster (*Crassostrea virginica*).
We analyze three generations (~2,400 oysters) genotyped at approximately **66,000 SNP loci** using a high-density SNP array. 
The goal is to train, optimize, and evaluate predictive models that support selective breeding for disease resistance.

---

## 📂 Repository Structure

TBA upon project completion

---

## ⚙️ Installation

All Python dependencies are listed in [`requirements.txt`](./requirements.txt).  
We recommend setting up a **Conda environment** for reproducibility.

### 1. Install Conda

You can install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution).

---

### 2. Create and Activate the Environment

```bash
# Clone the repository
git clone https://github.com/henrysun9074/gsAI.git
cd gsAI

# Create a new conda environment
conda create --name gsAI_env python=3.10

# Activate the environment
conda activate gsAI_env

# Install necessary packages
pip install -r requirements.txt
