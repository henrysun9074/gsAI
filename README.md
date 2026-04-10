# 🧬 Genomic Selection for Dermo Resistance in Eastern Oyster (*Crassostrea virginica*)

**Contact:** Henry Sun, hs325[at]duke.edu  
**DOI:** TBA  
**Cite:** TBA upon publication  

---

## 📖 Overview and Repository Structure

This project integrates data from three successive generations of lab-based dermo challenge and genotyping with a high-density SNP array for genomic selection and genome-wide association studies (GWAS) for dermo resistance in oysters, combining three generations into a single training population. We trained 9 different models to perform genomic selection and predict survival to dermo upon infection. We also evaluated the influence of rare variants, i.e. minor alleles, on genomic prediction accuracy, as well as strong-effect SNPs identified by GWAS. In the repository, please find the following folders.  

*/MLmodels* has code containing instructions for training, hyperparameter tuning, and cross-validation of LR, RF, GB genomic selection models, as well as .json files with optimal hyperparameter values for all tuned models.  
*/Rmodels* has code containing instructions for training and cross-validation of BayesB, BRR, LASSO, GBLUP, EGBLUP, RKHS genomic selection models.  
*/analysis* has code for statistical analyses comparing model performances and generating figures from the paper.  
**TBA**: Paul's code for GWAS

---

## ⚙️ Installation

All Python dependencies required for ML model training are listed in [`requirements.txt`](./requirements.txt).  
We recommend setting up a **Conda environment** for reproducibility.

```bash
# Clone the repository
git clone https://github.com/henrysun9074/gsAI.git
cd gsAI

# Create a new conda environment and install required packages with pip
conda create --name gsAI_env python=3.10
pip install -r requirements.txt

# Activate the environment
conda activate gsAI_env
