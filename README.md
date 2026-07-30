# 🧬 Multigenerational machine learning-based genomic prediction for dermo resistance in eastern oyster *Crassostrea virginica* 
Henry Sun, Paul Coyne, Zhenwei Wang, Sandra Casas, Jerome La Peyre, Mason L. Williams, Scott Rikard, David Bushek, Juliet Wong, Ximing Guo

---

## 📖 Overview and Repository Structure

This project integrates data from three successive generations of lab-based dermo challenge and genotyping with a high-density SNP array for genomic selection and genome-wide association studies (GWAS) for dermo resistance in oysters, combining three generations into a single training population. We trained 9 different genomic prediction models to predict breeding values for survivial against dermo challenge. We also evaluated the influence of alleles with rare variants on genomic prediction accuracy, as well as strong-effect SNPs identified by GWAS. In the repository, please find the following folders.  

  * `*/MLmodels*` has code containing instructions for training, hyperparameter tuning, and cross-validation of LR, RF, GB genomic selection models, as well as .json files with optimal hyperparameter values for all tuned models.  
  * `*/Rmodels*` has code containing instructions for training and cross-validation of BayesB, BRR, LASSO, GBLUP, EGBLUP, RKHS genomic selection models.  
  * `*/analysis*` has code for statistical analyses comparing model performances and generating figures from the paper, as well as code for a companion genome-wide association study (Coyne et al., 2026).

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
```  

All R libraries for training genomic prediction models or generating figures can be installed through CRAN or Bioconductor.   

---

## Contact

Raw genotype and phenotype data used in this study are available upon request to xguo[at]hsrl.rutgers.edu. Please contact hs325[at]duke.edu with any other questions about this project or the contents of this repository.  

**Citation:** TBA upon study publication. 
