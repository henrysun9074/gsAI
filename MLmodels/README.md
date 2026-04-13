**Directory Structure**  
01_hpt.py: Script for hyperparameter tuning and saving of model weights  
02_crossval.py: Script for 10 iterations of five-fold cross-validation, loads hyperparameter values  

The remaining directories have trained hyperparameter values stored in *.JSON* files:  
"MAF"[X] denotes model trained on minor allele frequency X  
"Old" denotes training without GSMs, "Extra" denotes training with GSMs  
"Default" contains default sklearn hyperparameter values loaded for models trained at all MAFs without GSMs  
