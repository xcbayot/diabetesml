# diabetesML  
**R and QIIME2 Scripts for Hypervariable sOTU Machine Learning Analysis**  

This repository contains scripts for analyzing 16S metabarcoding data to predict diabetes status using machine learning models. The analysis includes subsampling, normalization, and model training across different sOTU subsets (100, 500, 1000 hyper-variable sOTUs).  

---

##  **Pipeline Overview**  
1. **QIIME2 Analysis:**  
   - Demultiplexing and OTU calculation from 16S metabarcoding FASTQ files.  
2. **R Analysis:**  
   - Subsampling and normalization of hyper-variable sOTUs.  
   - Machine learning model training and evaluation for diabetes prediction.  

---

##  **Dependencies**  
### **QIIME2 Environment:**  
- QIIME2 version **2022.11**    
   - `dplyr` 
### **R Environment:**  
- **R version 4.3.3**  
- **Required R Packages:**  
   - `dplyr`  
   - `ggplot2`  
   - `vegan`  
   - `phyloseq`  
   - `qiime2R`  
   - `caret`  
   - `randomForest`  
   - `rpart`  
   - `glmnet`  
   - `kernlab`  
   - `pROC`
 
   - To install required R packages:  
```r  
install.packages(c("dplyr", "ggplot2", "vegan", "phyloseq", "caret", "randomForest", "rpart", "glmnet", "kernlab", "pROC"))  
