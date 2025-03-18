# TCRpred: incorporating T-cell receptor repertoire for clinical outcome prediction
[![License](https://img.shields.io/badge/license-LGPL--2.0-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.0.html)

## Overview
T-cell receptors (TCRs) play a crucial role in antigen recognition and adaptive immune response. High-throughput sequencing enables characterization of TCR repertoires at single-nucleotide resolution, offering valuable insights into the immune system and potential for clinical outcome prediction. However, integrating TCR data for prediction is challenging due to its unstructured and highly variable nature. We introduce **TCRpred**, a tool designed to incorporate both extracted and hidden TCR sequence features for clinical outcome prediction. Simulation studies show that TCRpred achieves strong predictive performance and outperforms alternative approaches. Application to real cancer datasets further demonstrates its practical utility in clinical research.


## Installation
You can install the TCRpred package by downloading the TCRpred_0.1.0.tar.gz file to your local system. Once downloaded, install it from your local file path using the following command in R:

```{r}
install.packages("path/to/TCRpred_0.1.0.tar.gz", repos = NULL, type = "source")
```
## Usage

```{r}
# Load the package
library(TCRpred)

data("TCRpred.data")                                                                                                                                                       
str(data)
 List of 5
  $ Y     : int [1:1000, 1] 0 0 1 0 0 0 0 1 1 0 ...
  $ X     : num [1:1000, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
   ..- attr(*, "dimnames")=List of 2
   .. ..$ : NULL
   .. ..$ : chr [1:2] "intercept" "X"
  $ K     : num [1:1000, 1:1000] 1.14 0.32 0.252 0.261 0.136 ...
   ..- attr(*, "dimnames")=List of 2
   .. ..$ : chr [1:1000] "TCGA-OR-A5JJ-01A-11R-A29S-07" "TCGA-OR-A5JB-01A-11R-A29S-07" "TCGA-PQ-A6FI-01A-11R-A31N-07" "TCGA-E7-A678-01A-11R-A30C-07" ...
   .. ..$ : chr [1:1000] "TCGA-OR-A5JJ-01A-11R-A29S-07" "TCGA-OR-A5JB-01A-11R-A29S-07" "TCGA-PQ-A6FI-01A-11R-A31N-07" "TCGA-E7-A678-01A-11R-A30C-07" ...
  $ Z     : num [1:1000, 1:5646] 0 0 0 0 0 0 0 0 0 0 ...
   ..- attr(*, "dimnames")=List of 2
   .. ..$ : chr [1:1000] "TCGA-OR-A5JJ-01A-11R-A29S-07" "TCGA-OR-A5JB-01A-11R-A29S-07" "TCGA-PQ-A6FI-01A-11R-A31N-07" "TCGA-E7-A678-01A-11R-A30C-07" ...
   .. ..$ : chr [1:5646] "AAA" "RAA" "NAA" "DAA" ...
  $ tcrdat:'data.frame': 24635 obs. of  3 variables:
   ..$ sid      : chr [1:24635] "TCGA-OR-A5JJ-01A-11R-A29S-07" "TCGA-OR-A5JJ-01A-11R-A29S-07" "TCGA-OR-A5JB-01A-11R-A29S-07" "TCGA-OR-A5JB-01A-11R-A29S-07" ...
   ..$ aaSeq    : chr [1:24635] "CASSLGSGDGYTF" "CASSPSWGYTF" "CASSYSGTDEQYF" "CASSLATGELFF" ...
   ..$ abundance: int [1:24635] 4 4 148 18 16 15 8 8 7 6 ...
```

`TCRpred.data` is a synthetic dataset designed to illustrate the usage of the **TCRpred** function. The dataset includes simulated TCR-related data, covariate, and outcome data for 1000 subjects.  

The dataset is structured as a **list containing five key components**. The **`Y` matrix** represents binary clinical outcomes, where each row corresponds to a subject, and the column indicates the outcome (0 or 1). The **`X` matrix** contains covariates, including an intercept term and an additional predictor variable, capturing clinical or demographic information for each subject. The **`K` matrix** is a similarity matrix that quantifies relationships between subjects based on T-cell receptor (TCR) repertoire features, with row and column names corresponding to subject IDs. The **`Z` matrix** encodes TCR sequence-derived features using k-mer representations, where each row represents a subject and each column corresponds to a specific k-mer feature. Additionally, the **`tcrdat` data frame** contains raw TCR sequence data, including subject identifiers (`sid`), amino acid sequences (`aaSeq`), and sequence abundance (`abundance`), providing essential information for further feature extraction and modeling.

```{r}
# Run TCRpred function with training data and specified parameters
out = TCRpred(Y = data$Y, X = data$X, K = data$K, Z = data$Z, ntrain = 500, maxiter = 10, tol = 0.01)

str(out)
 List of 2
  $ Y_true: int [1:500, 1] 0 0 0 0 0 1 1 0 0 1 ...
  $ Y_pred:Formal class 'dgeMatrix' [package "Matrix"] with 4 slots
   .. ..@ Dim     : int [1:2] 500 1
   .. ..@ Dimnames:List of 2
   .. .. ..$ : chr [1:500] "TCGA-OR-A5JJ-01A-11R-A29S-07" "TCGA-E7-A678-01A-11R-A30C-07" "TCGA-DK-A3WW-01A-22R-A23N-07" "TCGA-XF-A9T8-01A-11R-A39I-07" ...
   .. .. ..$ : chr "s0"
   .. ..@ x       : num [1:500] 0.2542 0.5479 0.0058 0.4412 0.578 ...
   .. ..@ factors : list()
```
## Interpretation of `TCRpred` Results

The **TCRpred** function was run using **500 training samples** with specified covariates, TCR features, and similarity matrices. The output is a **list containing two key components**:  

- **`Y_true`**: A numeric vector (`500 × 1`) representing the **true binary clinical outcomes** (0 or 1) for the training subjects.  
- **`Y_pred`**: A **sparse matrix (`dgeMatrix` from the `Matrix` package, size `500 × 1`)** containing the predicted probabilities for each subject.  
  - The `Dim` attribute confirms the dimensions (`500 × 1`).
  - The `Dimnames` attribute assigns **subject IDs** as row names (e.g., `"TCGA-OR-A5JJ-01A-11R-A29S-07"`).
  - The `x` slot contains the **predicted probability values** (ranging from 0 to 1), where **higher values** indicate a greater likelihood of the outcome being `1`.

### **Key Insights**
- The model **outputs probability scores** rather than binary classifications, allowing flexibility in choosing thresholds.
- A **subject with `Y_pred` > 0.5** is predicted as **likely belonging to the `1` class**, while **values < 0.5 suggest class `0`**.
- Sparse matrix storage (`dgeMatrix`) is used for efficiency in large datasets.


```{r}

# Evaluate prediction accuracy
eva_metric(out$Y_true, out$Y_pred)

     PPV        NPV clas_error        auc
0.6095238  0.7341772  0.2920000  0.7065732
```

## Evaluation of Prediction Accuracy

The `eva_metric` function was used to assess the predictive performance of **TCRpred** by comparing the predicted probabilities (`Y_pred`) to the true binary outcomes (`Y_true`). The results provide key classification performance metrics:

- **PPV (Positive Predictive Value)**: `0.6095`  
  - This represents **precision**, indicating that **60.95% of the predicted positive cases (class 1) are truly positive**.
  
- **NPV (Negative Predictive Value)**: `0.7342`  
  - This measures how well the model predicts negative cases, showing that **73.42% of the predicted negatives (class 0) are truly negative**.

- **Classification Error**: `0.2920`  
  - The overall error rate, meaning **29.2% of predictions were incorrect**.

- **AUC (Area Under the ROC Curve)**: `0.7066`  
  - The AUC score of **0.71** indicates **moderate predictive ability**, suggesting that the model performs better than random guessing but has room for improvement.
