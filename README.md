# GSCA.Reg_Prime

### Author:
Gyeongcheol Cho

## Description:
- The **GSCA.Reg_Prime** package enables users to estimate and evaluate Regularized GSCA models.

## Features:
- Estimate parameters of Regularized GSCA models and calculate their standard errors (SE) along with 95% confidence intervals (CI).
- Assess model performance based on both explanatory and predictive power.
- Compute the PET (Predictor Exclusion Threshold) statistic to evaluate the predictive power of individual predictor components.
- Enable parallel computing for bootstrap sampling.

## Installation:
To use this package in MATLAB:
1. Clone or download the repository:
   ```bash
   git clone https://github.com/PsycheMatrica/GSCA.Reg_Prime.git
   ```
2. Add the package to your MATLAB path:
   ```matlab
    addpath(genpath('GSCA.Reg_Prime'))
   ```

## Usage:
- For examples on how to use the package, refer to the `Run_Example_RegGSCA.m` file. This file demonstrates the implementation of `RegGSCA()` using the ACSI dataset.

## Compatibility:
- Tested on MATLAB R2023b.
- Likely compatible with earlier MATLAB versions.

## Future Update
- Lasso regularization will be accomodated in the near future.

### Citation (APA):
- If you use **GSCA.Reg_Prime** in your research or publications, please cite it in APA format as follows:

```plaintext
Cho, G. (2024). GSCA.Reg_Prime: A package for regularized generalized structured component analysis [Computer software]. GitHub. https://github.com/PsycheMatrica/GSCA.Reg_Prime.git
```
