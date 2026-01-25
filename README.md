# bartXViz: Visualization of BART and BARP using SHAP

### The **bartXViz** package provides SHAP-based model explanation tools for Bayesian Additive Regression Trees (BART) and Bayesian Additive Regression Trees with Post-Stratification (BARP).

The version uploaded on January 25, 2026 corresponds to **v1.0.10**, the same version currently released on CRAN.
([https://cran.r-project.org/web/packages/bartXViz/index.html](https://cran.r-project.org/web/packages/bartXViz/index.html))

---

Complex machine learning models are often difficult to interpret.  
**Shapley values** provide a principled framework for understanding why a model makes specific predictions by quantifying each variable's contribution.  

This package implements **permutation-based Shapley values** for BART and BARP models, enabling users to evaluate variable importance and contribution across Bayesian posterior samples obtained through MCMC.  
The SHAP approach follows the method proposed by **Strumbel and Kononenko (2014)** <doi:10.1007/s10115-013-0679-x>, adapted to the Bayesian tree ensemble framework introduced by **Chipman, George, and McCulloch (2010)** <doi:10.1214/09-AOAS285>.  

**bartXViz** is compatible with several popular R implementations of BART, including:  
- **BART** <doi:10.18637/jss.v097.i01>  
- **bartMachine** <doi:10.18637/jss.v070.i04>  
- **dbarts** (<https://CRAN.R-project.org/package=dbarts>)  

For gradient boosting and baseline comparisons, the package also considers the SHAP framework proposed by **Lundberg et al. (2020)** <doi:10.1038/s42256-019-0138-9>.  

The **BARP** model, originally proposed by **Bisbee (2019)** <doi:10.1017/S0003055419000480>, is implemented with reference to [jbisbee1/BARP](https://github.com/jbisbee1/BARP).  
BARP extends post-stratification to compute variable contributions within each stratum defined by stratifying variables, improving small-area estimation interpretability.  

The resulting Shapley values can be visualized through both **global** and **local** explanation methods, allowing users to explore model interpretability in Bayesian tree ensembles with intuitive visualizations.

---
