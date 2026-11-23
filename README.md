# teffects2
This repository contains companion code for the following paper:

Słoczyński, T., Uysal, S. D., & Wooldridge, J. M. (2025). Covariate Balancing and the Equivalence of Weighting and Doubly Robust Estimators of Average Treatment Effects

If you use this code please cite the paper using the bibtex reference below.

```
@article{słoczyński2025covariatebalancingequivalenceweighting,
  title={Covariate Balancing and the Equivalence of Weighting and Doubly Robust Estimators of Average Treatment Effects},
  author={Tymon Słoczyński, S. Derya Uysal and Jeffrey M. Wooldridge},
  url={https://arxiv.org/abs/2310.18563}, 
}
```

# Stata
It contains the Stata package teffects2. teffects2 estimates average treatment effects (ATEs) and average treatment effects on the treated (ATTs) using observational
    data.  As in Stata's official teffects command, inverse probability weighting (IPW), augmented inverse probability weighting
    (AIPW), and inverse probability weighted regression adjustment (IPWRA) estimators are supported.  However, unlike teffects,
    teffects2 supports covariate balancing estimation of the propensity score.

# R
It contains the R package to estimate AIPW, IPW, and IPWRA estimators using Covariate Balancing Methods: IPT, Graham et al. (2012), and CBPS, Imai and Ratkovic (2014) in addition to the standard approach using maximum likelihood. You can install it with the `devtools` (or `remotes`) package:

```r
devtools::install_github("deryauysal/teffects2", subdir = "R")
```

The `subdir = "R"` argument is needed because the R package is in the `R/` folder of this repository.

The package is based on the great work on treatment effect estimation in R by [Oliver Hines](https://github.com/ohines/teffectsR).

> **Status:** The R code is currently in **beta** and has **not been extensively tested**. Interfaces and results may change. Feedback, bug reports, and pull requests are very welcome.

Contact
-----------------
You can report issues and ask questions in the repository's issues page. If you choose to send an email instead, please direct it to Derya Uysal at derya.uysal@econ.lmu.de and include [teffects2] in the subject line.

## Acknowledgments
Special thanks to [**Lennart Ohly**](https://github.com/LennartOhly) for extensive work on the R implementation.  


Authors
-----------------
[Derya Uysal](https://www.deryauysal.com/)  
[Tymon Słoczyński](https://people.brandeis.edu/~tslocz/)  
[Jeffrey M. Wooldridge](https://econ.msu.edu/about/directory/Wooldridge-Jeffrey)  
