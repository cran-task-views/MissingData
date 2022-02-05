---
name: MissingData
topic: Missing Data
maintainer: Julie Josse, Imke Mayer, Nicholas Tierney, and Nathalie Vialaneix (r-miss-tastic team)
email: r-miss-tastic@clementine.wf
version: 2021-12-29
source: https://github.com/cran-task-views/MissingData/
---

Missing data are very frequently found in datasets. Base R provides a few 
options to handle them using computations that involve only observed data 
(`na.rm = TRUE` in functions `mean`, `var`, ... or 
`use = complete.obs|na.or.complete|pairwise.complete.obs` in functions `cov`, 
`cor`, ...). The base package `stats` also contains the generic function 
`na.action` that extracts information of the `NA` action used to create an 
object.

These basic options are complemented by many packages on CRAN. In this task 
view, we focused on the most important ones, which have been published more than
one year ago and are regularly updated. The task view is structured into main 
topics:

-   [Exploration of missing data](#exploration)
-   [Likelihood based approaches](#likelihood)
-   [Single imputation](#single)
-   [Multiple imputation](#multiple)
-   [Weighting methods](#weights)
-   [Specific types of data](#data)
-   [Specific tasks](#tasks)
-   [Specific application fields](#applications)

In addition to the present task view, this [reference website on missing
data](https://rmisstastic.netlify.com/) might also be helpful. Complementary
information might also be found in `r view("TimeSeries")`,
`r view("SpatioTemporal")`, and `r view("OfficialStatistics")`.

If you think that we missed some important packages in this list, please
e-mail the maintainers or submit an issue or pull request in the GitHub
repository linked above.

[**Exploration of missing data**]{#exploration}

-   *Manipulation of missing data* is implemented in the packages
    `r pkg("sjmisc")` and `r pkg("sjlabelled")`.
    `r pkg("memisc")` also provides defineable missing
    values, along with infrastructure for the management of survey data
    and variable labels.
-   *Missing data patterns* can be identified and explored using the
    packages `r pkg("mi")`, `r pkg("wrangle")`,
    `r pkg("DescTools")`, `r pkg("dlookr")` and
    `r pkg("naniar", priority = "core")`.
-   *Graphics that describe distributions and patterns of missing data*
    are implemented in `r pkg("VIM", priority = "core")` (which
    has a Graphical User Interface, VIMGUI, currently archived on CRAN)
    and `r pkg("naniar")` (which abides by
    [tidyverse](https://tidyverse.org) principles).
-   *Tests of the MAR assumption (versus the MCAR assumption):*
    `r pkg("RBtest")` proposes a regression based approach to
    test the missing data mechanism and `r pkg("samon")`
    performs sensitivity analysis in clinical trials to check the
    relevance of the MAR assumption.
-   *Evaluation* of the quality of imputation can be performed using the
    function `ampute` of `r pkg("mice", priority = "core")` through simulations,
    with the `r pkg("Iscores")` with a KL-based scoring rule, or with the 
    benchmark available in `r pkg("imputeTestbench")` (for univariate time
    series imputation). In addition, `r pkg("mi")` and `r pkg("VIM")` also 
    provide diagnostic plots to evaluate the quality of imputation.

[**Likelihood based approaches**]{#likelihood}

-   *Methods based on the Expectation Maximization (EM) algorithm* are
    implemented in `r pkg("norm")` (using the function
    `em.norm` for multivariate Gaussian data),
    `r pkg("norm2")` (using the function `emNorm`), in
    `r pkg("cat")` (function `em.cat` for multivariate
    categorical data), in `r pkg("mix")` (function `em.mix`
    for multivariate mixed categorical and continuous data). These
    packages also implement *Bayesian approaches* (with Imputation and
    Posterior steps) for the same models (functions `da.`XXX for
    `norm`, `cat` and `mix`) and can be used to obtain imputed complete
    datasets or multiple imputations (functions `imp.`XXX for `norm`,
    `cat` and `mix`), once the model parameters have been estimated.
    `r pkg("imputeR")` is a Multivariate
    Expectation-Maximization (EM) based imputation framework that offers
    several different algorithms, including Lasso, tree-based models or
    PCA. In addition, `r pkg("TestDataImputation")` implements
    imputation based on EM estimation (and other simpler imputation
    methods) that are well suited for dichotomous and polytomous tests
    with item responses.
-   *Full Information Maximum Likelihood* (also known as "direct
    maximum likelihood" or "raw maximum likelihood") is available in
    `r pkg("lavaan")` (and in its extension
    `r pkg("semTools")`), `r pkg("OpenMx")` and
    `r pkg("rsem")`, for handling missing data in structural
    equation modeling.
-   *Bayesian approaches* for handling missing values in model based
    clustering with variable selection is available in
    `r pkg("VarSelLCM")`. The package also provides imputation
    using the posterior mean.
-   *Missing values in generalized linear models* is provided in package 
    `r pkg("mdmb")` for various families. `r pkg("JointAI")` implements 
    Bayesian approaches for generalized linear mixed models.
-   *Missing data in item response models* is implemented in
    `r pkg("TAM")`, `r pkg("mirt")` and
    `r pkg("ltm")`.

[**Single imputation**]{#single}

-   The simplest method for missing data imputation is *imputation by
    mean* (or median, mode, \...). This approach is available in many
    packages among which `r pkg("ForImp")` and
    `r pkg("Hmisc")` that contain various proposals for
    imputing with the same value all missing instances of a variable.
-   *Generic packages*: The packages `r pkg("VIM")` and `r pkg("filling")` 
    contain several popular methods for missing value imputation (including some
    listed in the sections dedicated to specific methods as listed below). 
    In addition, `r pkg("simputation")` is a general package for imputation
    by any prediction method that can be combined with various regression 
    methods, and works well with the tidyverse.
-   *k-nearest neighbors* is a popular method for missing data imputation that 
    is available in many packages including the main packages
    `r pkg("yaImpute", priority = "core")` (with many different methods for 
    kNN imputation, including a CCA based imputation) and `r pkg("VIM")`. It
    is also available in `r bioc("impute")` (where it is oriented toward 
    microarray imputation).\
    `r pkg("isotree")` uses a similar approach to impute missing value, which is
    based on similarities between samples and isolation forests.
-   *hot-deck* imputation is implemented in the package 
    `r pkg("hot.deck", priority = "core")`, with various possible settings 
    (including multiple imputation). It is also available in `r pkg("VIM")` 
    (function `hotdeck`) and a fractional version (using weights) is provided
    in `r pkg("FHDI")`. `r pkg("StatMatch")` also uses hot-deck imputation
    to impute surveys from an external dataset.\
    Similarly, `r pkg("impimp")` uses the notion of "donor" to impute a set of 
    possible values, termed "imprecise imputation".
-   Imputation *based on random forest* is implemented in `r pkg("missForest")` 
    with a faster version in `r pkg("missRanger")`.
-   *Other regression based imputations* are implemented in `r pkg("VIM")` 
    (linear regression based imputation in the function `regressionImp`). 
    `r pkg("iai")` tunes optimal imputation based on knn, tree or SVM.
-   *Matrix completion* is implemented with iterative PCA/*SVD*-decomposition in
    the package `r pkg("missMDA", priority = "core")` for numerical, categorical 
    and mixed data (including imputation of groups). *NIPALS* (also based on SVD
    computation) is implemented in the packages `r bioc("mixOmics")` (for PCA 
    and PLS), `r pkg("ade4")`, `r pkg("nipals")` and `r pkg("plsRglm")` (for 
    generalized model PLS). `r pkg("cmfrec")` is also a large package dedicated
    to matrix factorization (for recommender systems), which includes
    imputation. Other PCA/factor based imputations are available in 
    `r bioc("pcaMethods")` (with a Bayesian implementation of PCA), in 
    `r pkg("primePCA")` (for heterogeneous missingness in high-dimensional PCA)
    and `r pkg("tensorBF")` (for 3-way tensor data).\
    *Low rank based imputation* is provided in 
    `r pkg("softImpute", priority = "core")`, which contains several
    methods for iterative matrix completion. This method is also available in
    the very general package `r pkg("rsparse")`, which contains various tools
    for sparse matrices. Variants based on low rank assumptions are available
    in `r pkg("denoiseR")`, in `r pkg("mimi")`, in `r pkg("ECLRMC")` (for
    ensemble matrix completion), and in `r pkg("ROptSpace")` (with a 
    computationally efficient approach).
-   Imputation *based on copula* is implemented in `r pkg("CoImp")` with a
    semi-parametric imputation procedure and in `r pkg("mdgc")` using Gaussian
    copula for mixed data types.

[**Multiple imputation**]{#multiple}

Some of the above mentioned packages can also handle multiple
imputations.

-   `r pkg("Amelia", priority = "core")` implements Bootstrap
    multiple imputation using EM to estimate the parameters, for
    quantitative data it imputes assuming a Multivariate Gaussian
    distribution. In addition, AmeliaView is a GUI for
    `r pkg("Amelia")`, available from the 
    [Amelia web page](https://gking.harvard.edu/amelia).
    `r pkg("NPBayesImputeCat")` also implements multiple
    imputation by joint modelling for categorical variables with a
    Bayesian approach.
-   `r pkg("mi")`, `r pkg("mice")` and
    `r pkg("smcfcs")` implement multiple imputation by Chained
    Equations. `r pkg("smcfcs")` extends the models covered by
    the two previous packages. `r pkg("miceFast")` provides an
    alternative implementation of mice imputation methods using object
    oriented style programming and C++. `r pkg("bootImpute")`
    performs bootstrap based imputations and analyses of these
    imputations to use with `r pkg("mice")` or
    `r pkg("smcfcs")`. `r pkg("miceRanger")`
    performs multiple imputation by chained equations using random
    forests.
-   `r pkg("missMDA")` implements multiple imputation based on
    SVD methods.
-   `r pkg("hot.deck")` implements hot-deck based multiple
    imputation.
-   *Multilevel imputation*: Multilevel multiple imputation is
    implemented in `r pkg("hmi")`,
    `r pkg("jomo", priority = "core")`,
    `r pkg("mice")`, `r pkg("miceadds")`,
    `r pkg("micemd")`, `r pkg("mitml")`, and
    `r pkg("pan")`.
-   `r pkg("Qtools")` and `r pkg("miWQS")` implement
    multiple imputation based on quantile regression.
-   `r pkg("lodi")` implements the imputation of observed
    values below the limit of detection (LOD) via censored likelihood
    multiple imputation (CLMI).
-   `r pkg("BaBooN")` implements a Bayesian bootstrap approach
    for discrete data imputation that is based on Predictive Mean
    Matching (PMM).

In addition, `r pkg("mitools")` provides a generic approach to
handle multiple imputation in combination with any imputation method.
And `r pkg("NADIA")` provides a uniform interface to compare
the performances of several imputation algorithms.

[**Weighting methods**]{#weights}

-   *Computation of weights* for observed data to account for unobserved
    data by *Inverse Probability Weighting (IPW)* is implemented in
    `r pkg("ipw")`. IPW is also used for quantile estimations and boxplots in 
    `r pkg("IPWboxplot")`.
-   *Doubly Robust Inverse Probability Weighted Augmented GEE Estimator
    with missing outcome* is implemented in `r pkg("CRTgeeDR")`.

[**Specific types of data**]{#data}

-   *Longitudinal data / time series and censored data*: Imputation for time 
    series is implemented in `r pkg("imputeTS", priority = "core")`. Other 
    packages, such as `r pkg("forecast")`, `r pkg("spacetime")`, 
    `r pkg("timeSeries")`, `r pkg("xts")`, `r pkg("prophet")`, 
    `r pkg("stlplus")` or `r pkg("zoo")`, are dedicated to time series but also 
    contain some (often basic) methods to handle missing data (see also 
    `r view("TimeSeries")`). Based on tidy principle, the `r pkg("padr")` and
    `r pkg("tsibble")` also provide methods for imputing missing values in
    time series.\
    More specific methods are implemented in other packages: imputation of time 
    series based on Dynamic Time Warping is implemented in the family of 
    packages `r pkg("DTWBI")`, `r pkg("DTWUMI")`, and `r pkg("FSMUMI")` for 
    univariate and multivariate time series. `r pkg("BMTAR")` provides an
    an estimation of the autoregressive threshold models with Gaussian noise 
    using a Bayesian approach in the presence of missing data in multivariate 
    time series. `r pkg("swgee")` implements an IPW approach for longitudinal 
    data with missing observations.
-   *Spatial data*: Imputation for spatial data is implemented in the package
    `r pkg("rtop")`, which performs geostatistical interpolation of irregular
    areal data, and in `r pkg("areal")`, which performs areal weighted 
    interpolation using a tidyverse data management.\
    Interpolation of spatial data based on genetic distances is also 
    available in `r pkg("phylin")`. 
-   *Spatio-temporal data* (see also `r view("SpatioTemporal")`): Imputation 
    for spatio-temporal data is implemented in the packages `r pkg("cutoffR")` 
    (using different methods as knn and SVD) and in in `r pkg("StempCens")` 
    with a SAEM approach that approximates EM when the E-step does not have an 
    analytic form.\
    `r pkg("gapfill")` is dedicated to imputation of satellite data observed at
    equally-spaced points in time and `r pkg("momentuHMM")` to the analysis of
    telemetry data using generalized hidden Markov models (including multiple 
    imputation for missing data). 
-   *Graphs/networks*: `r pkg("missSBM")` imputes missing edges in Stochastic
    Block models and `r pkg("cglasso")` implements an extension of the
    Graphical Lasso inference from censored and missing value measurements.
-   *Imputation for contingency tables* is implemented in
    `r pkg("lori")` that can also be used for the analysis of
    contingency tables with missing data.
-   *Imputation for compositional data (CODA)* is implemented in in
    `r pkg("zCompositions")` (various imputation methods for
    zeros, left-censored and missing data).
-   *Imputation for meta-analyses* of binary outcomes is provided in
    `r pkg("metasens")`.
-   *Experimental design*: `r pkg("experiment")` handles missing values in
    experimental design such as randomized experiments with missing covariate 
    and outcome data, matched-pairs design with missing outcome.
-   *Recurrent events*: `r pkg("dejaVu")` performs multiple imputation of
    recurrent event data based on a negative binomial regression model.
    
[**Specific tasks**]{#tasks}

-   *Regression and classification*: many different supervised methods can
    accommodate the presence of missing values. `r pkg("randomForest")`, 
    `r pkg("grf")`, and `r pkg("StratifiedRF")` handle missing values in 
    predictors in various random forest based methods. `r pkg("misaem")` 
    handles missing data in linear and logistic regression and allows for model 
    selection. `r pkg("psfmi")` also provides a framework for model selection
    for various linear models in multiply imputed datasets. 
    `r pkg("naivebayes")` provides an efficient implementation of the naive 
    Bayes classifier in the presence of missing data. `r pkg("plsRbeta")` 
    implements PLS for beta regression models with missing data in the 
    predictors. `r pkg("lqr")` provides quantile regression estimates based on 
    various distributions in the presence of missing values and censored data.
    `r pkg("eigenmodel")` handles missing values in regression models for 
    symmetric relational data.
-   *Clustering*: `r pkg("biclustermd")` handles missing data
    in biclustering. `r pkg("RMixtComp")`,
    `r pkg("MGMM")` and `r pkg("mixture")` fit
    various mixture models in the presence of missing data.
-   *Tests* for two-sample paired missing data are implemented in
    `r pkg("robustrank")`.
-   *Outlier detection* (and robust analysis) in the presence of missing values
    is implemented in `r pkg("GSE")` and `r pkg("rrcovNA")`

[**Specific application fields**]{#applications}

-   *Genetics*: Imputation of SNP data is implemented in `r pkg("alleHap")`
    (using solely deterministic techniques based on pedigree data), in
    `r pkg("QTLRel")` (using information on flanking SNPs), in 
    `r bioc("snpStats")` (using a nearest neighbor approach), in 
    `r pkg("HardyWeinberg")` (using multiple imputations with a multinomial
    model based on allele intensities and/or flanking SNPs).\
    EM algorithm is used to compute genetic statistics for population in the
    presence of missing SNP in `r pkg("StAMPP")`.
-   *Genomics*: Imputation for dropout events (*i.e.* , under-sampling of mRNA 
    molecules) in single-cell RNA-Sequencing data is implemented
    in `r pkg("DrImpute")`, `r pkg("Rmagic")`, and `r pkg("SAVER")` and ,
    based, respectively, on clustering of cells, Markov affinity graph, and an
    empirical Bayes approach. These packages are used and combined in
    `r bioc("scRecover")` and `r bioc("ADImpute")`\
    `r pkg("RNAseqNet")` uses hot-deck imputation to improve RNA-seq network 
    inference with an auxiliary dataset.
-   *Epidemiology*: `r pkg("idem")` implements a procedure for comparing
    treatments in clinical trials with missed visits or premature withdrawal. 
    `r pkg("InformativeCensoring")` implements multiple imputation for 
    informative censoring. `r pkg("pseval")` evaluates principal surrogates in
    a single clinical trial in the presence of missing counterfactual surrogate 
    responses. `r pkg("sievePH")` implements continuous, possibly multivariate, 
    mark-specific hazard ratio with missing values in multivariate marks using 
    an IPW approach. `r pkg("icenReg")` performs imputation for censored
    responses for interval data. 
-   *Proteomics*: Imputation of missing data in LC-MS/MS spectra for protein
    quantification is available in `r pkg("aLFQ")`.
-   *Health*: `r pkg("missingHE")` implements models for health economic 
    evaluations with missing outcome data. `r pkg("accelmissing")` provides 
    multiple imputation with the zero-inflated Poisson lognormal model for
    missing count values in accelerometer data.
-   *Environement*: `r pkg("AeRobiology")` imputes missing data in 
    aerobiological datasets imported from aerobiological public databases.
-   *Causal inference*: Causal inference with interactive fixed-effect
    models is available in `r pkg("gsynth")` with missing
    values handled by matrix completion. `r pkg("MatchThem")`
    matches multiply imputed datasets using several matching methods,
    and provides users with the tools to estimate causal effects in each
    imputed datasets. `r pkg("grf")` offers treatment effect
    estimation with incomplete confounders and covariates under modified
    unconfoundedness assumptions.
-   *Finance*: `r pkg("imputeFin")` handles imputation of missing values in 
    financial time series using AR models or random walk.
-   *Scoring*: Basic methods (mean, median, mode, \...) for imputing
    missing data in scoring datasets are proposed in
    `r pkg("scorecardModelUtils")`.
-   *Preference models*: Missing data in preference models are handled
    with a composite link approach that allows for MCAR and MNAR patterns to be 
    accounted for in `r pkg("prefmod")`.
-   *Administrative records / Surveys*: `r pkg("fastLink")`
    provides a Fellegi-Sunter probabilistic record linkage that allows
    for missing data and the inclusion of auxiliary information.
-   *Bibliometry*: `r pkg("robustrao")` computes the Rao-Stirling diversity
    index (a well-established bibliometric indicator to measure the
    interdisciplinarity of scientific publications) with data containing
    uncategorized references.



### Links
-   [Amelia II: A Program for Missing Data](https://gking.harvard.edu/amelia)
-   [R-miss-tastic: A resource website on missing data](https://rmisstastic.netlify.com/)
