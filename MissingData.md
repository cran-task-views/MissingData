---
name: MissingData
topic: Missing Data
maintainer: Julie Josse, Imke Mayer, Nicholas Tierney, and Nathalie Vialaneix (r-miss-tastic team)
email: r-miss-tastic@clementine.wf
version: 2021-12-29
source: https://github.com/cran-task-views/MissingData/
---


Missing data are very frequently found in datasets. Base R provides a
few options to handle them using computations that involve only observed
data ( `na.rm = TRUE` in functions `mean`, `var`, ... or
`use = complete.obs|na.or.complete|pairwise.complete.obs` in functions
`cov`, `cor`, ...). The base package stats also contains the generic
function `na.action` that extracts information of the `NA` action used
to create an object.

These basic options are complemented by many packages on CRAN, which we
structure into main topics:

-   [Exploration of missing data](#exploration)
-   [Likelihood based approaches](#likelihood)
-   [Single imputation](#single)
-   [Multiple imputation](#multiple)
-   [Weighting methods](#weights)
-   [Specific types of data](#data)
-   [Specific application fields](#applications)

In addition to the present task view, this [reference website on missing
data](https://rmisstastic.netlify.com/) might also be helpful. Complementary
information might also be found in `r view("TimeSeries")` and 
`r view("SpatioTemporal")`.

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
    function `ampute` of `r pkg("mice", priority = "core")`
    through with simulations or with the `r pkg("Iscores")`
    with a KL-based scoring rule.

[**Likelihood based approaches**]{#likelihood}

-   *Methods based on the Expectation Maximization (EM) algorithm* are
    implemented in `r pkg("norm")` (using the function
    `em.norm` for multivariate Gaussian data),
    `r pkg("norm2")` (using the function `emNorm`), in
    `r pkg("cat")` (function `em.cat` for multivariate
    categorical data), in `r pkg("mix")` (function `em.mix`
    for multivariate mixed categorical and continuous data). These
    packages also implement *Bayesian approaches* (with Imputation and
    Posterior steps) for the same models (functions `da.` XXX for
    `norm`, `cat` and `mix`) and can be used to obtain imputed complete
    datasets or multiple imputations (functions `imp.` XXX for `norm`,
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
-   *Missing values in mixed-effect models and generalized linear
    models* are supported in the packages `r pkg("mdmb")`,
    `r pkg("icdGLM")` and `r pkg("JointAI")`, the
    last one being based on a Bayesian approach.
    `r pkg("brlrmr")` also handles MNAR values in response
    variable for logistic regression using an EM approach.
    `r pkg("ui")` implements uncertainty intervals for linear
    and probit regressions when the outcome is missing not at random.
-   *Missing data in item response models* is implemented in
    `r pkg("TAM")`, `r pkg("mirt")` and
    `r pkg("ltm")`.
-   *Robust* covariance estimation is implemented in the package
    `r pkg("GSE")`. Robust location and scatter estimation and
    robust multivariate analysis with missing data are implemented in
    `r pkg("rrcovNA")`.

[**Single imputation**]{#single}

-   The simplest method for missing data imputation is *imputation by
    mean* (or median, mode, \...). This approach is available in many
    packages among which `r pkg("ForImp")` and
    `r pkg("Hmisc")` that contain various proposals for
    imputing with the same value all missing instances of a variable.
-   *k-nearest neighbors* is a popular method for missing data
    imputation that is available in many packages including
    `r bioc("impute")`, `r pkg("VIM")`,
    `r pkg("GenForImp")` and
    `r pkg("yaImpute", priority = "core")` (with many different
    methods for kNN imputation, including a CCA based imputation).
    `r pkg("wNNSel")` implements a kNN based method for
    imputation in large dimensional datasets.
    `r pkg("isotree")` uses a similar approach that is based
    on similarities between samples to impute missing data with
    isolation forests.
-   *hot-deck* imputation is implemented in
    `r pkg("hot.deck", priority = "core")`,
    `r pkg("FHDI")` and `r pkg("VIM")` (function
    `hotdeck`). `r pkg("StatMatch")` uses hot-deck imputation
    to impute surveys from an external dataset.
    `r pkg("impimp")` also uses the notion of "donor" to
    impute a set of possible values, termed "imprecise imputation".
-   *Other regression based imputations* are implemented in
    `r pkg("VIM")` (linear regression based imputation in the
    function `regressionImp`). In addition,
    `r pkg("simputation")` is a general package for imputation
    by any prediction method that can be combined with various
    regression methods, and works well with the tidyverse.
    `r pkg("WaverR")` imputes data using a weighted average of
    several regressions. `r pkg("iai")` tunes optimal
    imputation based on knn, tree or SVM.
-   *Based on random forest* in `r pkg("missForest")` with a
    faster version in `r pkg("missRanger")`.
-   *Based on copula* in `r pkg("CoImp")`, in
    `r pkg("mdgc")` (mixed type data), and in
    `r pkg("sbgcop")` (semi-parametric Bayesian copula
    imputation). The latter supports multiple imputation.
-   *PCA/Singular Value Decomposition/matrix completion* is implemented
    in the package `r pkg("missMDA", priority = "core")` for
    numerical, categorical and mixed data. Heterogeneous missingness in
    a high-dimensional PCA is also addressed in
    `r pkg("primePCA")`.
    `r pkg("softImpute", priority = "core")` contains several
    methods for iterative matrix completion, as well as
    `r pkg("filling")`, `r pkg("rsparse")` and
    `r pkg("denoiseR")` for numerical variables,
    `r pkg("mimi")` that uses low rank assumptions to impute
    mixed datasets, and `r pkg("ECLRMC")` performs ensemble
    correlation based low rank matrix completion that accounts for
    correlation among samples. The package
    `r bioc("pcaMethods")` offers some Bayesian implementation
    of PCA with missing data. *NIPALS* (based on SVD computation) is
    implemented in the packages `r bioc("mixOmics")` (for PCA
    and PLS), `r pkg("ade4")`, `r pkg("nipals")` and
    `r pkg("plsRglm")` (for generalized model PLS). As a
    generalization, `r pkg("tensorBF")` implements imputation
    in 3-way tensor data. `r pkg("ROptSpace")` proposes a
    matrix completion method under low-rank assumption and collective
    matrix factorization for imputation using Bayesian matrix completion
    for groups of variables (binary, quantitative, Poisson). Similarly,
    `r pkg("cmfrec")` proposes different version of matrix
    imputation based on matrix factorization. Finally, imputation for
    groups is also available in `r pkg("missMDA")` in the
    function `imputeMFA`.
-   *Imputation for non-parametric regression by wavelet shrinkage* is
    implemented in `r pkg("CVThresh")` using solely
    maximization of the h-likelihood.
-   `r pkg("mi")` and `r pkg("VIM")` also provide
    diagnostic plots to *evaluate the quality of imputation.*

[**Multiple imputation**]{#multiple}

Some of the above mentioned packages can also handle multiple
imputations.

-   `r pkg("Amelia", priority = "core")` implements Bootstrap
    multiple imputation using EM to estimate the parameters, for
    quantitative data it imputes assuming a Multivariate Gaussian
    distribution. In addition, AmeliaView is a GUI for
    `r pkg("Amelia")`, available from the [Amelia web
    page](https://gking.harvard.edu/amelia).
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
-   `r pkg("accelmissing")` provides multiple imputation with
    the zero-inflated Poisson lognormal model for missing count values
    in accelerometer data.

In addition, `r pkg("mitools")` provides a generic approach to
handle multiple imputation in combination with any imputation method.
And `r pkg("NADIA")` provides a uniform interface to compare
the performances of several imputation algorithms.

[**Weighting methods**]{#weights}

-   *Computation of weights* for observed data to account for unobserved
    data by *Inverse Probability Weighting (IPW)* is implemented in
    `r pkg("ipw")`. IPW is also used for quantile estimations
    and boxplots in `r pkg("IPWboxplot")`.
-   *Doubly Robust Inverse Probability Weighted Augmented GEE Estimator
    with missing outcome* is implemented in
    `r pkg("CRTgeeDR")`.

[**Specific types of data**]{#data}

-   *Longitudinal data / time series and censored data*: Imputation for
    time series is implemented in
    `r pkg("imputeTS", priority = "core")` and
    `r pkg("imputePSF")`. Other packages, such as
    `r pkg("forecast")`, `r pkg("spacetime")`,
    `r pkg("timeSeries")`, `r pkg("xts")`,
    `r pkg("prophet")`, `r pkg("stlplus")` or
    `r pkg("zoo")`, are dedicated to time series but also
    contain some (often basic) methods to handle missing data (see also
    `r view("TimeSeries")`). To help fill down missing values
    for time series, the `r pkg("padr")` and
    `r pkg("tsibble")` packages provide methods for imputing
    implicit missing values. Imputation of time series based on Dynamic
    Time Warping is implemented in `r pkg("DTWBI")` for
    univariate time series and in `r pkg("DTWUMI")` or in
    `r pkg("FSMUMI")` for multivariate ones.
    `r pkg("naniar")` also imputes data below the range for
    exploratory graphical analysis with the function `impute_below`.
    `r pkg("TAR")` implements an estimation of the
    autoregressive threshold models with Gaussian noise and of
    positive-valued time series with a Bayesian approach in the presence
    of missing data and `r pkg("BMTAR")` has a similar
    approach for multivariate time series. `r pkg("swgee")`
    implements a probability weighted generalized estimating equations
    method for longitudinal data with missing observations and
    measurement error in covariates based on SIMEX.
    `r pkg("icenReg")` performs imputation for censored
    responses for interval data. `r pkg("imputeTestbench")`
    proposes tools to benchmark missing data imputation in univariate
    time series. On a related topic, `r pkg("imputeFin")`
    handles imputation of missing values in financial time series using
    AR models or random walk.
-   *Spatial data*: Imputation for spatial data is implemented in
    `r pkg("phylin")` using interpolation with spatial
    distance weights or kriging. `r pkg("gapfill")` is
    dedicated to satellite data. Geostatistical interpolation of data
    with irregular spatial support is implemented in
    `r pkg("rtop")` and in `r pkg("areal")` that
    estimates values for overlapping but incongruent polygon features.
-   *Spatio-temporal data* (see also `r view("SpatioTemporal")`): Imputation 
    for spatio-temporal data is implemented in the package `r pkg("cutoffR")` 
    using different methods as knn and SVD and in `r pkg("CircSpaceTime")` for 
    circular data using kriging. 
    Similarly, `r pkg("reddPrec")` imputes missing values in daily precipitation
    time series accross different locations. 
    Estimation and prediction for spatio-temporal data with missing values is 
    implemented in `r pkg("StempCens")` with a SAEM approach that approximates 
    EM when the E-step does not have an analytic form.
-   *Graphs/networks*: Imputation for graphs/networks is implemented in
    the package `r pkg("dils")` to impute missing edges.
    `r pkg("PST")` provides a framework for analyzing
    Probabilistic Suffix Trees, including functions for learning and
    optimizing VLMC (variable length Markov chains) models from sets of
    individual sequences possibly containing missing values.
    `r pkg("missSBM")` imputes missing edges in Stochastic
    Block models and `r pkg("cassandRa")` predicts possible
    missing links with different stochastic network models.
    `r pkg("cglasso")` implements an extension of the
    Graphical Lasso with censored and missing values and
    `r pkg("mgm")` implements estimation of time-varying
    k-order mixed graphical models.
-   *Imputation for contingency tables* is implemented in
    `r pkg("lori")` that can also be used for the analysis of
    contingency tables with missing data.
-   *Imputation for compositional data (CODA)* is implemented in in
    `r pkg("zCompositions")` (various imputation methods for
    zeros, left-censored and missing data).
-   *Imputation for diffusion processes* is implemented in
    `r pkg("DiffusionRimp")` by imputing missing sample paths
    with Brownian bridges.
-   *Imputation for meta-analyses* of binary outcomes is provided in
    `r pkg("metasens")`.
-   `r pkg("experiment")` handles missing values in
    experimental design such as randomized experiments with missing
    covariate and outcome data, matched-pairs design with missing
    outcome.
-   `r pkg("cdparcoord")` provides tools to handle missing
    values in parallel coordinates settings.
-   `r pkg("dejaVu")` performs multiple imputation of
    recurrent event data based on a negative binomial regression model.

[**Specific application fields**]{#applications}

-   *Genetics*: Analyses of Case-Parent Triad and/or Case-Control Data
    with SNP haplotypes is implemented in `r pkg("Haplin")`,
    where missing genotypic data are handled with an EM algorithm.
    `r pkg("FamEvent")` and `r bioc("snpStats")`
    implement imputation of missing genotypes, respectively with an EM
    algorithm and a nearest neighbor approach. Imputation for genotype
    and haplotype is implemented in `r pkg("alleHap")` using
    solely deterministic techniques on pedigree databases; imputation of
    missing genotypes is also implemented in `r pkg("QTLRel")`
    that contains tools for QTL analyses. Tools for Hardy-Weinberg
    equilibrium for bi- and multi-allelic genetic marker data are
    implemented in `r pkg("HardyWeinberg")`, where genotypes
    are imputed with a multinomial logit model.
    `r pkg("StAMPP")` computes genomic relationship when SNP
    genotype datasets contain missing data and
    `r pkg("PSIMEX")` computes inbreeding depression or
    heritability on pedigree structures affected by missing paternities
    with a variant of the SIMEX algorithm.
-   *Genomics*: Imputation for dropout events ( *i.e.* , under-sampling
    of mRNA molecules) in single-cell RNA-Sequencing data is implemented
    in `r pkg("DrImpute")` and `r pkg("Rmagic")`.
    `r pkg("RNAseqNet")` uses hot-deck imputation to improve
    RNA-seq network inference with an auxiliary dataset.
-   *Epidemiology*: `r pkg("idem")` implements a procedure
    for comparing treatments in clinical trials with missed visits or
    premature withdrawal. `r pkg("InformativeCensoring")`
    implements multiple imputation for informative censoring.
    `r pkg("pseval")` evaluates principal surrogates in a
    single clinical trial in the presence of missing counterfactual
    surrogate responses. `r pkg("sievePH")` implements
    continuous, possibly multivariate, mark-specific hazard ratio with
    missing values in multivariate marks using an IPW approach.
-   *Causal inference*: Causal inference with interactive fixed-effect
    models is available in `r pkg("gsynth")` with missing
    values handled by matrix completion. `r pkg("MatchThem")`
    matches multiply imputed datasets using several matching methods,
    and provides users with the tools to estimate causal effects in each
    imputed datasets. `r pkg("grf")` offers treatment effect
    estimation with incomplete confounders and covariates under modified
    unconfoundedness assumptions.
-   *Scoring*: Basic methods (mean, median, mode, \...) for imputing
    missing data in scoring datasets are proposed in
    `r pkg("scorecardModelUtils")`.
-   *Preference models*: Missing data in preference models are handled
    with a *Composite Link* approach that allows for MCAR and MNAR
    patterns to be taken into account in `r pkg("prefmod")`.
-   *Health economy*: `r pkg("missingHE")` implements models
    for health economic evaluations with missing outcome data.
-   *Administrative records / Surveys*: `r pkg("fastLink")`
    provides a Fellegi-Sunter probabilistic record linkage that allows
    for missing data and the inclusion of auxiliary information.
-   *Regression and classification*: `r pkg("eigenmodel")`
    handles missing values in regression models for symmetric relational
    data. `r pkg("randomForest")`, `r pkg("grf")`
    and `r pkg("StratifiedRF")` handle missing values in
    predictors in various random forest based methods.
    `r pkg("misaem")` handles missing data in linear and
    logistic regression and allows for model selection.
    `r pkg("psfmi")` provides a framework for model selection
    for various linear models in multiply imputed datasets.
    `r pkg("naivebayes")` provides an efficient implementation
    of the naive Bayes classifier in the presence of missing data.
    `r pkg("plsRbeta")` implements PLS for beta regression
    models with missing data in the predictors. `r pkg("lqr")`
    provides quantile regression estimates based on various
    distributions in the presence of missing values and censored data.
-   *Clustering*: `r pkg("biclustermd")` handles missing data
    in biclustering. `r pkg("RMixtComp")`,
    `r pkg("MGMM")` and `r pkg("mixture")` fit
    various mixture models in the presence of missing data.
-   *Tests* for two-sample paired missing data are implemented in
    `r pkg("robustrank")`.
-   *Bibliometry*: `r pkg("robustrao")` computes the Rao-Stirling diversity
    index (a well-established bibliometric indicator to measure the
    interdisciplinarity of scientific publications) with data containing
    uncategorized references.



### Links
-   [Amelia II: A Program for Missing Data](https://gking.harvard.edu/amelia)
-   [R-miss-tastic: A resource website on missing data](https://rmisstastic.netlify.com/)
