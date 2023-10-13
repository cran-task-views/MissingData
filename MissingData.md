---
name: MissingData
topic: Missing Data
maintainer: Julie Josse, Imke Mayer, Nicholas Tierney, Nathalie Vialaneix
email: r-miss-tastic@clementine.wf
version: 2023-06-20
source: https://github.com/cran-task-views/MissingData/
---

Missing data are very frequently found in datasets. Base R provides a few
options to handle them using computations that involve only observed data
(`na.rm = TRUE` in functions `mean`, `var`, ... or
`use = complete.obs|na.or.complete|pairwise.complete.obs` in functions `cov`,
`cor`, ...). The base package `stats` also contains the generic function
`na.action` that extracts information of the `NA` action used to create an
object. In addition, the package `r pkg("ie2misc")` contains a dyadic operator
`+` that behaves differently than the original `+` operator regarding missing
data.

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
`r view("SpatioTemporal")`, `r view("Survival")`, and
`r view("OfficialStatistics")`. Note that most packages covering temporal, and
spatio-temporal interpolation and censored data are not covered by the Missing
Data task view.

If you think we have missed some important packages in this list, please
e-mail the maintainers or submit an issue or pull request in the GitHub
repository linked above.

[**Exploration of missing data**]{#exploration}

-   *Manipulation of missing data* is implemented in the packages
    `r pkg("sjmisc")`, `r pkg("sjlabelled")`, `r pkg("retroharmonize")`,
    `r pkg("mde")` (also providing basic functions to explore missingness
    patterns), `r pkg("tidyr")` (which abides by
    [tidyverse](https://tidyverse.org) principles), and `r pkg("declared")`. In 
    addition, `r pkg("memisc")` provides definable missing values, along with
    infrastructure for the management of survey data and variable labels. More
    specifically, `r pkg("fauxnaif")` converts given values to `NA` and 
    `r pkg("fillr")` fill missing values in vectors according to simple 
    predefined rules.\
    `r pkg("roperators")` provides string arithmetic, reassignment operators,
    logical operators that handle missing values.
-   *Missing data patterns* can be identified and explored using the
    packages `r pkg("mi")` (and its GUI `r pkg("migui")`), `r pkg("wrangle")`, 
    `r pkg("DescTools")`,
    `r pkg("dlookr")`, and `r pkg("naniar", priority = "core")`. 
    `r pkg("daqapo")` is a generic data quality toolbox that can also be used to
    identify missing data. More specifically, `r pkg("ggmice")` produces plots
    for the `r pkg("mice")` imputation workflow and can be used for missing data
    exploration and evaluation of imputation quality.
-   *Graphics that describe distributions and patterns of missing data*
    are implemented in `r pkg("VIM", priority = "core")` (which has a Graphical
    User Interface, VIMGUI, currently archived on CRAN) and `r pkg("naniar")`
    (which abides by [tidyverse](https://tidyverse.org) principles).
-   *Tests of the MAR assumption (versus the MCAR assumption)*: Little's test
    for the MCAR assumption is implemented in `r pkg("misty")`. Other approaches
    are also available elsewhere: `r pkg("RBtest")` proposes a regression based
    approach to test for missing data mechanisms and `r pkg("PKLMtest")`
    implements a KL-based test for MCAR.\
    In addition, `r pkg("isni")` tests sensitivity to the
    ignorability assumption by computing the index of local sensitivity to
    nonignorability.
-   *Evaluation*: `r pkg("missCompare")` and `r pkg("missMethods")` offer an
    entire framework to compare different imputation strategies (with
    diagnostics and visualizations). The package `r pkg("Iscores")` can also be
    useful to evaluate imputation quality using a KL-based scoring rule.\
    Simulations to evaluate imputation qualities can be performed using the
    function `ampute` of `r pkg("mice", priority = "core")`, the package
    `r pkg("simFrame")`, which proposes a very general framework for simulations,
    or the package `r pkg("simglm")`, which simulates data and missing values in
    simple and generalized linear regression models. Similarly,
    `r pkg("imputeTestbench")` provides a benchmark to evaluate univariate time
    series imputation. \
    In addition, `r pkg("mi")` and `r pkg("VIM")` also provide diagnostic plots
    that can help evaluate imputation quality.

[**Likelihood based approaches**]{#likelihood}

-   *Methods based on the Expectation Maximization (EM) algorithm* are
    implemented in `r pkg("norm")`, `r pkg("norm2")`, and `r pkg("mvnmle")` for
    multivariate normal datasets,
    in `r pkg("cat")` (function `em.cat` for multivariate categorical data), in
    `r pkg("mix")` (function `em.mix` for multivariate mixed categorical and
    continuous data). These packages also implement *Bayesian approaches* (with
    Imputation and Posterior steps) for the same models (functions `da.`XXX for
    `norm`, `cat` and `mix`) and can be used to obtain imputed complete
    datasets or multiple imputations (functions `imp.`XXX for `norm`,
    `cat` and `mix`), once the model parameters have been estimated.
    `r pkg("monomvn")` proposes similar methods for multivariate normal and
    Student distributions when the missingness pattern is monotonic.\
    `r pkg("CensMFM")`, `r pkg("imputeMulti")`, and `r pkg("MMDai")` extend 
    these methods by using an EM approach to fit different mixtures of 
    multivariate missing data for numeric or 
    categorical data. `r pkg("RMixtCompIO")` is a complete library of mixture 
    models that handles missing data and is based on the C++ library `MixtComp`.
    It can be used in combination with `r pkg("RMixtCompUtilities")`, which 
    provides various graphical, getters, and utility functions.\
    Hierarchical Gaussian and probit models with missing covariate values are
    implemented in `r pkg("ppmSuite")`. `r pkg("PReMiuM")` implements
    Dirichlet process mixture models (regression models linking the response to
    covariates through cluster membership) with missing covariate values.\
    `r pkg("imputeR")` is also using an EM based imputation framework that
    offers several different algorithms, including Lasso, tree-based models or
    PCA. In addition, `r pkg("TestDataImputation")` implements imputation based
    on EM estimation (and other simpler imputation methods) that are well suited
    for dichotomous and polytomous tests with item responses.
-   *Multiple imputation* is performed using Maximum Likelihood Multiple
    Imputation in `r pkg("mlmi")`.
-   *Full Information Maximum Likelihood* (also known as "direct maximum
    likelihood" or "raw maximum likelihood") is available in `r pkg("lavaan")`
    (and in its extension `r pkg("semTools")`), `r pkg("OpenMx")`,
    `r pkg("rsem")`, and `r pkg("simsem")` for handling missing data in
    structural equation modeling.
-   *Bayesian approaches* for handling missing values in model based
    clustering with variable selection is available in `r pkg("VarSelLCM")`.
    The package also provides imputation using the posterior mean.
-   *Missing values in generalized linear models* can be handled with package
    `r pkg("mdmb")` for various families. `r pkg("JointAI")` implements
    Bayesian approaches for generalized linear mixed models and `r pkg("bild")`
    implements logistic regression with mixed effects for binary longitudinal
    data allowing missing values. `r pkg("ClusPred")` also handles missing 
    values in mixed model with a fixed group effect, when the group variable is 
    missing.\
    `r pkg("brlrmr")` proposes a method to 
    reduce bias in estimating logistic regressions with missing response.
-   *Missing data in item response models* (including Rasch models and
    extensions) is implemented in `r pkg("TAM")`, `r pkg("mirt")`,
    `r pkg("eRm")`, and `r pkg("ltm")` for univariate or multivariate responses.
    `r pkg("LNIRT")` also addresses these models but allows missing values to
    be specified as "missing-by-design" and `r pkg("MLCIRTwithin")` includes
    latent-class models.
-   *Missing values in outcome of regression models* is handled in 
    `r pkg("mreg")`.

[**Single imputation**]{#single}

-   The simplest method for missing data imputation is *imputation by mean* (or
    median, mode, \...). This approach is available in many packages among which
    `r pkg("Hmisc")` that contains various proposals for
    imputing with the same value all missing instances of a variable.
-   *Generic packages*: The packages `r pkg("VIM")` and `r pkg("filling")`
    contain several popular methods for missing value imputation (including some
    listed in the sections dedicated to specific methods as listed below).
    In addition, `r pkg("simputation")` is a general package for imputation
    by any prediction method that can be combined with various regression
    methods, and works well with the [tidyverse](https://tidyverse.org).
-   *k-nearest neighbors* is a popular method for missing data imputation that
    is available in many packages including the main packages
    `r pkg("yaImpute", priority = "core")` (with many different methods for
    kNN imputation, including a CCA based imputation) and `r pkg("VIM")`. It
    is also available in `r bioc("impute")` (where it is oriented toward
    microarray imputation).\
    `r pkg("isotree")` uses a similar approach to impute missing values, which is
    based on similarities between samples and isolation forests.
-   *Hot-deck* imputation is implemented in the package
    `r pkg("hot.deck", priority = "core")`, with various possible settings
    (including multiple imputation). It is also available in `r pkg("VIM")`
    (function `hotdeck`) and a fractional version (using weights) is provided
    in `r pkg("FHDI")`. `r pkg("StatMatch")` also uses hot-deck imputation
    to impute surveys from an external dataset.\
    Similarly, `r pkg("impimp")` uses the notion of a "donor" to impute a set of
    possible values, termed "imprecise imputation".
-   Imputation *based on random forest* is implemented in `r pkg("missForest")`
    with a faster version in `r pkg("missRanger")`. `r pkg("Rforestry")` extend
    this method with variants of the original random forest method.
-   *Other regression based imputations* are implemented in `r pkg("VIM")`
    (linear regression based imputation in the function `regressionImp`).
    `r pkg("iai")` tunes optimal imputation based on knn, tree or SVM and
    `r pkg("SurrogateRegression")` uses bivariate regressions to perform 
    estimation and inference on partially missing target outcomes.
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
    and `r pkg("tensorBF")` (for 3-way tensor data).
    *Low rank based imputation* is provided in
    `r pkg("softImpute", priority = "core")`, which contains several
    methods for iterative matrix completion. `r pkg("eimpute")` implements 
    an efficient imputation methods based on low rank approximation of large
    matrices. Low rank imputation methods are also available in
    the very general package `r pkg("rsparse")`, which contains various tools
    for sparse matrices. Variants based on low rank assumptions are available
    in `r pkg("denoiseR")`, in `r pkg("mimi")`, in `r pkg("ECLRMC")` and
    `r pkg("CMF")` (for ensemble matrix completion), and in `r pkg("ROptSpace")`
    (with a computationally efficient approach).
-   Imputation for *categorical variables* is proposed in `r pkg("NIMAA")` based
    on data mining and simple rules. `r pkg("OTrecod")` can also impute 
    categorical variables by using information shared by two databases and a 
    method based on Optimal Transport.
-   Imputation *based on copula* is implemented in `r pkg("CoImp")` with a
    semi-parametric imputation procedure and in `r pkg("mdgc")` using Gaussian
    copula for mixed data types.
-   Imputation *based on self-organizing maps* is provided in 
    `r pkg("SOMbrero")` and `r pkg("missSOM")`.
-   Imputation *based on validation rules (deductive methods)* is implemented in
    `r pkg("deductive")`.

[**Multiple imputation**]{#multiple}

Some of the above mentioned packages can also handle multiple
imputations.

-   `r pkg("Amelia", priority = "core")` implements *Bootstrap multiple
    imputation using EM* to estimate the parameters, for quantitative data it
    imputes assuming a Multivariate Gaussian distribution. In addition,
    [AmeliaView](https://cran.r-project.org/web/packages/Amelia/vignettes/ameliaview.html)
    is a GUI for `r pkg("Amelia")`, available from the
    [Amelia web page](https://gking.harvard.edu/amelia). 
    `r pkg("FastImputation")` provides a fast approximation of the imputation 
    process used in `r pkg("Amelia")`.\
    `r pkg("NPBayesImputeCat")` also implements multiple imputation by joint
    modeling for categorical variables but using a Bayesian approach.
-   `r pkg("mi")`, `r pkg("mice")`, and `r pkg("smcfcs")` implement *multiple
    imputation by Chained Equations*. Other packages are based on or extend
    `r pkg("mice")`, like `r pkg("miceFast")`, which provides an
    alternative implementation of mice imputation methods using object
    oriented style programming and C++, `r pkg("bootImpute")`, which performs
    bootstrap based imputations and analyses of these imputations, and
    `r pkg("miceRanger")`, `r pkg("CALIBERrfimpute")`, and `r pkg("RfEmpImp")`,
    which all perform
    multiple imputation by chained equations using random forests.
-   *Multiple imputation based on Markov models* is proposed in 
`r pkg("niaidMI")`.
-   *Dealing with multiply imputed datasets*: `r pkg("mitools")` provides a 
    generic approach to handle multiple imputation in combination with any
    imputation method, `r pkg("cobalt")` computes balance tables and plots for 
    multiply imputed datasets, `r pkg("SynthTools")` provides confidence 
    intervals for multiply imputed datasets, `r pkg("miceafter")` allows
    different types of statistical analyses and pooling after multiple 
    imputation.\
    In addition, `r pkg("mitools")` provides a generic approach to handle 
    multiple imputation in combination with any imputation method, 
    `r pkg("cobalt")` computes balance tables and plots for multiply imputed 
    datasets, `r pkg("basecamb")` fits models from multiply imputed datasets 
    in a generic way, and `r pkg("bucky")` automatically computes summary 
    information for estimates computed on multiply imputed data sets in a social
    sciences context.
-   `r pkg("missMDA")` implements multiple imputation based on *SVD methods*.
-   `r pkg("hot.deck")` implements *hot-deck*-based multiple imputation.
-   `r pkg("rMIDAS")` implements multiple imputation based on *denoising
    auto-encoders*.
-   *Multilevel imputation*: Multilevel multiple imputation is implemented in
    `r pkg("jomo", priority = "core")`, `r pkg("mice")`,
    `r pkg("miceadds")`, `r pkg("micemd")`, `r pkg("mitml")`, and
    `r pkg("pan")`.
-   `r pkg("gerbil")` implements multiple imputation using latent joint 
    multivariate normal models.
-   `r pkg("Qtools")` and `r pkg("miWQS")` implement multiple imputation based
    on *quantile regression*.
-   `r pkg("lodi")` implements the *imputation of observed values below the
    limit of detection* (LOD) via censored likelihood multiple imputation
    (CLMI).
-   *Evaluation*: `r pkg("NADIA")` provides a uniform interface to compare the 
    performance of several imputation algorithms.

[**Weighting methods**]{#weights}

-   *Computation of weights* for observed data to account for unobserved
    data by *Inverse Probability Weighting (IPW)* is implemented in
    `r pkg("ipw")` and `r pkg("iWeigReg")`. `r pkg("nawtilus")` also proposes
    IPW computation but utilizing estimating equations suitable for a specific
    pre-specified parameter of interest.
-   IPW is also for *quantile estimations and boxplots* in
    `r pkg("IPWboxplot")`.
-   *Doubly Robust Inverse Probability Weighted Augmented GEE Estimator
    with missing outcome* is implemented in `r pkg("CRTgeeDR")`.
-   *IPW for time-course missing data* is implemented in `r pkg("MIIPW")`.

[**Specific types of data**]{#data}

-   *Longitudinal data / time series data*: Imputation for time series is
    implemented in `r pkg("imputeTS", priority = "core")`. Other packages, such
    as `r pkg("forecast")`, `r pkg("spacetime")`, `r pkg("timeSeries")`,
    `r pkg("xts")`, `r pkg("prophet")`, `r pkg("stlplus")`, or `r pkg("zoo")`,
    are dedicated to time series but also contain some (often basic) methods to
    handle missing data (see also `r view("TimeSeries")`). Based on tidy
    principle, the `r pkg("padr")` and `r pkg("tsibble")` also provide methods
    for imputing missing values in time series. Similarly, `r pkg("DTSg")`
    offers basic functionality for missing value description and imputation in
    time series based on the fast `data.table` framework.\
    More specific methods are implemented in other packages: imputation of time
    series based on Dynamic Time Warping is implemented in the family of
    packages `r pkg("DTWBI")`, `r pkg("DTWUMI")`, and `r pkg("FSMUMI")` for
    univariate and multivariate time series. `r pkg("BMTAR")` provides
    an estimation of the autoregressive threshold models with Gaussian noise
    using a Bayesian approach in the presence of missing data in multivariate
    time series. `r pkg("swgee")` implements an IPW approach
    for longitudinal data with missing observations. `r pkg("tsrobprep")`
    implements imputation of missing values using a robust decomposition of the
    time series. `r pkg("brokenstick")` handles missing at random data in 
    irregular time series with a brokenstick approach. `r pkg("TRMF")` uses
    temporally regularized matrix factorizations to impute values in 
    high-dimensional time series.\
    For more specific time series, `r pkg("cold")` fits longitudinal count  
    models from data with missing values.\
    Estimation of extremal indexes in time series is implemented in 
    `r pkg("exdex")` with K-gaps and D-gaps models that can accommodate with 
    missing values.
-   *Markov models*: `r pkg("hhsmm")` includes various methods for hidden hybrid
    Markov and semi-Markov models that can accomodate missing data.
-   *Spatial data*: Imputation for spatial data is implemented in the package
    `r pkg("rtop")`, which performs geostatistical interpolation of irregular
    areal data, and in `r pkg("areal")`, which performs areal weighted
    interpolation using a tidyverse data management. `r pkg("RcppCensSpatial")`
    estimates parameters in linear spatial models with missing data using EM, 
    SAEM, or MCEM.\
    Interpolation of spatial data based on genetic distances is also
    available in `r pkg("phylin")`.
-   *Spatio-temporal data* (see also `r view("SpatioTemporal")`): Imputation
    for spatio-temporal data is implemented in the package `r pkg("StempCens")`
    with a SAEM approach that approximates EM when the E-step does not have an
    analytic form.\
    From an application perspective, `r pkg("gapfill")` is dedicated to the
    imputation of satellite data observed at equally-spaced points in time and
    `r pkg("stfit")` uses Functional Principal Analysis by Conditional 
    Estimation to impute missing pixels in satellite data.
    `r pkg("momentuHMM")` is dedicated to the analysis of telemetry
    data using generalized hidden Markov models (including multiple imputation
    for missing data).
-   *Distance matrices*: Imputation for Euclidean distance matrix is implemented
    in `r pkg("edmcr")`, using different optimization approaches.
-   *Graphs/networks*: `r pkg("missSBM")` imputes missing edges in Stochastic
    Block models, `r pkg("cglasso")` implements an extension of the Graphical
    Lasso inference from censored and missing value measurements, and
    `r pkg("bnstruct")` provides an extension of various methods for Bayesian
    network inference from data with missing values. Oriented toward inference
    of species community networks, `r pkg("eicm")` uses an extension of
    binomial GLM that handles missing values and `r pkg("robber")` is based on
    stochastic block models and also handles missing values. `r pkg("rnmamod")` 
    includes functions to explore network meta-analysis with missing participant
    outcome data in clinical trials.
-   *Imputation for contingency tables* is implemented in `r pkg("lori")` that
    can also be used for the analysis of contingency tables with missing data.
-   *Imputation for compositional data (CODA)* is implemented in
    `r pkg("robCompositions")` and `r pkg("zCompositions")` (various imputation
    methods for zeros, left-censored and missing data).
-   *Rank models* with partially missing rankings are handled in
    `r pkg("BayesMallows")` with Bayesian methods, and in `r pkg("irrNA")` to
    compute inter-rater reliability and concordance.
-   *Experimental design*: `r pkg("experiment")` handles missing values in
    experimental design such as randomized experiments with missing covariate
    and outcome data, and matched-pairs design with missing outcome.
-   *Recurrent events*: `r pkg("dejaVu")` performs multiple imputation of
    recurrent event data based on a negative binomial regression model.

[**Specific tasks**]{#tasks}

-   *Regression and classification*: many different supervised methods can
    accommodate the presence of missing values. `r pkg("randomForest")`,
    `r pkg("grf")`, and `r pkg("StratifiedRF")` handle missing values in
    predictors in various random forest based methods.\
    `r pkg("misaem")`
    handles missing data in linear and logistic regression and allows for model
    selection. `r pkg("psfmi")` also provides a framework for model selection
    for various linear models in multiply imputed datasets and `r pkg("flare")` 
    accommodates missing values in some models related to Lasso regression.\
    `r pkg("naivebayes")` provides an efficient implementation of the naive
    Bayes classifier in the presence of missing data. `r pkg("plsRbeta")`
    implements PLS for beta regression models with missing data in the
    predictors. `r pkg("lqr")` provides quantile regression estimates based on
    various distributions in the presence of missing values and censored data.
    `r pkg("eigenmodel")` handles missing values in regression models for
    symmetric relational data. 
-   *Clustering*: `r pkg("biclustermd")` handles missing data in biclustering.
    `r pkg("RMixtComp")`, `r pkg("MGMM")`, `r pkg("mixture")`, and 
    `r pkg("MixtureMissing")` fit various
    mixture models in the presence of missing data. `r pkg("ClustImpute")` deals
    with missing values in k-means clustering. `r pkg("gbmt")` performs 
    clustering to identify similar trajectories in multivariate longitudinal 
    data containing missing values. `r pkg("LUCIDus")` performed clustering from
    multiple omics when some omics are missing. `r pkg("miclust")` handles 
    multiple imputation in clustering.
-   *Tests* for two-sample paired missing data are implemented in
    `r pkg("robustrank")`, `r pkg("IncomPair")`, and `r pkg("MKinfer")`, the 
    latter is based on multiple
    imputed datasets. Reliability of tests for data with missing values is
    assessed with a Bayesian approach in `r pkg("brxx")`.
-   *Meta-analysis*: `r pkg("metavcov")` offers a collection of functions,
    including multiple imputations for missing data, in multivariate
    meta-analyses. `r pkg("metansue")` can perform meta-analysis with some 
    missing (unreported) effects. 
    More specifically, imputation for *meta-analyses* of binary
    outcomes is provided in `r pkg("metasens")` and `r pkg("NMADiagT")` provides
    a Bayesian analysis using network meta-analysis of dose response studies
    in which MNAR missing values are accounted for.
-   *Sensitivity analysis* and confidence intervals with non-ignorable 
    missingness patterns are handled in `r pkg("ui")`.
-   *Outlier detection* (and robust analysis) in the presence of missing values
    is implemented in `r pkg("GSE")` and `r pkg("rrcovNA")`.
-   *ROC estimation* in the presence of missing values is available in
    `r pkg("bcROCsurface")` for ROC surface and in `r pkg("BLOQ")` for left
    censored data.
-   *Mediation analysis* in the presence of missing values is implemented in 
    `r pkg("bmem")` and `r pkg("bmemLavaan")`, the latter designed to handle 
    non-normal data. `r pkg("paths")` uses an imputation method for the 
    estimation of path specific effects in causal mediation analysis.
-   *Composite Indicator* can be imputed with the package `r pkg("COINr")`.
-   *Fuzzy logic*: `r pkg("lfl")`contains basic fuzzy-related algebraic 
    functions capable of handling missing values for fuzzy logic.
-   *ODE*: An implementation of the parameter cascade method for estimating 
    ordinary differential equation models with missing or complete observations
    is provided in the package `r pkg("pCODE")`.

[**Specific application fields**]{#applications}

-   *Genetics*: Imputation of SNP data is implemented in `r pkg("alleHap")`
    (using solely deterministic techniques based on pedigree data), in
    `r pkg("QTLRel")` (using information on flanking SNPs), in
    `r bioc("snpStats")` (using a nearest neighbor approach), in
    `r pkg("HardyWeinberg")` (using multiple imputations with a multinomial
    model based on allele intensities and/or flanking SNPs). In addition, 
    `r pkg("SNPassoc")` and `r pkg("SNPfiltR")` offer functions to explore 
    missing SNPs.\
    `r pkg("qgtools")` includes linear mixed models and resampling techniques
    for quantitative genetics analyses in the presence of missing data. EM
    algorithm is used to compute genetic statistics for population in the
    presence of missing SNP in `r pkg("StAMPP")` and to fit 
    genotype-to-phenotype models in `r pkg("FamEvent")`, `r pkg("hapassoc")`,
    and `r pkg("Haplin")`.\
    Finally, `r pkg("FILEST")` is used to simulate SNP datasets with outlying
    individuals and missing values.
-   *Phylogeny*: Imputation of missing data for phylogeny is implemented in
    `r pkg("Rphylopars")` with different evolutionary models. Simulation of
    incomplete phylogeny can be performed with `r pkg("TreeSim")`.
-   *Genomics*: Imputation for dropout events (i.e., under-sampling of mRNA
    molecules) in single-cell RNA-Sequencing data is implemented
    in `r pkg("DrImpute")`, `r pkg("SAVER")`, and
    `r pkg("iCellR")`, and is based, respectively, on clustering of cells,
    Markov affinity graph, an empirical Bayes approach, and k-nearest neighbors.
    The first three packages are used and combined in `r bioc("scRecover")` and
    `r bioc("ADImpute")` and the last one can also handle other types of
    single-cell data, such as scATAC-Seq or CITE-Seq.\
    `r pkg("RNAseqNet")` uses hot-deck imputation to improve RNA-seq network
    inference with an auxiliary dataset.
-   *Chemometrics*: `r pkg("imp4p")`, `r pkg("wrProteo")`, `r pkg("mi4p")`,
    `r pkg("imputeLCMD")` and `r pkg("aLFQ")` use imputation for protein 
    quantification from LC-MS/MS data. The first three use multiple imputation 
    and `r pkg("imp4p")`, `r pkg("wrProteo")`, and `r pkg("imputeLCMD")` can 
    work under an MNAR mechanism. Imputation for quantified metabolomics data is
    implemented in `r pkg("lilikoi")` with a k-NN approach.\
    Imputation of data under detection limit for NIR spectra is provided in
    `r pkg("NIRStat")` for standard analyses of NIR time series.
-   *Epidemiology*: `r pkg("bayesCT")` implements various methods for simulation
    and analysis of clinical trials in a Bayesian framework that allows for
    handling and imputation of missing data. `r pkg("sanon")` implements a
    method for analysis of randomized clinical trials with strata that can
    handle MCAR data. `r pkg("didimputation")` implements treatment effect 
    estimation and pre-trend testing in diff-in-diff designs with an imputation 
    approach. `r pkg("diyar")` implements record linkage and epidemiological case
    definitions while addressing missing data across linkage stages.\
    More specifically, 
    `r pkg("InformativeCensoring")` implements multiple imputation for
    informative censoring. `r pkg("pseval")` evaluates principal surrogates in
    a single clinical trial in the presence of missing counterfactual surrogate
    responses. `r pkg("sievePH")` implements continuous, possibly multivariate,
    mark-specific hazard ratio with missing values in multivariate marks using
    an IPW approach. `r pkg("icenReg")` performs imputation for censored
    responses for interval data.
-   *Health*: `r pkg("missingHE")` implements models for health economic
    evaluations with missing outcome data. `r pkg("accelmissing")` provides
    multiple imputation with the zero-inflated Poisson lognormal model for
    missing count values in accelerometer data. `r pkg("CGManalyzer")` provides
    tools for the analysis of continuous glucose monitoring that can handle 
    missing data. `r pkg("qpNCA")` implements imputation for noncomportmental
    pharmacokinetic longitudinal data mostly using interpolation methods.
-   *Morphometry*:  `r pkg("LOST")` can be used to simulate missing morphometric
    data randomly, with taxonomic bias and with anatomical biases.
-   *Environment*: `r pkg("AeRobiology")` imputes missing data in
    aerobiological datasets imported from aerobiological public databases.
    `r pkg("climatol")` implements functions for missing data filling of 
    climatological series. 
    `r pkg("QUALYPSO")` can handle missing data and provides unbiased estimates
    of climate change responses for incomplete ensembles of climate projections.
-   *Social sciences*: `r pkg("coefficientalpha")` computes coefficients Alpha,
    social, behavioral and education sciences, in the presence of missing data.
-   *Causal inference*: Various methods for causal inference with missing data
    are implemented in `r pkg("targeted")`, using augmented IPW estimators.
    Causal inference with interactive fixed-effect models is available in
    `r pkg("gsynth")`, with missing values handled by matrix completion, and in
    `r pkg("dosearch")`, via extension of do-calculus to missing data. 
    `r pkg("R6causal")` implements R6 class for structural causal models where
    the missing data mechanism can be specified.
    `r pkg("MatchThem")` matches multiply imputed datasets using several
    matching methods, and provides users with tools to estimate causal effects
    in each imputed dataset. `r pkg("grf")` offers treatment effect estimation
    with incomplete confounders and covariates under modified unconfoundedness
    assumptions and `r pkg("RCAL")` implements regularized calibrated estimation
    for causal inference with missing values and high dimension.
-   *Finance*: `r pkg("imputeFin")` handles imputation of missing values in
    financial time series using AR models or random walk.
-   *Finance*: Basic methods (mean, median, mode, \...) for imputing missing
    data in scoring datasets are proposed in `r pkg("scorecardModelUtils")`.
    `r pkg("creditmodel")` can handle missing values treatment for credit 
    modeling.
-   *Preference models*: Missing data in preference models are handled
    with a composite link approach that allows for MCAR and MNAR patterns to be
    accounted for in `r pkg("prefmod")`.
-   *Administrative records / Surveys*: `r pkg("BIFIEsurvey")` is a very general
    package that contains tools for survey statistics and that can handle
    multiply imputed datasets. More specifically, `r pkg("fastLink")` provides
    a Fellegi-Sunter probabilistic record linkage that allows for missing data
    and the inclusion of auxiliary information. `r pkg("eatRep")` implements 
    replication methods in complex survey designs comprising multiple imputed 
    variables, and `r pkg("modi")` provides multivariate outlier detection and
    `r pkg("convergEU")` can process data from Eurostat data and impute missing
    values to monitor convergence between EU countries. `r pkg("eechidna")` has
    similar features for Australian election and public census datasets.
-   *Bibliometry*: `r pkg("robustrao")` computes the Rao-Stirling diversity
    index (a well-established bibliometric indicator to measure the
    interdisciplinarity of scientific publications) with data containing
    uncategorized references. `r pkg("metagear")` provides hot-deck imputation
    in bibliographic data for systematic reviews and meta-analysis.
-   *Agriculture*: `r pkg("geneticae")` implements imputation techniques for
    multi-environment agronomic trials.



### Links
-   [Amelia II: A Program for Missing Data](https://gking.harvard.edu/amelia)
-   [R-miss-tastic: A resource website on missing data](https://rmisstastic.netlify.com/)
