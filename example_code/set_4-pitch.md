Introduction
------------

This document is a supplement to “Evaluating generalised additive mixed
modelling strategies for dynamic speech analysis,” relating specifically
to the contents of the “real pitch” columns of Table 4 in Section 3.4.2.
It presents code that illustrates (i) how the resampled data were
generated and (ii) the models whose performance is summarised in the
“real pitch” columns of Table 4.

Preliminaries
-------------

The code below loads the relevant libraries.

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 3.5.2

``` r
library(mgcv)
```

    ## Warning: package 'mgcv' was built under R version 3.5.2

``` r
library(itsadug)
library(stringr)
```

    ## Warning: package 'stringr' was built under R version 3.5.2

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 3.5.2

Data generation
---------------

The code in this section can be only be used to process the data for
type I simulations. Note that the paths in this file will only work if
the whole GitHub repository is downloaded and this markdown file is kept
in its folder.

The data for this set of simulations consist of real pitch trajectories
representing contrastive focus in Standard German. All contours (max.
20) are used from each of the 27 speakers in the data set. Each speaker
is randomly assigned to one of two groups (A and B). The difference
between males and females has been residualised out of the f0
measurements, which are also logged (but not normalised otherwise).

``` r
dat <- readRDS("../data/final_data/pitch_contr_27_speakers.rds")

# we now add randomly assigned category labels

ids <- unique(dat$speaker)
group.Bs <- sample(ids, round(length(ids)/2))
dat$group <- "A"
dat$group[dat$speaker %in% group.Bs] <- "B"

# setting up different types of grouping factors
dat$group.factor <- as.factor(dat$group)
dat$group.ordered <- as.ordered(dat$group)
contrasts(dat$group.ordered) <- "contr.treatment"
dat$group.bin <- as.numeric(dat$group.factor) - 1

# ids ought to be factors  
dat$traj <- as.factor(dat$traj)
dat$speaker <- as.factor(dat$speaker)

# dat$start has already been added at data prep stage (for AR.start, i.e. for autoregressive error models)
```

Here is what the data set looks like.

``` r
ggplot(dat, aes(x=measurement.no, y=f0_log_norm, group=traj, col=group)) +
  geom_line(alpha=0.6) +
  facet_wrap(~speaker) +
  theme_bw()
```

![](set_4-pitch_files/figure-markdown_github/unnamed-chunk-3-1.png)

Methods of significance testing
-------------------------------

All the models (and sets of models) from Table 4 are shown below in the
same order as in the table. The numbers in the section headers
correspond to the row numbers. Note that all models contain AR1
components to deal with dependencies within trajectories. The rho value
used for these AR1 components is taken from a single model fitted
without any random structures. This model is estimated below.

``` r
# thin
rho_mod <- bam(f0_log_norm ~ group.ordered + 
                     s(measurement.no, bs = "tp", k = 20) + 
                     s(measurement.no, by = group.ordered, bs = "tp", k = 20), 
                   data = dat, method = "fREML", discrete = T, nthreads = 1)

rho <- start_value_rho(rho_mod)
```

### MODEL SUMMARY: 1. Looking at model summary

``` r
modsum <- bam(f0_log_norm ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "fREML", discrete = T, nthreads = 1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
summary(modsum)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f0_log_norm ~ group.ordered + s(measurement.no, bs = "tp", k = 20) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 20) + 
    ##     s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10)
    ## 
    ## Parametric coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)    0.0081691  0.0361880   0.226    0.821
    ## group.orderedB 0.0004903  0.0501399   0.010    0.992
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf  Ref.df      F p-value    
    ## s(measurement.no)                 16.271  17.697  6.831  <2e-16 ***
    ## s(measurement.no):group.orderedB   6.011   7.861  0.455   0.881    
    ## s(measurement.no,speaker)        206.942 266.000 14.761  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.803   Deviance explained = 80.4%
    ## fREML = -38493  Scale est. = 0.016669  n = 20820

### MODEL SUMMARY: 2. Looking at model summary + Bonferroni correction

For the Bonferroni correction, the alpha-level of the parametric and
smooth terms is lowered to 0.025. This does not require fitting a
separate model.

### MODEL SUMMARY: 3. Binary smooth

``` r
binsmooth <- bam(f0_log_norm ~ s(measurement.no, bs = "tp", k = 20) + 
                      s(measurement.no, by = group.bin, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "fREML", discrete=T, nthreads=1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
summary(binsmooth)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## f0_log_norm ~ s(measurement.no, bs = "tp", k = 20) + s(measurement.no, 
    ##     by = group.bin, bs = "tp", k = 20) + s(measurement.no, speaker, 
    ##     bs = "fs", m = 1, xt = "cr", k = 10)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept) 0.008169   0.036188   0.226    0.821
    ## 
    ## Approximate significance of smooth terms:
    ##                                 edf  Ref.df      F p-value    
    ## s(measurement.no)            16.271  17.697  6.831  <2e-16 ***
    ## s(measurement.no):group.bin   7.011   8.861  0.462   0.914    
    ## s(measurement.no,speaker)   206.942 267.000 14.419  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.803   Deviance explained = 80.4%
    ## fREML = -38493  Scale est. = 0.016669  n = 20820

### LIKELIHOOD RATIO TESTS: 4. Likelihood Ratio Test using models fitted with ML

Please note that these models may take quite a while to fit (5-10
minutes).

``` r
lrt_ML_full <- bam(f0_log_norm ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "ML")
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
lrt_ML_nested <- bam(f0_log_norm ~ # group.ordered + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      # s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "ML")
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
compareML(lrt_ML_full, lrt_ML_nested)
```

    ## lrt_ML_full: f0_log_norm ~ group.ordered + s(measurement.no, bs = "tp", k = 20) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 20) + 
    ##     s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10)
    ## 
    ## lrt_ML_nested: f0_log_norm ~ s(measurement.no, bs = "tp", k = 20) + s(measurement.no, 
    ##     speaker, bs = "fs", m = 1, xt = "cr", k = 10)
    ## 
    ## Chi-square test of ML scores
    ## -----
    ##           Model     Score Edf Difference    Df p.value Sig.
    ## 1 lrt_ML_nested -38490.44   6                              
    ## 2   lrt_ML_full -38491.46   9      1.024 3.000   0.562     
    ## 
    ## AIC difference: 0.83, model lrt_ML_nested has lower AIC.

    ## Warning in compareML(lrt_ML_full, lrt_ML_nested): AIC might not be
    ## reliable, as an AR1 model is included (rho1 = 0.959661, rho2 = 0.959661).

    ## Warning in compareML(lrt_ML_full, lrt_ML_nested): Only small difference in ML...

### LIKELIHOOD RATIO TESTS: 5. Likelihood Ratio Test using models fitted with fREML

As noted in the main text of the paper, the results of this model
comparison are meaningless.

``` r
lrt_fREML_full <- bam(f0_log_norm ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "fREML", discrete=T, nthreads=1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
lrt_fREML_nested <- bam(f0_log_norm ~ # group.ordered + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      # s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "fREML", discrete=T, nthreads=1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
compareML(lrt_fREML_full, lrt_fREML_nested)
```

    ## lrt_fREML_full: f0_log_norm ~ group.ordered + s(measurement.no, bs = "tp", k = 20) + 
    ##     s(measurement.no, by = group.ordered, bs = "tp", k = 20) + 
    ##     s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10)
    ## 
    ## lrt_fREML_nested: f0_log_norm ~ s(measurement.no, bs = "tp", k = 20) + s(measurement.no, 
    ##     speaker, bs = "fs", m = 1, xt = "cr", k = 10)
    ## 
    ## Model lrt_fREML_nested preferred: lower fREML score (3.959), and lower df (3.000).
    ## -----
    ##              Model     Score Edf Difference    Df
    ## 1   lrt_fREML_full -38492.60   9                 
    ## 2 lrt_fREML_nested -38496.56   6     -3.959 3.000
    ## 
    ## AIC difference: -1.97, model lrt_fREML_full has lower AIC.

    ## Warning in compareML(lrt_fREML_full, lrt_fREML_nested): AIC might not be
    ## reliable, as an AR1 model is included (rho1 = 0.959661, rho2 = 0.959661).

    ## Warning in compareML(lrt_fREML_full, lrt_fREML_nested): Only small difference in fREML...

### LIKELIHOOD RATIO TESTS: 6. Likelihood Ratio Test with fREML trick

For this model, fixed effects are estimated as random effects, which
makes model comparison based on models fitted with (f)REML valid in
principle.

``` r
lrt_fREML_trick_full <- bam(f0_log_norm ~ s(group.ordered, bs="re") + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "fREML", discrete=T, nthreads=1,
                    select=T)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
lrt_fREML_trick_nested <- bam(f0_log_norm ~ # s(group.ordered, bs="re") + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      # s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "fREML", discrete=T, nthreads=1,
                    select=T)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
compareML(lrt_fREML_trick_full, lrt_fREML_trick_nested)
```

    ## lrt_fREML_trick_full: f0_log_norm ~ s(group.ordered, bs = "re") + s(measurement.no, 
    ##     bs = "tp", k = 20) + s(measurement.no, by = group.ordered, 
    ##     bs = "tp", k = 20) + s(measurement.no, speaker, bs = "fs", 
    ##     m = 1, xt = "cr", k = 10)
    ## 
    ## lrt_fREML_trick_nested: f0_log_norm ~ s(measurement.no, bs = "tp", k = 20) + s(measurement.no, 
    ##     speaker, bs = "fs", m = 1, xt = "cr", k = 10)
    ## 
    ## Chi-square test of fREML scores
    ## -----
    ##                    Model     Score Edf Difference    Df p.value Sig.
    ## 1 lrt_fREML_trick_nested -38497.72   6                              
    ## 2   lrt_fREML_trick_full -38498.31   9      0.590 3.000   0.758     
    ## 
    ## AIC difference: -0.35, model lrt_fREML_trick_full has lower AIC.

    ## Warning in compareML(lrt_fREML_trick_full, lrt_fREML_trick_nested): AIC
    ## might not be reliable, as an AR1 model is included (rho1 = 0.959661, rho2 =
    ## 0.959661).

    ## Warning in compareML(lrt_fREML_trick_full, lrt_fREML_trick_nested): Only small difference in fREML...

### Model comparison with AIC (models fitted with fREML)

``` r
aic_fREML_full <- bam(f0_log_norm ~ group.ordered + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "fREML", discrete=T, nthreads=1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
aic_fREML_nested <- bam(f0_log_norm ~ # group.ordered + 
                      s(measurement.no, bs = "tp", k = 20) + 
                      # s(measurement.no, by = group.ordered, bs = "tp", k = 20) +
                      s(measurement.no, speaker, bs = "fs", m = 1, xt = "cr", k = 10), 
                    data = dat, 
                    AR.start = dat$start, rho = rho, 
                    method = "fREML", discrete=T, nthreads=1)
```

    ## Warning in gam.side(sm, X, tol = .Machine$double.eps^0.5): model has
    ## repeated 1-d smooths of same variable.

``` r
AIC(aic_fREML_full, aic_fREML_nested)
```

    ##                        df       AIC
    ## aic_fREML_full   234.2038 -77413.29
    ## aic_fREML_nested 230.5276 -77411.32

### VISUAL TESTS

This is not a detailed implementation of the percentage-cut-off-based
reasoning examined in the paper, simply some example code that can be
used to generate (i) prediction plots with confidence intervals for the
two groups and (ii) plots of the estimated difference between the
groups. Note also that these pred

``` r
plot_smooth(modsum, view="measurement.no", plot_all="group.ordered", rm.ranef=T)
```

    ## Summary:
    ##  * group.ordered : factor; set to the value(s): A, B. 
    ##  * measurement.no : numeric predictor; with 30 values ranging from 0.000000 to 1.000000. 
    ##  * speaker : factor; set to the value(s): BP_91. (Might be canceled as random effect, check below.) 
    ##  * NOTE : The following random effects columns are canceled: s(measurement.no,speaker)
    ## 

![](set_4-pitch_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
plot_diff(modsum, view="measurement.no", comp=list(group.ordered=c("A","B")), rm.ranef=T)
```

    ## Summary:
    ##  * measurement.no : numeric predictor; with 100 values ranging from 0.000000 to 1.000000. 
    ##  * speaker : factor; set to the value(s): BP_91. (Might be canceled as random effect, check below.) 
    ##  * NOTE : The following random effects columns are canceled: s(measurement.no,speaker)
    ## 

![](set_4-pitch_files/figure-markdown_github/unnamed-chunk-11-2.png)

    ## 
    ## Difference is not significant.
