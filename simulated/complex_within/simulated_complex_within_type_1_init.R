#-----------------------------
# code for creating data set 
#-----------------------------

# type of data: simulated f2 trajectories for /eI/ (logistic curve)
# 2000 trajectories, 50 words
# parameterised variation in
#   - start of trajectory
#   - end of trajectory
#   - transition point
#   - transition steepness
# + a bit of random noise
#   - hierarchical structure (measurements within trajectories within speakers)
#   - across-speaker variation

# setting time dimension

xs = seq(0,1,0.1)

# population parameters: individual speakers come from this dist
f2_start_mean = 1300
f2_end_mean = 1650
f2_start_sd.speaker = 90
f2_end_sd.speaker = 90
# correlation between pairs of parameter values sampled for two groups within same speaker;
# rho = 7/8 means that the sd of the difference between the two groups within
#   the same speaker is 0.5 * the sd of the target values
# in mixed model lingo: the sd of the random slope is half of the sd of the intercept
f2_start_rho.speaker = 7/8
f2_end_rho.speaker = 7/8

# expected value & sd for transition point
x0_mean = 0.6
x0_sd.speaker = 0.020
x0_rho.speaker = 7/8 # same as above
# expected value & sd for steepness (higher -> more steep)
k_mean = 15
k_sd.speaker = 4
k_rho.speaker = 7/8 # same as above

# how much variation within speaker-group pairs? (unchanged from before)
f2_start_sd.traj = 150
f2_end_sd.traj = 150
x0_sd.traj = 0.015
k_sd.traj = 3

# amount of random noise

noise_sd <- 40

n_speakers <- 50
n_trajectories_per_speaker <- 40

# assembling trajectories

make_cov_matrix <- function (a.sd, b.sd, rho) {
  matrix(c(a.sd**2, rho*a.sd*b.sd, rho*a.sd*b.sd, b.sd**2), nrow=2)
}

ys_m <- matrix(0, nrow=length(xs), ncol=n_speakers*n_trajectories_per_speaker)
for (i in 1:n_speakers) {
  f2_start.speaker <- mvrnorm(1, rep(f2_start_mean, 2), 
                           make_cov_matrix(f2_start_sd.speaker, f2_start_sd.speaker, f2_start_rho.speaker))
  f2_end.speaker <- mvrnorm(1, rep(f2_end_mean, 2), 
                         make_cov_matrix(f2_end_sd.speaker, f2_end_sd.speaker, f2_end_rho.speaker))
  x0.speaker <- mvrnorm(1, rep(x0_mean, 2), 
                     make_cov_matrix(x0_sd.speaker, x0_sd.speaker, x0_rho.speaker))
  k.speaker <- mvrnorm(1, rep(k_mean, 2), 
                     make_cov_matrix(k_sd.speaker, k_sd.speaker, k_rho.speaker))
  for (j in 1:(n_trajectories_per_speaker/2)) {
    # group A
    f2_start <- rnorm(1, f2_start.speaker[1], f2_start_sd.traj)
    f2_end <- rnorm(1, f2_end.speaker[1], f2_end_sd.traj)
    x0 <- rnorm(1, x0.speaker[1], x0_sd.traj)
    k <- rnorm(1, k.speaker[1], k_sd.traj)
    ys_m[,(i-1)*n_trajectories_per_speaker + j*2 - 1] <- ((f2_end - f2_start) / (1 + exp(-k*(xs-x0)))) + f2_start + rnorm(length(xs), 0, noise_sd)
    # group B
    f2_start <- rnorm(1, f2_start.speaker[2], f2_start_sd.traj)
    f2_end <- rnorm(1, f2_end.speaker[2], f2_end_sd.traj)
    x0 <- rnorm(1, x0.speaker[2], x0_sd.traj)
    k <- rnorm(1, k.speaker[2], k_sd.traj)
    ys_m[,(i-1)*n_trajectories_per_speaker + j*2] <- ((f2_end - f2_start) / (1 + exp(-k*(xs-x0)))) + f2_start + rnorm(length(xs), 0, noise_sd)
  }
}

# assembling data set (randomly assigned to categories)
dat <- data.frame(traj=paste("traj_", rep(1:(n_speakers*n_trajectories_per_speaker), each=length(xs)), sep=""),
                  speaker=paste("speaker_", rep(1:n_speakers, each=length(xs)*n_trajectories_per_speaker), sep=""),
                  group=rep(c("A","B"), each=length(xs), times=n_speakers*n_trajectories_per_speaker / 2),
                  measurement.no=xs, 
                  f2=c(ys_m),
                  stringsAsFactors = F
                 )

# setting up different types of grouping factors
dat$group.factor <- as.factor(dat$group)
dat$group.ordered <- as.ordered(dat$group)
contrasts(dat$group.ordered) <- "contr.treatment"
dat$group.bin <- as.numeric(dat$group.factor) - 1

# ids ought to be factors  
dat$traj <- as.factor(dat$traj)
dat$speaker <- as.factor(dat$speaker)
dat$speakerGroup <- interaction(dat$speaker, dat$group)

# add dat$start for AR.start (for autoregressive error models)

dat$start <- dat$measurement.no == 0

#ggplot(dat, aes(x=measurement.no, y=f2, group=traj, col=group)) + geom_line(alpha=0.5) + facet_wrap(~ speaker)
#ggplot(dat, aes(x=measurement.no, y=f2, col=group)) + geom_smooth(alpha=0.5) + facet_wrap(~ speaker)

#-----------------------------
# function for assembling bam
#-----------------------------

### utility function for assembling random structure
# INPUT: rintcpt OR rslope OR rsmooth_(tp/cr/..)_(k), e.g. rsmooth_cr_3
#        name of grouping variable
assemble_rstruct <- function (r, grouping) {
  out <- ''
  if (r[1] == "rintcpt") {
    out <- paste(out, ' + s(',  grouping, ', bs="re")', sep="")
  } else if (r[1] == "rslope") {
    out <- paste(out, ' + s(',  grouping, ', bs="re") + s(',  grouping, ', measurement.no, bs="re")', sep="")
  } else if (r[1] == "rsmooth") {
    out <- paste(out, 
                 ' + s(measurement.no, ',  grouping, ', bs="fs", m=1, xt="',
                 r[2],
                 '", k=',
                 r[3],
                 ')',
                 sep="")
  } else if (r[1] == "rsmoothslope") {
    out <- paste(out, 
                 ' + s(measurement.no, ',  grouping[1], ', bs="fs", m=1, xt="',
                 r[2],
                 '", k=',
                 r[3],
                 ') + s(',
                 grouping[1],
                 ', group.ordered, bs="re")',
                 sep="")
  } else if (r[1] == "rsmoothcrossed") {
    out <- paste(out, 
                 ' + s(measurement.no, ',  grouping[1], 'Group, bs="fs", m=1, xt="',
                 r[2],
                 '", k=',
                 r[3],
                 ')',
                 sep="")
  } else if (r[1] == "rsmoothdiff") {
    out <- paste(out, 
                 ' + s(measurement.no, ', grouping[1], ', bs="fs", m=1, xt="',
                 r[2],
                 '", k=',
                 r[3],
                 ') + s(measurement.no, ', grouping[1], ', bs="fs", m=1, xt="',
                 r[2],
                 '", k=',
                 r[3],
                 ', by=group.ordered)',
                 sep="")
  } else if (r[1] == "rsmoothsep") {
    out <- paste(out, 
                 ' + s(measurement.no, ', grouping[1], ', bs="fs", m=1, xt="',
                 r[2],
                 '", k=',
                 r[3],
                 ', by=group.factor)',
                 sep="")
  } else if (r[1] == "noranef") {
  } else {
    stop("Invalid random effect string: has to start with one of 'rintcpt', 'rslope', 'rsmooth' or 'noranef'")
  }
  return(out)
}

assemble_bam <- function (fixefs, ranefs, AR, method, dataset) {
  
  # extract parameters from fixef and ranef string specification
  fixef_specs <- str_split(fixefs, "_")[[1]]
  ranef_specs <- str_split(unlist(str_split(ranefs, "[+]")[[1]]), "_")
  
  # initialise fixed effect part of formula
  fixef_formula <- "f2 ~ "
  nested_formula <- "f2 ~ "
  
  if (fixef_specs[1] == "diff") {
    # name of grouping variable (part of output)
    grouping_var <- "group.ordered"
    
    # if a difference smooth is used, assemble relevant fixef formula
    # using bs and k specifications from fixefs
    fixef_formula <- paste(fixef_formula,
                           'group.ordered + s(measurement.no, bs="',
                           fixef_specs[2],
                           '", k=',
                           fixef_specs[3],
                           ') + s(measurement.no, by=group.ordered, bs="',
                           fixef_specs[2],
                           '", k=',
                           fixef_specs[3],
                           ')',
                           sep="")
    nested_formula <- paste(nested_formula,
                            's(measurement.no, bs="',
                            fixef_specs[2],
                            '", k=',
                            fixef_specs[3],
                            ')',
                            sep="")
  } else if (fixef_specs[1] == "bin") {
    # name of grouping variable (part of output)
    grouping_var <- "group.bin"
    
    # if a binary smooth is used, assemble relevant fixef formula
    # using bs and k specifications from fixefs (basically, no parametric term)
    fixef_formula <- paste(fixef_formula,
                           's(measurement.no, bs="',
                           fixef_specs[2],
                           '", k=',
                           fixef_specs[3],
                           ') + s(measurement.no, by=group.bin, bs="',
                           fixef_specs[2],
                           '", k=',
                           fixef_specs[3],
                           ')',
                           sep="")
    nested_formula <- paste(nested_formula,
                            's(measurement.no, bs="',
                            fixef_specs[2],
                            '", k=',
                            fixef_specs[3],
                            ')',
                            sep="")
  } else {
    stop("Invalid 'by' specification: has to be either 'diff' or 'bin'.")
  }
  
  # code for generating formula for random intercept / random intercept + slope
  # and random smooth (with user-specified k and basis type)
  
  ranef_formula <- ''
  groupers <- c("traj", "speaker")
  for (j in 1:length(ranef_specs)) {
    ranef_formula <- paste(ranef_formula, assemble_rstruct(ranef_specs[[j]], groupers[j]), sep="")
  }
  
  final_formula <- paste(fixef_formula, ranef_formula, sep="")
  final_nested_formula <- paste(nested_formula, ranef_formula, sep="")
  
  # setting AR parameters based on string from job file
  
  AR_str <- ""
  if (AR == "AR_est") {
    AR_str <- paste("AR.start=dat$start, rho=rho.est, ", sep="")
  } else if (AR != "noAR") {
    AR_str <- paste("AR.start=dat$start, rho=", str_split(AR, "_")[[1]][2], sep="")
  }
  
  # setting method parameter(s) based on string from job file
  
  method_str <- ""
  if (method == "ML") {
    method_str <- 'method="ML"'
  } else if (method == "REML") {
    method_str <- 'method="REML"'
  } else if (method == "fREML") {
    method_str <- 'method="fREML"'
  } else if (method == "discrete") {
    method_str <- 'method="fREML", discrete=T, nthreads=1'
  } else if (method == "fREML+select") {
    method_str <- 'method="fREML", select=TRUE'
  } else if (method == "discrete+select") {
    method_str <- 'method="fREML", discrete=T, nthreads=1, select=TRUE'
  }
  
  # assembling bam command
  bam_str <- paste("bam(", final_formula, ", data=dat, ", AR_str, method_str, ")", sep="")
  bam_nested_str <- paste("bam(", final_nested_formula, ", data=dat, ", AR_str, method_str, ")", sep="")
  return(list(full=bam_str, nested=bam_nested_str, grouping_var=grouping_var))
}
