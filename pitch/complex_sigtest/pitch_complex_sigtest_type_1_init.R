#-----------------------------
# code for loading / resampling data set 
#-----------------------------

# type of data: resampled pitch trajectories

# path previously determined via config file

dat <- readRDS(file.path(dirname(config.file.curr), data.file))

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
  } else if (r[1] == "gamcheck") {
    out <- paste(out, 
                 ' + s(measurement.no, ',  grouping, ', bs="fs", m=1, xt="',
                 r[2],
                 '", k=',
                 as.character(seq(as.numeric(r[3]),as.numeric(r[4]),as.numeric(r[5]))),
                 ')',
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
  fixef_formula <- "f0_log_norm ~ "
  nested_formula <- "f0_log_norm ~ "
  
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
    final_formula <- gsub("~ group.ordered", '~ s(group.ordered, bs="re")', final_formula)
    method_str <- 'method="fREML", select=TRUE'
  } else if (method == "discrete+select") {
    final_formula <- gsub("~ group.ordered", '~ s(group.ordered, bs="re")', final_formula)
    method_str <- 'method="fREML", discrete=T, nthreads=1, select=TRUE'
  }
  
  # assembling bam command
  bam_str <- paste("bam(", final_formula, ", data=dat, ", AR_str, method_str, ")", sep="")
  bam_nested_str <- paste("bam(", final_nested_formula, ", data=dat, ", AR_str, method_str, ")", sep="")
  return(list(full=bam_str, nested=bam_nested_str, grouping_var=grouping_var))
}

#------------------------------------------------
# Function for extracting p-value from gam.check
#------------------------------------------------

gam.check.p.value <- function (mod, which.line) {
  str.out <- capture.output(gam.check(mod))
  relevant.line <- str.out[grep(which.line, str.out)]
  p.value <- as.numeric(str_match(relevant.line, "([0-9.]*)[ *.]*$")[[2]])
  return(p.value)
}
