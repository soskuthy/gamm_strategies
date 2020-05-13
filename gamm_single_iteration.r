#-------------------------------
# initialisation: reading files
#-------------------------------

# no warnings

options(warn=-1)

# read in command arguments
all.args <- commandArgs()
args <- commandArgs(trailingOnly=TRUE)

# find path to script (*only works when script is run from 
# command line or using source()*)

script_name <- sub("--file=", "", all.args[grep("--file=", all.args)])
script_dir <- dirname(script_name)


# first argument should be config file
config.file.curr <- args[1]


# second argument should be iteration number
iteration <- as.numeric(args[2])

#script_dir <- dirname(sys.frame(1)$ofile)

# load libraries
t <- system.time(source(file.path(script_dir, "libraries.r")))
cat("Iteration ", iteration, "; libraries loaded in ", t[[3]], " seconds.\n", sep="")

# source config file
# this provides the following variables:
#   - name: name of simulation
#   - output_dir: name of directory with output files (script assumes this already exists)
#   - init_file: name of init file
#   PLUS either
#     - fixed_effects, random_effects, AR, method
#     - mod_comp
#        OR
#     - to_fit (data frame specifying details of each model separately)

t <- system.time(source(config.file.curr))
cat("Iteration ", iteration, "; config file loaded in ", t[[3]], " seconds.\n", sep="")

# source init file (assumed to be in same directory as job file)
t <- system.time(source(file.path(dirname(config.file.curr), init.file)))
cat("Iteration ", iteration, "; data set created in ", t[[3]], " seconds.\n", sep="")

# at this point, we have the data and a function for specifying bam()'s from
# string parameters

#------------------------------------
# create list of models to be fitted
#------------------------------------

# if to_fit is not provided, the elements of fixed_effects, 
# random_effects, AR, method and modcomp are crossed
if (!exists("to_fit")) {
  to_fit <- expand.grid(fixed_effects=fixed_effects,
                        random_effects=random_effects,
                        AR=AR,
                        method=method,
                        mod_comp=mod_comp)
}

# columns of to_fit should be strings

for (cl in 1:ncol(to_fit)) {
  to_fit[,cl] <- as.character(to_fit[,cl])
}

#-------------
# main loop
#-------------

# create output list

output <- list()

for (r in 1:nrow(to_fit)) {
  pars <- to_fit[r,]
  # simulation name (to be used as ID in output list)
  sim_name <- paste(pars, collapse=".")
  
  cat("fitting", sim_name, iteration, "\n")
  
  # generate code for running bams
  bams_code <- assemble_bam(pars$fixed_effects, 
                            pars$random_effects, 
                            pars$AR, 
                            pars$method,
                            pars$dataset)
  full_code <- parse(text=paste("full_mod <-", bams_code$full))
  nested_code <- parse(text=paste("nested_mod <-", bams_code$nested))
  grouping_var <- bams_code$grouping_var
  
  # we start collecting timing data here
  t <- system.time({
    
    # if AR, then first fit model with fREML to get rho.est; 
    if (pars$AR == "AR_est") {
      if (substring(pars$random_effects, 1, nchar("gamcheck")) == "gamcheck") {
        stop("Combining random effect string 'gamcheck' with AR models not yet implemented")
      }
      bam_AR_code <- parse(text=paste("AR_mod <-", assemble_bam(pars$fixed_effects, 
                                                                pars$random_effects, 
                                                                "noAR", 
                                                                "discrete",
                                                                pars$dataset)$full
                                      )
      )
      eval(bam_AR_code)
      rho.est <- start_value_rho(AR_mod)
    }
    
    # if gamcheck, fit increasingly complex models and run gam.check at each
    # point to see if there is undersmoothing; break loop if model is
    # sufficiently complex
    
    if (grepl("gamcheck", pars$random_effects)) {
      if (grepl("[+]gamcheck", pars$random_effects)) {
        random_groups <- "speaker"
      } else {
        random_groups <- "traj"
      }
      for (counter in 1:length(full_code)) {
        curr_full_code <- full_code[counter]

        # fit full model, record peak memory consumption
        
        zz <- textConnection("memory_use", "w")
        sink(zz, type="message")
        gcinfo(T)
        gc(verbose=T)
        gc(verbose=T)
        eval(curr_full_code)
        gc(verbose=T)
        gcinfo(F)
        sink(file=NULL, type="message")
        close(zz)
        
        mems <- as.numeric(str_extract(memory_use[grep("Mbytes of vectors", memory_use)], "^[0-9]*[.]*[0-9]*"))
        mem.peak <- max(mems) - mems[1]
        
        # break loop if model is sufficiently complex or maximum complexity
        # has been reached
        
        if (gam.check.p.value(full_mod, random_groups) >= 0.05 | counter == length(full_code)) {
          full_code <- curr_full_code
          bams_code$full <- bams_code$full[counter]
          nested_code <- nested_code[counter]
          bams_code$nested <- bams_code$nested[counter]
          break
        }
        
        # (note: k for final model can be read out of full_code
        # slot of output object - so not recorded separately)
      }
    } else {
    
    # if not gamcheck:
    
      # fit full model, record peak memory consumption
      
      zz <- textConnection("memory_use", "w")
      sink(zz, type="message")
      gcinfo(T)
      gc(verbose=T)
      gc(verbose=T)
      eval(full_code)
      gc(verbose=T)
      gcinfo(F)
      sink(file=NULL, type="message")
      close(zz)
      
      mems <- as.numeric(str_extract(memory_use[grep("Mbytes of vectors", memory_use)], "^[0-9]*[.]*[0-9]*"))
      mem.peak <- max(mems) - mems[1]
    
    }
    
    # get model summary
    
    full_summary <- summary(full_mod)
    out_summary <- list(smooth=full_summary$s.table,
                        parametric=full_summary$p.table,
                        ml=full_summary$sp.criterion)
    
    # get visual data: 
    # (1) preds & confints for both smooths
    # (2) preds & confints for difference smooth
    
    fixef_specs <- str_split(pars$fixed_effects, "_")[[1]]
    if (pars$visual == "vis") {
      pdf(file=NULL)
      
      # overlap between individual smooths
      cond_list_1 <- list()
      cond_list_2 <- list()
      if (fixef_specs[1]=="diff") {
        cond_list_1[[grouping_var]] <- "A"
        cond_list_2[[grouping_var]] <- "B"
      } else if (fixef_specs[1]=="bin") {
        cond_list_1[[grouping_var]] <- 0
        cond_list_2[[grouping_var]] <- 1
      }
      d1 <- plot_smooth(full_mod, view="measurement.no", cond=cond_list_1, 
                        print.summary=F, rm.ranef=T, n.grid=100, rug=F)$fv
      d2 <- plot_smooth(full_mod, view="measurement.no", cond=cond_list_2, 
                        print.summary=F, rm.ranef=T, n.grid=100, rug=F)$fv
      
      visual_no_overlap <- mean( (d1$ul < d2$ll) | (d1$ll > d2$ul) )
      
      # difference smooth
      comp_list <- list()
      comp_list[[grouping_var]] <- c(cond_list_2[[1]], cond_list_1[[1]])
      diff_confint <- plot_diff(full_mod, view="measurement.no", comp=comp_list, 
                               print.summary=F, rm.ranef=T, n.grid=100)
      visual_excludes_0 <- mean( !((0 > diff_confint$est - diff_confint$CI) & 
                                   (0 < diff_confint$est + diff_confint$CI)) )
      
      dev.off()
    } else {
      visual_no_overlap <- NULL
      visual_excludes_0 <- NULL
    }
    
    # and if there is model comparison, fit nested model & run anova
    
    model_comparison <- NULL
    
    if (pars$mod_comp == "modcomp") {
      eval(nested_code)
      model_comparison <- compareML(full_mod, nested_mod, print.output=F)
      model_comparison$AIC_2 <- AIC(full_mod, nested_mod)
    }
    
    # generate output
    
    
  })
  output[[sim_name]] <- list(parameters=pars,
                             full_code=bams_code$full,
                             nested_code=bams_code$nested,
                             summary=out_summary,
                             visual_no_overlap=visual_no_overlap,
                             visual_excludes_0=visual_excludes_0,
                             model_comparison=model_comparison,
                             time=t,
                             memory=mem.peak,
                             speakers=as.character(unique(dat$speaker)))
}

#-------------
# save output
#-------------

saveRDS(output, file=file.path(output.dir, paste(overall_name, "_", iteration, ".rds", sep="")))
