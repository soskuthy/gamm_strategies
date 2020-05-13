# name and paths
overall_name = "gamm_simulated_complex_within_type_1"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_strategies/simulated/complex_within"
path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/simulated/complex_within"
output.dir = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/simulated/complex_within/output_type_1"
r.script.path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "simulated_complex_within_type_1_init.r"

# name of config file
config.file = "simulated_complex_within_type_1_config_file.R"

# name of data file
data.file = NULL

# details of models to be run 
to_fit <- expand.grid(fixed_effects=c("diff_tp_10"),
                      random_effects=c("noranef+rsmooth_cr_3","noranef+rsmooth_cr_5","noranef+rsmooth_cr_10",
                                       "noranef+rsmoothslope_cr_3","noranef+rsmoothslope_cr_5","noranef+rsmoothslope_cr_10",
                                       "noranef+rsmoothcrossed_cr_3","noranef+rsmoothcrossed_cr_5","noranef+rsmoothcrossed_cr_10",
                                       "noranef+rsmoothsep_cr_3","noranef+rsmoothsep_cr_5","noranef+rsmoothsep_cr_10",
                                       "noranef+rsmoothdiff_cr_3","noranef+rsmoothdiff_cr_5","noranef+rsmoothdiff_cr_10"),
                      AR=c("noAR", "AR_est"),
                      method="discrete",
                      mod_comp="nomodcomp",
                      dataset=c("dense"),
                      visual="noVis")


