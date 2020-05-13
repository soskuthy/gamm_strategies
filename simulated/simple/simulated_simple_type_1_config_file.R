# name and paths
overall_name = "gamm_simulated_simple_type_1"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_strategies/simulated/simple"
path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/simulated/simple"
output.dir = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/simulated/simple/output_type_1"
r.script.path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "simulated_simple_type_1_init.r"

# name of config file
config.file = "simulated_simple_type_1_config_file.R"

# name of data file
data.file = NULL

# details of models to be run 
to_fit <- rbind(
                expand.grid(fixed_effects=c("diff_tp_10"),
                      random_effects=c("noranef","rintcpt","rslope",
                                       "rsmooth_tp_3","rsmooth_tp_5","rsmooth_tp_10",
                                       "rsmooth_cr_3","rsmooth_cr_5","rsmooth_cr_10",
                                       "gamcheck_tp_4_10_2", "gamcheck_cr_4_10_2"),
                      AR=c("noAR"),
                      method="discrete",
                      mod_comp="nomodcomp",
                      dataset=c("dense"),
                      visual="noVis"),
                expand.grid(fixed_effects=c("diff_tp_10"),
                            random_effects=c("noranef","rintcpt","rslope",
                                             "rsmooth_tp_3","rsmooth_tp_5","rsmooth_tp_10",
                                             "rsmooth_cr_3","rsmooth_cr_5","rsmooth_cr_10"),
                            AR=c("AR_est"),
                            method="discrete",
                            mod_comp="nomodcomp",
                            dataset=c("dense"),
                            visual="noVis")
)

