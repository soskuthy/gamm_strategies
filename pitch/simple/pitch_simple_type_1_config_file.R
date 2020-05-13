# name and paths
overall_name = "gamm_pitch_real_simple_type_1"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_strategies/pitch/simple"
path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/pitch/simple"
output.dir = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/pitch/simple/output_type_1"
r.script.path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "pitch_simple_type_1_init.r"

# name of config file
config.file = "pitch_simple_type_1_config_file.R"

# name of data file
data.file = "pitch_contr_27_speakers.rds"

# details of models to be run
to_fit <- rbind(
                expand.grid(fixed_effects=c("diff_tp_20"),
                      random_effects=c("noranef","rintcpt","rslope",
                                       "rsmooth_tp_3","rsmooth_tp_5","rsmooth_tp_10",
                                       "rsmooth_tp_15","rsmooth_tp_20",
                                       "rsmooth_cr_3","rsmooth_cr_5","rsmooth_cr_10","rsmooth_cr_15",
                                       "rsmooth_cr_15","rsmooth_cr_20",
                                       "gamcheck_tp_4_16_3", "gamcheck_cr_4_16_3"),
                      AR=c("noAR"),
                      method="discrete",
                      mod_comp="nomodcomp",
                      dataset=c("dense"),
                      visual="noVis"),
                expand.grid(fixed_effects=c("diff_tp_20"),
                            random_effects=c("noranef","rintcpt","rslope",
                                             "rsmooth_tp_3","rsmooth_tp_5","rsmooth_tp_10",
                                             "rsmooth_tp_15","rsmooth_tp_20",
                                             "rsmooth_cr_3","rsmooth_cr_5","rsmooth_cr_10","rsmooth_cr_15",
                                             "rsmooth_cr_15","rsmooth_cr_20"),
                            AR=c("AR_est"),
                            method="discrete",
                            mod_comp="nomodcomp",
                            dataset=c("dense"),
                            visual="noVis")
)
