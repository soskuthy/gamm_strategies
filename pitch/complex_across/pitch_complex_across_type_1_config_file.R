# name and paths
overall_name = "gamm_pitch_real_complex_across_type_1"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_strategies/pitch/complex_across"
path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/pitch/complex_across"
output.dir = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/pitch/complex_across/output_type_1"
r.script.path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "pitch_complex_across_type_1_init.r"

# name of config file
config.file = "pitch_complex_across_type_1_config_file.R"

# name of data file
data.file = "pitch_contr_27_speakers.rds"

# details of models to be run 
to_fit <- rbind(expand.grid(fixed_effects=c("diff_tp_20"),
                            random_effects=c("noranef+noranef","noranef+rintcpt","noranef+rslope",
                                             "noranef+rsmooth_tp_3","noranef+rsmooth_tp_5","noranef+rsmooth_tp_10",
                                             "noranef+rsmooth_tp_15",
                                             "noranef+rsmooth_cr_3","noranef+rsmooth_cr_5","noranef+rsmooth_cr_10",
                                             "noranef+rsmooth_cr_15",
                                             "noranef+gamcheck_tp_4_16_3",
                                             "noranef+gamcheck_cr_4_16_3"),
                            AR=c("AR_est"),
                            method="discrete",
                            mod_comp="nomodcomp",
                            dataset=c("dense"),
                            visual="noVis"),
                expand.grid(fixed_effects=c("diff_tp_20"),
                            random_effects=c("noranef+rsmooth_tp_10","noranef+rsmooth_tp_15",
                                             "noranef+rsmooth_cr_10","noranef+rsmooth_cr_15",
                                             "noranef+gamcheck_tp_4_16_3",
                                             "noranef+gamcheck_cr_4_16_3"),
                            AR=c("noAR"),
                            method="discrete",
                            mod_comp="nomodcomp",
                            dataset=c("dense"),
                            visual="noVis")
)


