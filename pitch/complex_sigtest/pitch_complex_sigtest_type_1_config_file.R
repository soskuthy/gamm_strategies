# name and paths
overall_name = "gamm_pitch_real_complex_sigtest_type_1"
path.on.local = "~/documents/research/projects/dynamic-gam/gamm_strategies/pitch/complex_sigtest"
path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/pitch/complex_sigtest"
output.dir = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/pitch/complex_sigtest/output_type_1"
r.script.path.on.server = "/home/soskuthy/projects/def-soskuthy/soskuthy/gamm_strategies/gamm_single_iteration.r"

# name of init file (containing details of data & models)
init.file = "pitch_complex_sigtest_type_1_init.r"

# name of config file
config.file = "pitch_complex_sigtest_type_1_config_file.R"

# name of data file
data.file = "pitch_contr_27_speakers.rds"

# - fREML par/smooth/both + modcomp + vis
# - ML par /smooth/both + modcomp + vis
# - fREML bin smooth
# - fREML par/smooth/both select


# details of models to be run 
to_fit <- rbind(expand.grid(fixed_effects=c("diff_tp_20"),
                            random_effects=c("noranef+rsmooth_cr_10"),
                            AR=c("noAR"),
                            method=c("discrete","ML"),
                            mod_comp="modcomp",
                            dataset=c("dense"),
                            visual="vis"),
                expand.grid(fixed_effects=c("bin_tp_20"),
                            random_effects=c("noranef+rsmooth_cr_10"),
                            AR=c("noAR"),
                            method=c("discrete"),
                            mod_comp="nomodcomp",
                            dataset=c("dense"),
                            visual="noVis"),
                expand.grid(fixed_effects=c("diff_tp_10"),
                            random_effects=c("noranef+rsmooth_cr_10"),
                            AR=c("noAR"),
                            method=c("discrete+select"),
                            mod_comp="modcomp",
                            dataset=c("dense"),
                            visual="noVis")
)


