library(tidyverse)
library(mgcv)
library(itsadug)
library(gridExtra)
library(lme4)
library(optimx)


#goose <- readRDS("~/Documents/Research/projects/north_fronting/topics/goose_fin_topics.rdata")
price <- readRDS("~/Documents/Research/projects/dynamic-gam/gamm_strategies/data/raw_data/price_full.rds")
price_spkrs <- read_delim("~/Documents/Research/projects/dynamic-gam/gamm_strategies/data/raw_data/speakers_to_use.csv", delim="\t", col_names=F)
price_to_save <- price %>%
  filter(following_voiceless==F,
         speaker %in% price_spkrs$X2,
         dur > 0.1) %>%
  group_by(speaker, id) %>%
  mutate(start=measurement_no==min(measurement_no)) %>%
  ungroup() %>%
  mutate(speaker=as.factor(paste("speaker", as.numeric(as.factor(speaker)))),
         traj=as.factor(id)) %>%
  group_by(speaker) %>%
  mutate(n=length(unique(id))) %>%
  ungroup() %>%
  rename(measurement.no="measurement_no") %>%
  dplyr::select(speaker, measurement.no, f1, f2, traj, start, n)

saveRDS(price_to_save, "~/Documents/Research/projects/dynamic-gam/gamm_strategies/data/final_data/price_vd_30_speakers.rds")

ggplot(price_to_save, aes(x=measurement.no, y=f2, group=traj)) +
  geom_line(alpha=0.5) +
  facet_wrap(~ speaker)

# generating some descriptive stats
# running by-speaker mixed effects models to fit lines to initial components of curves

start.speaker.mods <- list()
for (s in 1:length(unique(price_to_save$speaker))) {
  spe <- unique(price_to_save$speaker)[s]
  d <- filter(price_to_save, speaker==spe & measurement.no > 2 & measurement.no <= 5)
  initial.mod <- lmer(f2 ~ measurement.no + (1 + measurement.no | traj),
                       data=d,
                       control=lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb",
                       starttests = FALSE, kkt = FALSE))
                  )
  start.speaker.overall.mean <- mean(predict(initial.mod, 
                                        newdata=data.frame(measurement.no=3, f2=0, traj=d[1,"traj"]),
                                        re.form=NA))
  start.speaker.overall.sd <- sd(predict(initial.mod, 
                                         newdata=data.frame(measurement.no=3, f2=0, traj=unique(d$traj)),
                                         re.form=NULL))
  start.speaker.residual.sd <- sd(resid(initial.mod))
  start.speaker.mods[[s]] <- list(speaker.mean=start.speaker.overall.mean,
                            traj.sd=start.speaker.overall.sd,
                            resid.sd=start.speaker.residual.sd)
}

# population mean (across speakers, each speaker equally weighted)
mean(unlist(lapply(start.speaker.mods, function(x) x$speaker.mean)))
# population sd (across speakers, each speaker equally weighted)
sd(unlist(lapply(start.speaker.mods, function(x) x$speaker.mean)))
# standard deviation across trajs
mean(unlist(lapply(start.speaker.mods, function(x) x$traj.sd)))
# residual sd
mean(unlist(lapply(start.speaker.mods, function(x) x$resid.sd)))

#############################
# now the final components
#############################

end.speaker.mods <- list()
for (s in 1:length(unique(price_to_save$speaker))) {
  spe <- unique(price_to_save$speaker)[s]
  d <- filter(price_to_save, speaker==spe & measurement.no >= 8 & measurement.no <= 10)
  final.mod <- lmer(f2 ~ measurement.no + (1 + measurement.no | traj),
                      data=d,
                      control=lmerControl(optimizer = "optimx", calc.derivs = FALSE, optCtrl = list(method = "nlminb",
                                                                                                    starttests = FALSE, kkt = FALSE))
  )
  end.speaker.overall.mean <- mean(predict(final.mod, 
                                             newdata=data.frame(measurement.no=9, f2=0, traj=d[1,"traj"]),
                                             re.form=NA))
  end.speaker.overall.sd <- sd(predict(final.mod, 
                                         newdata=data.frame(measurement.no=9, f2=0, traj=unique(d$traj)),
                                         re.form=NULL))
  end.speaker.residual.sd <- sd(resid(final.mod))
  end.speaker.mods[[s]] <- list(speaker.mean=end.speaker.overall.mean,
                                  traj.sd=end.speaker.overall.sd,
                                  resid.sd=end.speaker.residual.sd)
}

# population mean (across speakers, each speaker equally weighted)
mean(unlist(lapply(end.speaker.mods, function(x) x$speaker.mean)))
# population sd (across speakers, each speaker equally weighted)
sd(unlist(lapply(end.speaker.mods, function(x) x$speaker.mean)))
# standard deviation across trajs
mean(unlist(lapply(end.speaker.mods, function(x) x$traj.sd)))
# residual sd
mean(unlist(lapply(end.speaker.mods, function(x) x$resid.sd)))