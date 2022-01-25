# RSD & Similarity Code Part 3/5

# 3. inherent_noise_encode_02_sd_rsd.r
#     a. activity, rsd, sd calculations
#     b. banding
#     c. band plots

# 3.
###############################################################################
################################# sd, rsd #####################################
###############################################################################

# Issues: asympotote in delta rsd plot


# 3.a Activity, rsd, sd

# i_scores
noise[["i_scores"]][["activity"]]  <-  apply(noise[["i_scores"]][["whole_dataset"]][,idx], 1, mean)
noise[["i_scores"]][["rsd"     ]]  <-  apply(noise[["i_scores"]][["whole_dataset"]][,idx], 1, rsd)
noise[["i_scores"]][["sd"      ]]  <-  apply(noise[["i_scores"]][["whole_dataset"]][,idx], 1, sd)

# i_scores-cor
noise[["i_scores-cor"]][["activity"]]  <-  apply(noise[["i_scores-cor"]][["whole_dataset"]][,idx], 1, mean)
noise[["i_scores-cor"]][["rsd"     ]]  <-  apply(noise[["i_scores-cor"]][["whole_dataset"]][,idx], 1, rsd)
noise[["i_scores-cor"]][["sd"      ]]  <-  apply(noise[["i_scores-cor"]][["whole_dataset"]][,idx], 1, sd)

# z_scores
noise[["z_scores"]][["rsd"]]   <- apply(noise[["z_scores"]][["whole_dataset"]],1,rsd)
noise[["z_scores"]][["sd" ]]   <- apply(noise[["z_scores"]][["whole_dataset"]],1,sd)

# inh-units
noise[["inh-units"]][["activity"]]   <- apply(noise[["inh-units"]][["whole_dataset"]][,idx],1, mean)
noise[["inh-units"]][["rsd"     ]]   <- apply(noise[["inh-units"]][["whole_dataset"]][,idx],1, rsd)
noise[["inh-units"]][["sd"      ]]   <- apply(noise[["inh-units"]][["whole_dataset"]][,idx],1, sd)


# 3.b banding

score_names  <- c("z_scores", "i_scores", "inh-units")

for (score_name in score_names){
  if (score_name == "z_scores"){
    activity <- noise[[score_name]][["activity"]]
    activity <- as.data.frame(activity)
    activity <- activity %>% rowwise() %>% mutate(bands = band_code(activity))
    noise[[score_name]][["band"]] <- activity$bands
  }
  else if (score_name =="i_scores"){
    activity <- noise[[score_name]][["activity"]]
    activity <- as.data.frame(activity)
    activity <- activity %>% rowwise() %>% mutate(bands = band_code(activity))
    noise[[score_name]][["band"]] <- activity$bands
  }
  else if (score_name == "inh-units"){
    activity <- noise[[score_name]][["activity"]]
    activity <- as.data.frame(activity)
    activity <- activity %>% rowwise() %>% mutate(bands = band_code(activity))
    noise[[score_name]][["band"]] <- activity$bands
  }
}

sum(str_count(noise[["i_scores"]][["band"]], "A"))
sum(str_count(noise[["i_scores"]][["band"]], "B"))
sum(str_count(noise[["i_scores"]][["band"]], "C"))
sum(str_count(noise[["i_scores"]][["band"]], "D"))

# 3.c band plots

# delta sd
plot_data <- data.frame( activity = abs( noise[["i_scores"]][["activity"]]) ,
                         sd_diff = c( noise[["i_scores"]][["sd"]] - noise[["z_scores"]][["sd"]]) ,
                         band = noise[["i_scores"]][["band"]])
ggplot(plot_data) +
  geom_point( aes(activity, sd_diff, color=band), size=0.5, alpha=0.05)  +
  theme_classic() +
  xlab("Activity") +
  xlim(-0.1,2) +
  #scale_colour_manual(values = c("#758994","#0091c3","#a47c83","#f44746","#da6500")) +
  scale_colour_manual(values = c("#0091c3","#a47c83","#f44746","#da6500")) +
  theme(legend.position = "none") +
  ylab( expression(SD[i_scores] - SD[z_scores])) +
  ylim(-.75,0.1)  +
  geom_hline(yintercept= 0, color="red", lty=2, lwd=0.5)


# delta rsd (asymptote)
plot_data <- data.frame( activity=  abs( noise[["i_scores"]][["activity"]] ) ,
                         rsd_diff=    c( noise[["i_scores"]][["rsd"]] - noise[["z_scores"]][["rsd"]]),
                         band = noise[["i_scores"]][["band"]])

ggplot(plot_data) +
  geom_point( aes(activity, rsd_diff, color=band), size=0.05, alpha=0.05)  +
  theme_classic() +
  xlab("Activity") +
  xlim(0,1) +
  scale_colour_manual(values = c("#758994","#0091c3","#a47c83","#f44746","#da6500")) +
  theme(legend.position = "none") +
  ylab( expression(RSD[i_scores] - RSD[z_scores])) +
  ylim(-1500,1000)  +
  geom_hline(yintercept= 0, color="red", lty=2, lwd=0.5)
