# RSD & Similarity Code Part 2/5

# 2. inherent_noise_encode_01_normalization.r
#     a. Initialize data
#     b. Jitter and Fang plots for raw data
#     c. Subset data
#     d. Create and initialize noise list
#     e. iscore correction
#     f. Density plots

# 2.
################################################################################
############################### Normalization ##################################
################################################################################

# 2.a Initialize data

load("/ENCODE-3_replicates_bigWig_hg38.RData")
load("/ENCODE-3_replicates_bigWig_metadata_hg38.RData")

idx <- colnames(DNase)[-(1:3)]
sample_metadata <- sample_metadata[match(idx,sample_metadata$file_accession),]
no_stem_idx <- sample_metadata[ sample_metadata$state %in% c("Primary", "Immortalized") , "file_accession" ]

scaffolds  <- list()
scaffolds[["continuous"]]  <- read.table(
  "D:/Axiotl/RSD_and_Similarity.R/encode-3_scaffold_sorted_hg38.bed",
  as.is = T)
colnames(scaffolds[["continuous"]]) <- c("chr","start","end","clust")

sample_metadata$id <- as.factor( c(
  rep("Fibroblasts", 14), rep("Stromal", 2), rep("Fibroblasts", 4),
  rep("Lymphoblastoid", 6), rep("ESC", 6), rep("Astrocytes", 6),
  rep("Epithelial", 2), rep("Endothelial",2), rep("Smooth Muscle", 2),
  rep("Fibroblasts", 4), rep("Cardiomyocytes", 2), rep("Fibroblasts", 2),
  rep("Epithelial", 4), rep("Fibroblasts", 6), rep("Epithelial", 4),
  rep("Fibroblasts", 2), rep("Endothelial", 14), rep("Epithelial", 2),
  rep("Fibroblasts", 6), rep("Epithelial", 4), rep("Endothelial", 2),
  rep("Epithelial", 2), rep("Myoblasts", 4), rep("Endothelial", 2),
  rep("Fibroblasts", 2), rep("Myoblasts", 4), rep("Astrocytes", 2),
  rep("Epithelial", 2), rep("Fibroblasts", 4), rep("Keratinocytes", 2),
  rep("Fibroblasts", 2), rep("Epithelial", 2), rep("Tubular", 2),
  rep("Epithelial", 2), rep("Fibroblasts", 4), rep("iPSC", 2),
  rep("ESC", 2), rep("SC", 6), rep("ESC", 2),
  rep("Tubular", 2), rep("Fibroblasts", 2), rep("Immortalized", 2),
  rep("ESC",2), rep("iPSC", 2), rep("SC", 2), rep("Adipocytes", 2),
  rep("Immortalized", 2), rep("SC", 2), rep("iPSC", 2),
  rep("ESC", 2), rep("Immortalized", 2), rep("ESC", 4) ) )


# 2.b Jitter and Fang plots for raw data

# fang plots for whole dataset
fang.raw <- TRUE

if (fang.raw){
  fang(data=DNase, idx=idx, path_out="D:/R/test/grid_fang/raw.pdf")
}

# fang plots for random subset w/ n=4, used to track changes through normalization
fang.tracker.raw <- TRUE

if(fang.tracker.raw){
  fang(data=DNase.track.raw, idx=sample.ids, path_out="D:/R/test/grid_fang/raw.track.pdf")
}

# 2.c Subset data: may wish to use subset for faster processing aka sanity check
# Subset is random and enriched for raw elements high activity

flag_subset <- TRUE

# changing activity threshold will change subset size. Lower threshold = more data.
# Threshold of 0.1 gives n = 663,140
# Threshold of 0.5 gives n = 155,952
# Threshold of 1 gives n = 85,241

activity_threshold <- .5

if (flag_subset){

  rand <- DNase[sample(nrow(DNase), 50000, replace = FALSE, prob = NULL),]
  rand.rows <- rownames(rand)
  activity <- as.data.frame(apply(DNase[,-c(1:4)], 1, mean))
  colnames(activity) <- "raw.activity"
  activity <- subset(activity, raw.activity > activity_threshold)
  act.rows <- rownames(activity)
  rand.rows <- rownames(rand)
  sub <- sort(unique(c(act.rows, rand.rows)))
  DNase <- DNase %>% filter(row_number() %in% sub)
  scaffolds[["continuous"]] <- scaffolds[["continuous"]] %>% filter(row_number() %in% sub)
  rm(activity, rand, rand.rows, act.rows, sub)
}

# 2.d Create and initialize noise list

noise <- list()

noise[["z_scores"     ]]   <- list()
noise[["i_scores"     ]]   <- list()
noise[["i_scores-cor"  ]]   <- list()
noise[["log2-i_scores"]]   <- list()
noise[["log2-z_scores"]]   <- list()
noise[["inh-units"    ]]   <- list()

noise[["z_scores"     ]][["whole_dataset"]]   <-     scale( DNase[,idx], center = T, scale = T)
noise[["i_scores"     ]][["whole_dataset"]]   <-  inh_norm( DNase, scaffolds[["continuous"]], iscores = T)
noise[["inh-units"    ]][["whole_dataset"]]   <-  inh_norm( DNase, scaffolds[["continuous"]], iscores = F)

# 2.e iscore correction: # adjust i-scores to make delta rsd comparable to z-scores
# (iscores-cor) by subtracting the mean of the mean and dividing by the mean of the sds

correct_iscores <- TRUE

if (correct_iscores){
  iscores <- noise[["i_scores"]][["whole_dataset"]][,-c(1:4)]
  iscores.means<-(rowMeans(iscores))
  meanofmeans <- mean(iscores.means)
  iscores.sds <- apply(iscores, 1, sd)
  meanofsds <- mean(iscores.sds)
  iscores <- (iscores - meanofmeans) / meanofsds
  head(iscores)
  rm(iscores.means, iscores.sds, meanofmeans, meanofsds)
  noise[[ "i_scores-cor" ]][["whole_dataset"]] <- cbind(noise[["i_scores"]][["whole_dataset"]][,c(1:4)],iscores)
  rm(iscores)
}

# 2.f Density plots
#Raw
raw.activity <- apply(DNase[,-c(1:4)], 1, mean)
raw.activity <- cbind(DNase[,c(1:4)], raw.activity)
raw.activity <- raw.activity %>% filter(clust=="pro_on" | clust=="enh_off")
raw.density <- density(raw.activity$raw.activity, bw = 0.07)
rm(raw.activity)

#i-scores
i.activity <- noise$i_scores$whole_dataset
ia <- apply(i.activity[,-c(1:4)], 1, mean)
i.activity <- cbind(i.activity[,c(1:4)], ia)
i.activity <- i.activity %>% filter(clust=="pro_on" | clust=="enh_off")
i.density <- density(i.activity$ia)

#inh-units
inh.activity <- noise$`inh-units`$whole_dataset
inh.a <- apply(inh.activity[,-c(1:4)], 1, mean)
inh.activity <- cbind(inh.activity[,c(1:4)], inh.a)
inh.activity <- inh.activity %>% filter(clust=="pro_on" | clust=="enh_off")
inhdens <- density(inh.activity$inh.a)

plot(raw.density, lwd = 2, xlim=c(-0.1,2), ylim=c(0,4.5),
       main="Pro_on & Enh_off Density", col="#2cccff")
lines(i.density, type ="l", col = "#eeb266", lwd = 2)
lines(inhdens, type = "l", lwd = 2, col="#e4778f")
legend(1.5, 4, legend=c("raw data", "i-scores", "inh-units"),
       col=c("#2cccff", "#eeb266", "#e4778f"), lty=1, cex=0.8)

sum(DNase$clust == "pro_on")
sum(DNase$clust == "enh_off")

rm(raw.activity, raw.density, i.activity, ia, i.density, inh.activity, inh.a, inhdens)
