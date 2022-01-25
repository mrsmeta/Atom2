# RSD & Similarity Code - Full


########################### Table of Contents #############################
# 1. inherent_noise_encode_00_functions.r
#     a. librarires
#     b. functions
#         1. geometric mean: gives geometric mean for various score types
#         2. fang: gives density distributions of pro_on, pro_off, enh_on, and
#                  enh_off elements.
#         3. inh-norm: gives inh-units and i-scores
#         4. rsd: gives residual standard devation
#         5. band codes: assigns bands A-D based on activity
#         6. entropy: caculates entropy
#         7. jaccard: gives us jaccard similarity index values
#         8. mus: gives mu values for z-scores, i-scores, and inh-units
# 2. inherent_noise_encode_01_normalization.r
#     a. Initialize data
#     b. Jitter and Fang plots for raw data
#     c. Subset data
#     d. Create and initialize noise list
#     e. iscore correction
#     f. Density plots
# 3. inherent_noise_encode_02_sd_rsd.r
#     a. activity, rsd, sd calculations
#     b. banding
#     c. band plots
# 4. inherent_noise_encode_03_similarity.r
#     a. entropy, jaccard, and similarity calculations
#     b. entropy, jaccard, and similarity plots
# 5. inherent_noise_encode_04_jitter.r
#     a. jitter plots for normalized data
#     b. fang plots for normalized data



################################################################################
############################ librarires & functions ############################
################################################################################
# 1.a librarires

library(mixtools)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gplots)
library(reshape2)
library(network)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(tibble)
library(tidyr)
options(scipen=333)

# 1.b functions
# 1.b.1 geometric mean: calculates geometric means
gm_mean <- function(x, na.rm=F){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# 1.b.2 fang: gives density distributions of pro_on, pro_off, enh_on, and
#             enh_off elements

fang <- function(data, idx, path_out){

  module_ids         <- c("on","off","bg")
  module_ids_pro_on  <- "pro_on"
  module_ids_pro_off <- "pro_off"
  module_ids_enh_on  <- "enh_on"
  module_ids_enh_off <- "enh_off"
  module_ids_bg      <- "bg"

  X <- list();                                                            C <- list()
  X[["pro_on" ]] <- data[data[,"clust"] %in% module_ids_pro_on,  -(1:4)]; C[["pro_on" ]] <- "tomato2"
  X[["pro_off"]] <- data[data[,"clust"] %in% module_ids_pro_off, -(1:4)]; C[["pro_off"]] <- "deepskyblue1"
  X[["enh_on" ]] <- data[data[,"clust"] %in% module_ids_enh_on,  -(1:4)]; C[["enh_on" ]] <- "gold"
  X[["enh_off"]] <- data[data[,"clust"] %in% module_ids_enh_off, -(1:4)]; C[["enh_off"]] <- "turquoise"
  X[["bg"     ]] <- data[data[,"clust"] %in% module_ids_bg,      -(1:4)]; C[["bg"     ]] <- "#888888"

  mus <- matrix(nrow = length(idx), ncol = 2, data = 0)
  colnames(mus) <- c("pro_on","enh_off")
  rownames(mus) <- idx

  for (i in 1:length(X)) {
    if (length(X[[i]]) == 0) {
      stop(sprintf("Input does not contain module: \"%s\"", names(X)[i]))
    }
  }
  mus              <- matrix( nrow = length( idx ), ncol = 4, data = 0 )
  colnames(mus)  <- c("pro_off","pro_on","enh_off","enh_on")
  rownames(mus)  <- idx

  for(i in 1:length(idx)){

    id <- idx[i]
    cat(sprintf( "%-7s \n", id ))

    error <- FALSE

    ilam <- c(   2/3,  1/3  )
    imus <- c(     0,    1  )
    isig <- c(     1,    1  )

    mix <- list()

    mix_or_no_mix <- function(data, lambda, mu, sigma, maxit){

      result <-  tryCatch({

        message("Trying to fit Gaussian Mixture Models to determine on and off means")
        normalmixEM( log2(data), lambda, mu, sigma, verb = FALSE, maxit)

      },

      error=function(cond) {
        message(paste("This sample returns error while fitting GMMs ", id, sep="- "))
        message("Proceeding without GMMs for this sample")

        error <<- TRUE
      })

      if(error){

        result <- list()
        result[["lambda"]]   <- c( log2(gm_mean(data))        , log2(gm_mean(data))        )
        result[["mu"    ]]   <- c( log2(gm_mean(data))        , log2(gm_mean(data))        )
        result[["sigma" ]]   <- c( sd(log2((data + 2^(-10)))) , sd(log2((data + 2^(-10)))) )
      }

      return(result)
    }

    mix[["pro_on"  ]]  <- mix_or_no_mix(data = X[["pro_on" ]][,i] + 2^(-10), lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)
    mix[["enh_off" ]]  <- mix_or_no_mix(data = X[["enh_off"]][,i] + 2^(-10), lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)

    mus[id,"pro_on" ]  <- max(mix[["pro_on" ]]$mu)
    mus[id,"enh_off"]  <- min(mix[["enh_off" ]]$mu)

    mus[id,"pro_off" ] <- log2(gm_mean(X[["pro_off"  ]][,id]))
    mus[id,"enh_on"  ] <- log2(gm_mean(X[["enh_on"   ]][,id]))
  }

  highest_var_unit <- names(sort(apply(mus[,1:2],2,var))[2])
  mus_ordered      <- mus[order(mus[,highest_var_unit], decreasing = T),]
  idx_main         <- rownames(mus_ordered)
  num_samples      <- length(idx)

  scale <-  6

  mean_on_raw  <- mean( apply( log2( X[["pro_on"  ]] + 0.01 ), 2, mean) )
  mean_off_raw <- mean( apply( log2( X[["enh_off" ]] + 0.01 ), 2, mean) )
  xlimL_raw    <- mean_off_raw - (mean_on_raw - mean_off_raw) - 1
  xlimR_raw    <- mean_on_raw  + (mean_on_raw - mean_off_raw) + 1

  size_each_fang <- 3
  size_of_pdf    <- ceiling(sqrt(num_samples))

  pdf(path_out, onefile = T,
      width  = size_of_pdf * size_each_fang,
      height = size_of_pdf * size_each_fang, useDingbats = F)

  par(mfrow=c(size_of_pdf,size_of_pdf), mar=c(2,2,3,2))

  for ( i in (1:length(idx_main))){

    id <- idx_main[i]
    cat(sprintf( "%-7s ", id ))

    plot( NA, ylim = c(0,1), xlim = c(xlimL_raw,xlimR_raw),
          ylab = NA, xlab = NA,
          axes = FALSE )

    mtext(side=3, line= 0.9, adj=0, cex=  1,    id, col="gray28")

    box(lty="solid", col="gray50")

    abline( h =                     2/3, col = "lightgray")
    # abline( h =                  1.35/3, col =    "gray90")
    # abline( h =                  1.65/3, col =    "gray90")
    abline( h =                     1/3, col = "lightgray")
    abline( h =                       0, col = "lightgray")
    abline( v =   mean(mus[,"enh_off"]), col =    "gray28")
    abline( v =   mean(mus[, "pro_on"]), col =    "gray28")

    for ( test in c( "bg", "pro_off", "pro_on" )) {
      dens <- density( log2( X[[test]][,id] ) ) #[,i]
      points( dens$x, (dens$y/scale) + 2/3, type="l", col = C[[test]], lwd = 2 )
    }

    for ( test in c( "enh_off","enh_on" )) {
      dens <- density( log2( X[[test]][,id] ) ) #[,i]
      points( dens$x, (dens$y/scale)      , type="l", col = C[[test]], lwd = 2 )
    }

    for ( test in c( "pro_off",  "pro_on") ){

      m <- mus_ordered[id,test]
      lines( c( m, m ), c(0,0.0015) + 1.5/3, col = C[[test]], lwd = 9 )
      text( m, 0.035 + 1.65/3, sprintf("%4.2f", m ), col = C[[test]], cex = 1.5)
    }

    for ( test in c( "enh_off", "enh_on" ) ){

      m <- mus_ordered[id,test]
      lines( c( m, m ), c(0,0.0015) + 1.3/3, col = C[[test]], lwd = 9 )
      text( m, 0.035 + 1/3, sprintf("%4.2f", m ), col = C[[test]], cex = 1.5)
    }
  }
  dev.off()
}


# 1.b.3 inh-norm: gives inh-units and i-scores

# issues:
# this function uses at least one variable, iscores,
# that isn't passed in the arguments

inh_norm <- function(raw_data, scaffold, iscores){

  clust <- scaffold[,"clust"]
  raw_data <- cbind( DNase[ ,c("chr","start","end") ], clust, DNase[ ,idx ] ) ; rm(clust)
  module_ids_pro_on  <- "pro_on"; module_ids_enh_off <- "enh_off"

  X <- list()
  X[["pro_on"]] <- raw_data[ raw_data[,"clust"] %in% module_ids_pro_on, -(1:4)]
  X[["enh_off"]] <- raw_data[ raw_data[,"clust"] %in% module_ids_enh_off, -(1:4)]

  C <- list()
  C[["pro_on"]] <- "tomato2"
  C[["enh_off"]] <- "turquoise"

  mus <- matrix(nrow = length(idx), ncol = 2, data = 0)
  colnames(mus) <- c("pro_on","enh_off")
  rownames(mus) <- idx

  no_mix_mean <- function(x){
    2^( mean( log2( x + 2^(-10))))
  }

  for ( i in 1:length(idx)){
    id <- idx[i]
    cat(sprintf( "%-7s \n", id ))

    error <- FALSE

    ilam <- c(   2/3,  1/3  )
    imus <- c(     0,    1  )
    isig <- c(     1,    1  )

    mix <- list()

    mix_or_no_mix <- function(data, lambda, mu, sigma, maxit){

      result <-  tryCatch({

        if(iscores){

          message("Trying to fit Gaussian Mixture Models to determine on and off means")
          normalmixEM( data, lambda, mu, sigma, verb = FALSE, maxit)
        }
        else{
          message("Trying to fit Gaussian Mixture Models to determine on and off means")
          normalmixEM( log2(data), lambda, mu, sigma, verb = FALSE, maxit)
        }
      },

      error=function(cond) {
        message(paste("This sample returns error while fitting GMMs ", id, sep="- "))
        message("Proceeding without GMMs for this sample")

        error <<- TRUE
      })

      if(error){
        if(iscores){
          result <- list()
          result[["lambda"]] <- c( (no_mix_mean(data))        , (no_mix_mean(data))        )
          result[["mu"    ]] <- c( (no_mix_mean(data))        , (no_mix_mean(data))        )
          result[["sigma" ]] <- c( sd(((  data     )))        ,  sd(((  data    )))        )
        }
        else{
          result <- list()
          result[["lambda"]] <- c( log2(no_mix_mean(data))        , log2(no_mix_mean(data))        )
          result[["mu"    ]] <- c( log2(no_mix_mean(data))        , log2(no_mix_mean(data))        )
          result[["sigma" ]] <- c( sd(log2((data + 2^(-10))))     , sd(log2((data + 2^(-10)))) )
        }
      }
      return(result)
    }

    if(iscores){

      mix[["pro_on"  ]]  <- mix_or_no_mix(data = X[["pro_on" ]][,i] , lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)
      mix[["enh_off" ]]  <- mix_or_no_mix(data = X[["enh_off"]][,i] , lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)

      mus[id,"enh_off"]  <- min(mix[["enh_off"]]$mu)
      mus[id,"pro_on" ]  <- max(mix[["pro_on" ]]$mu)
    }
    else{
      mix[["pro_on"  ]]  <- mix_or_no_mix(data = X[["pro_on" ]][,i] + 2^(-10), lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)
      mix[["enh_off" ]]  <- mix_or_no_mix(data = X[["enh_off"]][,i] + 2^(-10), lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)

      mus[id,"enh_off"]  <- min(mix[["enh_off"]]$mu)
      mus[id,"pro_on" ]  <- max(mix[["pro_on" ]]$mu)
    }
  }

  if(iscores){

    norm_data <- raw_data
    for (id in idx ){
      cat(sprintf("%s ", id))
      mu0 <- mus[id,"enh_off"]
      mu1 <- mus[id,"pro_on" ]
      norm_data[,id] <-  (raw_data[,id]  - mu0 )/( mu1 - mu0)
    }
  }
  else{
    norm_data <- raw_data
    for (id in idx ){
      cat(sprintf("%s ", id))
      mu0 <- mus[id,"enh_off"]
      mu1 <- mus[id,"pro_on" ]
      norm_data[,id] <- 2^((log2(raw_data[,id] + 2^(-10)) - mu0 )/( mu1 - mu0 ))
    }
  }
  rm(module_ids_enh_off, module_ids_pro_on,
     raw_data, X, no_mix_mean, id, ilam, imus, isig, error, mix,
     mix_or_no_mix, mus, mu0, mu1)
  return(norm_data)
}


# 1.b.4 rsd: gives us residiual standard deviation
rsd <- function(vec){
  mean_vec <- mean(vec)
  rsd_vec <- (sd(vec)/abs(mean_vec))*100
  return(rsd_vec)
}

# 1.b.5 band code: assigns bands A-D based on i-score activity, band_code.inh
# assigns bands A-D based on inh-activity in log2 space (?)
band_code <- function(x){
  x <- x[1]
  if(x < (-0.25)){color <- "A"}
  else if(x > (-0.25) && x < (0.25)){color <- "B"}
  else if(x > (0.25) && x < (0.75)){color <- "C"}
  else if(x > (0.75)){color <- "D"}
  return(color)
}

band_code.inh <- function(x){
  x <- x[1]
  if(x < 2^( 0 - 1/4 )){color <- "A"}
  else if(x > 2^( 0 - 1/4 ) && x < 2^( 0 + 1/4 )){color <- "B"}
  else if(x > 2^( 0 + 1/4 ) && x < 2^( 1 - 1/4 )){color <- "C"}
  else if(x > 2^( 1 - 1/4 )){color <- "D"}
  return(color)
}

# 1.b.6 entropy: calculates entropy
entropy <- function(vec){
  p <- sum(vec)/length(vec)
  return( -((p)*log2(p) + (1-p)*log2(1-p)))
}

# 1.b.7. jaccard: gives us jaccard similarity index values
jaccard <- function(x,y){
  # https://www.learndatasci.com/glossary/jaccard-similarity/
  a <- length( which( (x==1) & (y==1)))
  b <- length( which( (x==0) & (y==1)))
  c <- length( which( (x==1) & (y==0)))
  sim <- a/(a+b+c)
  return(sim)
}

# 1.b.8. mus: gives mu values for z-scores, i-scores, and inh-units

mus <- function(raw_data, scaffold, iscores, z_pivots){

  if(z_pivots){
    pivots <- data.frame( mean = apply( raw_data[,idx], 2, mean),
                          sd = apply( raw_data[,idx], 2, sd))
    rownames(pivots) <- idx
    return(pivots)
    break()
  }
  clust <- scaffolds[["continuous"]][["clust"]]
  raw_data <- DNase

  module_ids_pro_on  <- "pro_on"; module_ids_enh_off <- "enh_off"

  X <- list()
  X[["pro_on"]] <- raw_data[ raw_data[,"clust"] %in% module_ids_pro_on, -(1:4)]
  X[["enh_off"]] <- raw_data[ raw_data[,"clust"] %in% module_ids_enh_off, -(1:4)]

  C <- list()
  C[["pro_on"]] <- "tomato2"
  C[["enh_off"]] <- "turquoise"

  mus <- matrix(nrow = length(idx), ncol = 2, data = 0)
  colnames(mus) <- c("pro_on","enh_off")
  rownames(mus) <- idx

  no_mix_mean <- function(x){
    2^( mean(log2( x + 2^(-10))))
  }

  for ( i in 1:length(idx)){
    id <- idx[i]
    cat(sprintf( "%-7s \n", id ))
    error <- FALSE

    ilam <- c(   2/3,  1/3  )
    imus <- c(     0,    1  )
    isig <- c(     1,    1  )

    mix <- list()

    mix_or_no_mix <- function(data, lambda, mu, sigma, maxit){

      result <-  tryCatch({

        if(iscores){

          message("Trying to fit Gaussian Mixture Models to determine on and off means")
          normalmixEM( data, lambda, mu, sigma, verb = FALSE, maxit)
        }
        else{
          message("Trying to fit Gaussian Mixture Models to determine on and off means")
          normalmixEM( log2(data), lambda, mu, sigma, verb = FALSE, maxit)
        }
      },
      error=function(cond) {
        message(paste("This sample returns error while fitting GMMs ", id, sep="- "))
        message("Proceeding without GMMs for this sample")

        error <<- TRUE
      })

      if(error){

        if(iscores){

          result <- list()
          result[["lambda"]] <- c((no_mix_mean(data)), (no_mix_mean(data)))
          result[["mu"    ]] <- c((no_mix_mean(data)), (no_mix_mean(data)))
          result[["sigma" ]] <- c(sd(((data))),  sd(((data))))
        }
        else{

          result <- list()
          result[["lambda"]] <- c( log2(no_mix_mean(data))        , log2(no_mix_mean(data)))
          result[["mu"    ]] <- c( log2(no_mix_mean(data))        , log2(no_mix_mean(data)))
          result[["sigma" ]] <- c( sd(log2((data + 2^(-10))))     , sd(log2((data + 2^(-10)))))
        }
      }
      return(result)
    }
    if(iscores){

      mix[["pro_on"  ]]  <- mix_or_no_mix(data = X[["pro_on" ]][,i] , lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)
      mix[["enh_off" ]]  <- mix_or_no_mix(data = X[["enh_off"]][,i] , lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)

      mus[id,"enh_off"]  <- min(mix[["enh_off"]]$mu)
      mus[id,"pro_on" ]  <- max(mix[["pro_on" ]]$mu)
    }
    else{

      mix[["pro_on"  ]]  <- mix_or_no_mix(data = X[["pro_on" ]][,i] + 2^(-10), lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)
      mix[["enh_off" ]]  <- mix_or_no_mix(data = X[["enh_off"]][,i] + 2^(-10), lambda = ilam[1], mu = imus, sigma = isig^2, maxit = maxit)

      mus[id,"enh_off"]  <- min(mix[["enh_off"]]$mu)
      mus[id,"pro_on" ]  <- max(mix[["pro_on" ]]$mu)
    }
  }

  if(iscores){

    norm_data <- raw_data
    for (id in idx ){
      cat(sprintf("%s ", id))
      mu0 <- mus[id,"enh_off"]
      mu1 <- mus[id,"pro_on" ]
      norm_data[,id] <- (raw_data[,id]  - mu0)/(mu1 - mu0)
    }
  }

  else{
    norm_data <- raw_data
    for (id in idx ){
      cat(sprintf("%s ", id))
      mu0 <- mus[id,"enh_off"]
      mu1 <- mus[id,"pro_on" ]
      norm_data[,id] <- 2^( ( log2(raw_data[,id] + 2^(-10) ) - mu0 )/( mu1 - mu0 ) )
    }
  }
  # rm(module_ids_enh_off, module_ids_pro_on,
  #    raw_data, X, no_mix_mean, id, ilam, imus, isig, error, mix,
  #    mix_or_no_mix, norm_data, mu0, mu1)
  return(mus)
}

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


# 4.
################################################################################
############################## similarity ######################################
################################################################################

## issues:
# 1. Is stand-alone bands list redundant? Can we merge it with noise list?
# WIP, notes below.

# Attempted to add bands list with noise list to simplify, but it seems to make things more complicated.
# All entropy/jaccard/similarity loops in next section are based on bands list and data returns to bands list.
# This could be changed, but will require re-constructing the next section as well, which is currently working.

# 2. Banding could be determined by other score types. Is this significant?

# In this section, element IDs are sorted into bands based on activity levels of
# inh-units, whereas in previous section (rsd and sd), bands were sorted using i_score
# activity data.

# 3. Thresholds could be calculated by other score types. Is this significant?

# Thresholds are calculated with z_scores and are used in all similarity calculations.
# Thresholds could be calculated for inh-units and iscores, but will require
# more time and loops and re-structuring of similarity calculations. In the noise
# list, z-score whole dataset is double whereas iscore and inh-norm are lists so
# code must be indexed differenty

# Basic outline of attempt at merging band and noise lists and using alternate
# score type to calculate thresholds:

# noise[["inh-units"]][["band_ids"]][["A"]][["element_ids"]] <-
#                which(noise[["inh-units"]][["activity"]] < 2^( 0 - 1/4 ))
# noise[["inh-units"]][["band_ids"]][["B"]][["element_ids"]] <-
#                which(noise[["inh-units"]][["activity"]] > 2^( 0 - 1/4 )
#                & noise[["inh-units"]][["activity"]] < 2^( 0 + 1/4 ))
# noise[["inh-units"]][["band_ids"]][["C"]][["element_ids"]] <-
#                which(noise[["inh-units"]][["activity"]] > 2^( 0 + 1/4 )
#                & noise[["inh-units"]][["activity"]] < 2^( 1 - 1/4 ))
# noise[["inh-units"]][["band_ids"]][["D"]][["element_ids"]] <-
#                which(noise[["inh-units"]][["activity"]] > 2^( 1 - 1/4 ))

# Q <- quantile(noise[["inh-units"]][["whole_dataset"]][noise[["band_ids"]][["inh-units"]][["D"]][["element_ids"]], idx], seq(0,1,0.01))
# noise[["inh-units"]][["thresholds"]] <- round(seq(Q[2], Q[100], (Q[100] - Q[2])/199), 3)
#
# Q works if z_scores are used, gives error as is with inh-units



bands        <- list()
bands[["A"]] <- list()
bands[["B"]] <- list()
bands[["C"]] <- list()
bands[["D"]] <- list()

bands[["A"]][["element_ids"]] <-  which(noise[["inh-units"]][["activity"]] < 2^( 0 - 1/4 ))
bands[["B"]][["element_ids"]] <-  which(noise[["inh-units"]][["activity"]] > 2^( 0 - 1/4 )
                                  & noise[["inh-units"]][["activity"]] < 2^( 0 + 1/4 ) )
bands[["C"]][["element_ids"]] <-  which(noise[["inh-units"]][["activity"]] > 2^( 0 + 1/4 )
                                  & noise[["inh-units"]][["activity"]] < 2^( 1 - 1/4 ) )
bands[["D"]][["element_ids"]] <-  which(noise[["inh-units"]][["activity"]] > 2^( 1 - 1/4 ))


Q <- quantile(noise[["z_scores"]][["whole_dataset"]][bands[["D"]][["element_ids"]], idx], seq(0,1,0.01))
bands[["D"]][["thresholds"]] <- round(seq(Q.og[2], Q.og[100], (Q.og[100] - Q.og[2])/199), 3)


replicate_idx <- seq( 1, length(idx), 2)


# 4.a Jaccard, entropy, and similarity scores for i-scores, z-scores, and inh-units

score_names  <- c( "z_scores", "i_scores", "inh-units"  )

for (score_name in score_names){
  for(band in c("D")){

    cat(sprintf("%s ", band))

    binary_thresholds  <- bands[[band]][["thresholds"]]
    #binary_thresholds  <- noise[["z_scores"]][["thresholds"]]

    matrix_jaccard  <- matrix( nrow = length( replicate_idx ), ncol = length( binary_thresholds ), data = 0 )
    colnames(matrix_jaccard) <- as.character(binary_thresholds)
    rownames(matrix_jaccard) <- as.character(replicate_idx)

    matrix_entropy  <- matrix( nrow = length( replicate_idx ), ncol = length( binary_thresholds ), data = 0 )
    colnames(matrix_entropy) <- as.character(binary_thresholds)
    rownames(matrix_entropy) <- as.character(replicate_idx)

    matrix_similarity  <- matrix( nrow = length( replicate_idx ), ncol = length( binary_thresholds ), data = 0 )
    colnames(matrix_similarity) <- as.character(binary_thresholds)
    rownames(matrix_similarity) <- as.character(replicate_idx)

    for( threshold in binary_thresholds ){

      cat(sprintf("%s ", threshold))

      for( i in replicate_idx ){

        rep_1 <- noise[[score_name]][["whole_dataset"]][bands[[band]][["element_ids"]] , idx[i]   ]
        rep_2 <- noise[[score_name]][["whole_dataset"]][bands[[band]][["element_ids"]] , idx[i+1] ]

        rep_1 <- as.numeric( rep_1 > threshold )
        rep_2 <- as.numeric( rep_2 > threshold )

        rep_1 <- c( 1, rep_1, 0 )
        rep_2 <- c( 1, rep_2, 0 )

        jaccard_similarity <- jaccard( rep_1, rep_2 )

        rep_1_entropy <- entropy(rep_1)
        rep_2_entropy <- entropy(rep_2)

        matrix_jaccard[ as.character(i), as.character(threshold)] <- jaccard_similarity
        matrix_entropy[ as.character(i), as.character(threshold)] <- mean( rep_1_entropy , rep_2_entropy)
        matrix_similarity[ as.character(i), as.character(threshold)] <- jaccard_similarity * mean( rep_1_entropy , rep_2_entropy)

      }

      rm(rep_1, rep_2, rep_1_entropy, rep_2_entropy, jaccard_similarity)
    }
    bandsText = paste0(score_name, "_matrix_jaccard")
    bands[[band]][[bandsText]] <- matrix_jaccard
    bandsText = paste0(score_name, "_matrix_entropy")
    bands[[band]][[bandsText]] <- matrix_entropy
    bandsText = paste0(score_name, "_matrix_similarity")
    bands[[band]][[bandsText]] <- matrix_similarity
    rm( matrix_jaccard, matrix_entropy, matrix_similarity )
  }
}
rm(score_name, score_names)


# 4.b Similarity plots

score_names  <- c( "z_scores", "i_scores", "inh-units"  )
metric_names <- c( "jaccard", "entropy", "similarity" )
matrix_types <- c("_matrix_jaccard", "_matrix_entropy", "_matrix_similarity")

metric_colors <- list()
metric_colors[[ "jaccard"    ]] <- "#d35d7a"
metric_colors[[ "entropy"    ]] <- "#6483c5"
metric_colors[[ "similarity" ]] <- "#bd62b6"

path_out <- "D:/R/test/"


# 9 line plots with color

for (metric_name in metric_names){
  if (metric_name == "jaccard"){
    colVar<<- "#d35d7a"
    colType <<- "jaccard"
    matrixShifter <<- "_matrix_jaccard"
  }else if (metric_name == "entropy"){
    colVar<<- "#6483c5"
    colType <<- "entropy"
    matrixShifter <<- "_matrix_entropy"
  }else if (metric_name == "similarity"){
    colVar<<- "#bd62b6"
    colType <<- "similarity"
    matrixShifter <<- "_matrix_similarity"
  }
  for (score_name in score_names ){
    pdf(paste0(path_out, score_name, "_", metric_name, "_line.pdf"))
    print(ggplot(
      melt(bands[["D"]][[ paste0(score_name, matrixShifter) ]]) ) +
        geom_line( aes( col=colVar, x=as.factor(Var2) , y=value , group=as.factor(Var1)), alpha=0.5 ) +
        scale_color_manual(values = metric_colors[[ colType ]]) +
        theme_classic() +
        theme(legend.position = "none") +
        scale_x_discrete(breaks = ((bands[["D"]][["thresholds"]][c(T,F,F,F)]))) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.ticks.x = element_blank()) +
        #add.expr = abline(v=12.160) +
        xlab("Thresholds") +
        ylab( metric_name ) +
        ggtitle(paste(score_name, metric_name)) +
        theme(plot.title = element_text(hjust = 0.5)))
    dev.off()

    print(metric_name)
    print(score_name)
  }
}


# heatmaps

# z-score, i-score, inh-units entropy heatmaps
for (metric_name in metric_names){
  if (metric_name == "jaccard"){
    colVar<<- "#d35d7a"
    colType <<- "jaccard"
    matrixShifter <<- "_matrix_jaccard"
  }else if (metric_name == "entropy"){
    colVar<<- "#6483c5"
    colType <<- "entropy"
    matrixShifter <<- "_matrix_entropy"
  }else if (metric_name == "similarity"){
    colVar<<- "#bd62b6"
    colType <<- "similarity"
    matrixShifter <<- "_matrix_similarity"
  }
  for ( score_name in score_names ){
    mainText = paste("Band-D ", score_name, metric_name)
    pdf(paste0( path_out, score_name, " ", metric_name, " heatmap.pdf"))
    heatmap.2(bands[["D"]][[paste0(score_name, matrixShifter) ]],
      scale="none", dendrogram ="none",
      col = colorRampPalette(c("#ffffff",colVar))(100),
      trace="none", Rowv=NA, Colv=NA,
      main = mainText,
      xlab = "Thresholds", ylab = "Pairs")
dev.off()
  }
}

###############################################################################
################################## jitter #####################################
###############################################################################

# issues:

# 1. Plots not working correctly, needs attention
# 2. Need to get working for raw data and normalized data for tracking changes

noise[["z_scores"]][["pivots"]] <- mus(DNase, z_pivots = T)
noise[["i_scores"]][["pivots"]] <- mus(DNase, scaffolds[["continuous"]], iscores = T, z_pivots = F)
noise[["inh-units"]][["pivots"]] <- mus(DNase, scaffolds[["continuous"]], iscores = T, z_pivots = F)

num_samples <- length(idx)

xlimL <-  - 2
xlimR <-    5

# 1. i-scores

pdf("D:/R/test/jitter.pdf")
plot( NA, ylim = c(0,308), xlim = c(xlimL,xlimR), ylab = NA, xlab = NA, axes = FALSE ) #380
abline( h = seq(0,272,4), col = "lightgray", lwd = 0.5)                                #344
abline( h = 298, col = "lightgray", lwd = 0.5)                                         #370
abline( v = -1, col = "lightgray", lwd = 0.75)
text( -0.95, 0.05 + 303.5, sprintf("%10.1f", -1 ), col = "black", cex = 0.8)           #375.5
abline( v =  4, col = "lightgray", lwd = 1)
text(  3.95, 0.05 + 303.5, sprintf("%10.1f",  4 ), col = "black", cex = 0.8)
dev.off()

plot.new()
plotting_order <- sample_metadata[order(sample_metadata$id),"file_accession"]

for ( i in 1:length(plotting_order)){

  id <- plotting_order[i]

  x <- sort(c(rev(seq(0,272,4)),rev(seq(0,272,4))), decreasing = T)[i]

  C <- list()
  C[["pro_on"]] <- "tomato2" ; C[["enh_off"]] <- "turquoise"

  cat(sprintf( "%-7s ", id ))

  if(i %% 2 == 0){C[["pro_on"]] <- "tomato4" ; C[["enh_off"]] <- "turquoise4"}

  for ( test in c( "enh_off",  "pro_on") ){
    m <- noise[["i_scores"]][["pivots"]][id,test]
    lines( c( m, m ), c(0,0.5) + x, col = C[[test]], lwd = 2.5 )
  }


}

pdf("D:/R/test/jitter.pdf")
plot( NA, ylim = c(0,308), xlim = c(xlimL,xlimR), ylab = NA, xlab = NA, axes = FALSE ) #380
abline( h = seq(0,272,4), col = "lightgray", lwd = 0.5)                                #344
abline( h = 298, col = "lightgray", lwd = 0.5)                                         #370
abline( v = -1, col = "lightgray", lwd = 0.75)
text( -0.95, 0.05 + 303.5, sprintf("%10.1f", -1 ), col = "black", cex = 0.8)           #375.5
abline( v =  4, col = "lightgray", lwd = 1)
text(  3.95, 0.05 + 303.5, sprintf("%10.1f",  4 ), col = "black", cex = 0.8)
dev.off()



# 1. z-scores


plot( NA, ylim = c(0,308), xlim = c(xlimL,xlimR), ylab = NA, xlab = NA, axes = FALSE ) #380
abline( h = seq(0,272,4), col = "lightgray", lwd = 0.5)                                #344
abline( h = 298, col = "lightgray", lwd = 0.5)                                         #370
abline( v = -1, col = "lightgray", lwd = 0.75)
text( -0.95, 0.05 + 303.5, sprintf("%10.1f", -1 ), col = "black", cex = 0.8)           #375.5
abline( v =  4, col = "lightgray", lwd = 1)
text(  3.95, 0.05 + 303.5, sprintf("%10.1f",  4 ), col = "black", cex = 0.8)
abline( v =  0, col = "lightgray", lwd = 0.75)
abline( v =  1, col = "lightgray", lwd = 1)



plotting_order <- sample_metadata[order(sample_metadata$id),"file_accession"]

for ( i in 1:length(plotting_order)){

  id <- plotting_order[i]

  x <- sort(c(rev(seq(0,272,4)),rev(seq(0,272,4))), decreasing = T)[i]

  C <- list()
  C[["sd"]] <- "tomato2" ; C[["mean"]] <- "turquoise"

  cat(sprintf( "%-7s ", id ))

  if(i %% 2 == 0){C[["sd"]] <- "tomato4" ; C[["mean"]] <- "turquoise4"}

  for ( test in c( "mean",  "sd") ){
    m <- noise[["z_scores"]][["pivots"]][id,test]
    lines( c( m, m ), c(0,0.5) + x, col = C[[test]], lwd = 2.5 )
  }


}
