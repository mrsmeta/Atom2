.# RSD & Similarity Code Part 1/5

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
