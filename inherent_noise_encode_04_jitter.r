# RSD & Similarity Code Part 5/5

# 5. inherent_noise_encode_04_jitter.r
#     a. jitter plots for normalized data
#     b. fang plots for normalized data

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
