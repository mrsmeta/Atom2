# RSD & Similarity Code Part 4/5

# 4. inherent_noise_encode_03_similarity.r
#     a. entropy, jaccard, and similarity calculations
#     b. entropy, jaccard, and similarity plots

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
