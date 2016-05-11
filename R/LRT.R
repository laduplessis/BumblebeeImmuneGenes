# Functions for performing LRTs and for extracting data

#############################################################################################
# Do LRT, extract significant, non-significant

# Do LRT with specified number of degrees of freedom
# Datasets need to be in the same order!
doLRT <- function(alt, null, df, alpha) {
  
  if (sum(null$Alignment == alt$Alignment) == length(null$Alignment)) {
    
    lkratio <- 2*(alt$Lk - null$Lk)
    
    # To calculate mixtures of chi squared distributions
    if (length(df) > 1) {
       ptemp <- matrix(nrow=length(null$Alignment), ncol=length(df))
       for (i in 1:length(df)) {
         ptemp[,i] <- 1 - pchisq(2*(alt$Lk - null$Lk), df = df[i])
       }
       pvals <- rowSums(ptemp)/length(df)
    } else    
      pvals    <- 1 - pchisq(2*(alt$Lk - null$Lk), df = df)
       
    sign     <- pvals < alpha
    mtrpvals <- p.adjust(pvals, method="BH")
    mtrsign  <- mtrpvals < alpha
    
    res <- cbind(pvals, sign, mtrpvals, mtrsign, lkratio)
    colnames(res) <- c("P-val", "Significant", "BH P-val", "BH significant", "Lk-ratio")
    rownames(res) <- null$Alignment
    
    cat(sprintf("\t\t\t%d/%d significant (no correction)\n", sum(res[,2]), length(null$Alignment)))
    cat(sprintf("\t\t\t%d/%d significant after BH-correction\n\n", sum(res[,4]), length(null$Alignment)))
    
    return(res)
  } else
    stop("Datasets missing")
}

getSignificant <- function (test) {
  return(test[,4] == 1)
}

getNotSignificant <- function(test) {
  return(test[,4] == 0)
}

getSignificantNames <- function (test) {
  return(rownames(test)[test[,4] == 1])
}

getNotSignificantNames <- function(test) {
  return(rownames(test)[test[,4] == 0])
}


#############################################################################################
# Summarize the data in a model

# Summarize Branch Site Alternative such that the only fields are:
# - Alignment
# - Lk
# - Length
# - w2
# - p2
# - predicted number of sites
summarizeBSAlt <- function(BSAlt, Alignments) {
  BSsummary <- matrix(nrow=nrow(BSAlt), ncol=5)
  for (i in 1:nrow(BSAlt)) {
      prop <- BSAlt$p2a[i]+BSAlt$p2b[i]
      len  <- Alignments[which(Alignments[,1] == BSAlt[i,1]),2]
      BSsummary[i,] <- c(BSAlt$Lk[i], len, BSAlt$w2[i], prop, prop*len)
  }
  rownames(BSsummary) <- BSAlt$Alignment
  colnames(BSsummary) <- c('Lk', 'Length', 'w2', 'p2', 'nrsites')
  return(BSsummary)
}


# Summarize 2 Clade model such that the only fields are:
# - Alignment
# - Lk
# - Length
# - w2
# - w3
# - p2
# - predicted number of sites in p2
summarize2Clade <- function(CladeC, Alignments) {
  CladeCsummary <- matrix(nrow=nrow(CladeC), ncol=6)
  for (i in 1:nrow(CladeC)) {
      prop <- CladeC$p2[i]
      len  <- Alignments[which(Alignments[,1] == CladeC[i,1]),2]
      CladeCsummary[i,] <- c(CladeC$Lk[i], len, CladeC$w2[i], CladeC$w3[i], prop, prop*len)
  }
  rownames(CladeCsummary) <- CladeC$Alignment
  colnames(CladeCsummary) <- c('Lk', 'Length', 'w2', 'w3', 'p2', 'nrsites')
  return(CladeCsummary)
}

# Summarize M8 model such that the only fields are:
# - Alignment
# - Lk
# - Length
# - w10
# - p10
# - predicted number of sites in p10
summarizeM8 <- function(M8, Alignments) {
  M8summary <- matrix(nrow=nrow(M8), ncol=5)
  for (i in 1:nrow(M8)) {
      prop <- M8$p10[i]
      len  <- Alignments[which(Alignments[,1] == M8[i,1]),2]
      M8summary[i,] <- c(M8$Lk[i], len, M8$w10[i], prop, prop*len)
  }
  rownames(M8summary) <- M8$Alignment
  colnames(M8summary) <- c('Lk', 'Length', 'w10', 'p10', 'nrsites')
  return(M8summary)
}

# Summarize M3 model such that the only fields are:
# - Alignment
# - Lk
# - Length
# - w10
# - p10
# - predicted number of sites in p10
summarizeM3 <- function(M3, Alignments) {
  M3summary <- matrix(nrow=nrow(M3), ncol=5)
  for (i in 1:nrow(M3)) {
      prop <- M3$p2[i]
      len  <- Alignments[which(Alignments[,1] == M3[i,1]),2]
      M3summary[i,] <- c(M3$Lk[i], len, M3$w2[i], prop, prop*len)
  }
  rownames(M3summary) <- M3$Alignment
  colnames(M3summary) <- c('Lk', 'Length', 'w2', 'p2', 'nrsites')
  return(M3summary)
}

#############################################################################################
# Summarize the data in a model and extract significant and p-vals from a test


# Like summarizeBSAlt, but with p-values and automatically extracting significant
extractBSAlt <- function(BSAlt, Alignments, Test) {
  BSsummary <- matrix(nrow=sum(Test[,"BH significant"]), ncol=7)
  rows      <- matrix(nrow=nrow(BSsummary))
  j <- 1
  for (i in 1:nrow(BSAlt)) {
      prop <- BSAlt$p2a[i]+BSAlt$p2b[i]
      testidx <- which(rownames(Test) == BSAlt[i,1])
      len  <- Alignments[which(Alignments[,1] == BSAlt[i,1]),2]
      if (Test[testidx,"BH significant"] == 1) {
          BSsummary[j,] <- c(Test[testidx,"Lk-ratio"], len, BSAlt$w2[i], prop, prop*len,
                             Test[testidx,"P-val"], Test[testidx,"BH P-val"])
          rows[j]       <- rownames(Test)[testidx]
          j <- j + 1
      }
  }
  rownames(BSsummary) <- rows
  colnames(BSsummary) <- c('Lk-ratio', 'Length', 'w2', 'p2', 'nrsites', 'P-val', 'BH P-val')
  return(BSsummary)
}


# Like summarize 2 Clade, but with p-values and automatically extracting significant
extract2Clade <- function(Clade, Alignments, Test) {
  Cladesummary <- matrix(nrow=sum(Test[,"BH significant"]), ncol=8)
  rows         <- matrix(nrow=nrow(Cladesummary))
  j <- 1
  for (i in 1:nrow(Clade)) {
      prop <- Clade$p2[i]
      testidx <- which(rownames(Test) == Clade[i,1])
      len  <- Alignments[which(Alignments[,1] == Clade[i,1]),2]
      if (Test[testidx,"BH significant"] == 1) {
          Cladesummary[j,] <- c(Test[testidx,"Lk-ratio"], len, Clade$w2[i], Clade$w3[i], prop, prop*len,
                                Test[testidx,"P-val"], Test[testidx,"BH P-val"])
          rows[j]       <- rownames(Test)[testidx]
          j <- j + 1
      }                             
  }
  rownames(Cladesummary) <- rows
  colnames(Cladesummary) <- c('Lk-ratio', 'Length', 'w2', 'w3', 'p2', 'nrsites', 'P-val', 'BH P-val')
  return(Cladesummary)
}

# Like summarize M8 but with p-values and automatically extracting significant
extractM8 <- function(M8, Alignments, Test) {
  M8summary <- matrix(nrow=sum(Test[,"BH significant"]), ncol=7)
  rows      <- matrix(nrow=nrow(M8summary), ncol=1)
  j <- 1
  for (i in 1:nrow(M8)) {
      prop <- M8$p10[i]
      testidx <- which(rownames(Test) == M8[i,1])
      len  <- Alignments[which(Alignments[,1] == M8[i,1]),2]
      if (Test[testidx,"BH significant"] == 1) {
          M8summary[j,] <- c(Test[testidx,"Lk-ratio"], len, M8$w10[i], prop, prop*len, 
                            Test[testidx,"P-val"], Test[testidx,"BH P-val"])
          rows[j]       <- rownames(Test)[testidx]
          j <- j + 1
      }
  }
  rownames(M8summary) <- rows
  colnames(M8summary) <- c('Lk-ratio', 'Length', 'w10', 'p10', 'nrsites', 'P-val', 'BH P-val')
  return(M8summary)
}







