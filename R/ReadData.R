 # Functions for reading data for a run
# for 4 species and 5 species cases2
# More functions for doing some plots/tables
  



splitGroupnames <- function(names,groups) {  
  
  newnames <- matrix(nrow = length(names), ncol=1)  
  for (i in 1:length(names)) {      
      namesplit <- strsplit(names[i], '_')      
      
      if (length(namesplit[[1]]) > 1) {
          
          orthodbid <- strsplit(namesplit[[1]][2], ' ')[[1]][1]
          gid <- which(groups[,1] == orthodbid)
                            
          newnames[i] <- paste(namesplit[[1]][2], '&',groups[gid,4], '&', groups[gid,3])
      } else
          newnames[i] <- names[i]
  }  
  return(newnames)
}

getOrthoDB <- function(names) {
  names    <- as.character(names)
  newnames <- matrix(nrow = length(names), ncol=1)
  for (i in 1:length(names)) {
    namesplit <- strsplit(strsplit(names[i],' ')[[1]][1], '_')
    newnames[i] <- namesplit[[1]][2]    
  }
  return(newnames) 
}

getGroupInfo <- function(names, groups, field) {
  
  col <- switch(field, 
               orthodb    = 1,
               shortclass = 2,
               class      = 3,
               name       = 4,
               dmel       = 5,
               amell      = 6)
  
  names    <- as.character(names)
  newnames <- matrix(nrow = length(names), ncol=1)
  for (i in 1:length(names)) {
    namesplit <- strsplit(strsplit(names[i],' ')[[1]][1], '_')
    
    orthodbid <- strsplit(namesplit[[1]][2], ' ')[[1]][1]
    gid <- which(groups[,1] == orthodbid)  
    
    newnames[i] <- groups[gid,col]
  }
  return(newnames) 
}

# Expand OrthoDB ids to Id, Gene name, Classification
# Assumes ids are in the first column
expandLabels <- function(table, groups) {
  
  orthodb <- getGroupInfo(table[,1], groups, "orthodb")
  names   <- getGroupInfo(table[,1], groups, "name")
  class   <- getGroupInfo(table[,1], groups, "class")
  
  newtable <- cbind(orthodb, names, class, table[,-1])
  colnames(newtable) <- c("OrthoDB group", "Gene", "Classification", colnames(table[-1]))
  return(newtable)
}



getBEB <- function(basepath, run, model, alignment) {
  return(read.table(paste(basepath,"PAML_",run,"_Results/BEB/",alignment,'.',model,'.BEB', sep=''), header=TRUE)) 
}


# Print sites under positive selection and AA's from AMELL
getPositiveSites <- function(path, run, model, alignments, limit1=0.95, limit2=0.99, markup=FALSE) {
  
  if (ncol(alignments) > 1 && nrow(alignments) > 1) {
    modeldata <- alignments
    alignments <- rownames(modeldata)
  }
  
  possites <- matrix(nrow=length(alignments), ncol=1)  
  nrsites  <- matrix(nrow=length(alignments), ncol=1)
  for (i in 1:length(alignments)) {
    BEB <- getBEB(path, run, model, alignments[i])
    
    if (model == 'M8') {
      sites1 <- which(BEB$p10 > limit1)              
      sites2 <- which(BEB$p10 > limit2)              
    } else
      if (substr(model,1,2) == 'BS') {
        sites1 <- which((BEB$p2+BEB$p3) > limit1)
        sites2 <- which((BEB$p2+BEB$p3) > limit2)
      }
    else 
      sites1 <- sites2 <- c()          
    
    sites <- as.character(sites1)
    nrsites[i] <- length(sites1)
    if (length(sites1) > 0) {
      
      # Add Amino Acid
      for (j in 1:length(sites1))
        sites[j] <- paste(sites[j], BEB$AA[sites1[j]],sep='')
      
      if (markup==TRUE) {
          # Add underline for M8 if mean-var > 1.25
          if (model == 'M8') {          
            for (j in 1:length(sites1)) {
              if ((BEB$mean_w[sites1[j]]-BEB$std_w[sites1[j]]) > 1.25)
                sites[j] <- paste("[",sites[j],"]", sep='')
                #sites[j] <- paste("underline(",sites[j], ')', sep='')
            }
          }
          
          # Add bold font for limit2
          for (j in 1:length(sites1))
            if (sites1[j] %in% sites2) {
              sites[j] <- paste("<",sites[j],">", sep='')
              #sites[j] <- paste('bold(',sites[j],')', sep='')
            }
      }
      
      possites[i] <- paste(sites, collapse=', ')            
    } else
      possites[i] <- 'None'
    
  }
  rownames(possites) <- alignments
  colnames(possites) <- c("Positively selected sites")
  
  return(cbind(nrsites, possites))
}

getPositiveSelectionTable <- function(path, run, model, alignments, groups, limit1=0.95, limit2=0.99) {
  
  possites <- getPositiveSites(path, run, model, alignments, limit1, limit2, markup=TRUE)
  perm     <- order(as.numeric(possites[,1]), decreasing=TRUE)
  possel   <- cbind(expandLabels(rownames(alignments), groups)[perm,], alignments[perm,c(2,6,7)], possites[perm,2])
  
  return(possel)
}

# Print sites under positive selection and AA's from AMELL
getPositiveSitesLatex <- function(path, run, model, alignments, limit1=0.95, limit2=0.99, latex=FALSE, caption='', label='') {
  
  if (ncol(alignments) > 1 && nrow(alignments) > 1) {
      modeldata <- alignments
      alignments <- rownames(modeldata)
  }
  
  possites <- matrix(nrow=length(alignments), ncol=1)  
  nrsites  <- matrix(nrow=length(alignments), ncol=1)
  for (i in 1:length(alignments)) {
      BEB <- getBEB(path, run, model, alignments[i])
   
      if (model == 'M8') {
          sites1 <- which(BEB$p10 > limit1)              
          sites2 <- which(BEB$p10 > limit2)              
      } else
      if (substr(model,1,2) == 'BS') {
          sites1 <- which((BEB$p2+BEB$p3) > limit1)
          sites2 <- which((BEB$p2+BEB$p3) > limit2)
      }
      else 
          sites1 <- sites2 <- c()          
     
      sites <- as.character(sites1)
      nrsites[i] <- length(sites1)
      if (length(sites1) > 0) {
      
          # Add Amino Acid
          for (j in 1:length(sites1))
              sites[j] <- paste(sites[j], BEB$AA[sites1[j]],sep='')
                    
          # Add underline for M8 if mean-var > 1.25
          if (model == 'M8') {          
              for (j in 1:length(sites1)) {
                  if ((BEB$mean_w[sites1[j]]-BEB$std_w[sites1[j]]) > 1.25)
                      sites[j] <- paste("\\underline{",sites[j], '}', sep='')
              }
          }
          
          # Add bold font for limit2
          for (j in 1:length(sites1))
              if (sites1[j] %in% sites2) {
                  sites[j] <- paste('\\textbf{',sites[j],'}', sep='')
              }
          
          possites[i] <- paste(sites, collapse=', ')            
      } else
          possites[i] <- 'None'
      
  }
  rownames(possites) <- alignments
  colnames(possites) <- c("Positively selected sites")
  
  if (latex == TRUE && !is.null(modeldata)) {
    require(xtable)
    
    perm <- order(nrsites, decreasing=TRUE)
    
    table <- data.frame(modeldata[,c("Length", "Lk-ratio", "P-val", "BH P-val")], possites)    
    colnames(table) <- c("N", "Likelihood Ratio",  "$p$-value", "BH-corrected $p$-value", "Positively selected sites")

        
    return(xtable(table[perm,], align="lrrrrp{3cm}", digits=c(0,0,3,5,5,0), label=label, caption=caption))
    
    
  } else
      return(possites)  
}

plotM8BEB <- function(path, run, alignment, groups, limit1=0.95, limit2=0.99, outfile='') {
    
  old.par <- par()
  par(mar = c(5,5,4,2)+0.1)
  
  capt <- paste(strsplit(splitGroupnames(alignment, groups), '&')[[1]], collapse='-')
  
  BEB   <- getBEB(path, run, 'M8', alignment)
  high  <- BEB$mean_w+BEB$std_w
  low   <- BEB$mean_w-BEB$std_w
  ymax  <- 5 #max(high)+0.25
  low[which(low < 0)] <- 0
  plot(1, type="n", xlab='Codon', ylab=expression(E(omega)), ylim=c(0,ymax), xlim=c(0,length(BEB$Site)*1.12), main=capt)
    
  polygon(c(0, length(BEB$Site), length(BEB$Site), 0), c(1.25, 1.25, 0.75, 0.75), col=pal.trans[6], border=NA)  
  #polygon(c(0, length(BEB$Site), length(BEB$Site), 0), c(ymax, ymax, 1.25, 1.25), col=paste(pal.light[2], '88', sep= ''), border=NA)
  #polygon(c(0, length(BEB$Site), length(BEB$Site), 0), c(0.75, 0.75, 0, 0), col=paste(pal.light[3], '88', sep= ''), border=NA)
  
  polygon(c(BEB$Site, rev(BEB$Site)), c(high, rev(low)), col=pal.trans[1], border=NA)
  lines( BEB$Site,BEB$mean_w, las=1, col=pal.dark[1], lwd=2)
  points(BEB$Site,BEB$mean_w, las=1, col=pal.dark[1], pch=20)
  
  #lines(BEB$Site, rep(1,length(BEB$Site)),col=pal.light[3], lty=2, lwd=2)
  lines(BEB$Site, rep(1.25,length(BEB$Site)),col=pal.light[3], lty=5, lwd=2)
  lines(BEB$Site, rep(0.75,length(BEB$Site)),col=pal.light[3], lty=5, lwd=2)
  
  # Sites where posterior mean - variance is above 1.25
  possites3 <- which(low > 1.25)
  points(possites3, rep(ymax-0.30, length(possites3)), pch=18, col=pal.light[5])
  
  # Sites with Pr > limit1 to be positively selected
  possites1 <- which(BEB$p10 > limit1)
  points(possites1, rep(ymax-0.15, length(possites1)), pch=18, col=pal.light[2])
  
  # Sites with Pr > limit2 to be positively selected
  possites2 <- which(BEB$p10 > limit2)
  points(possites2, rep(ymax, length(possites2)), pch=18, col=pal.light[3])
  
  lines(BEB$Site, rep(ymax-0.30, length(BEB$Site)), col=pal.trans[5], lty=1, lwd=1)
  lines(BEB$Site, rep(ymax-0.15, length(BEB$Site)), col=pal.trans[2], lty=1, lwd=1)
  lines(BEB$Site, rep(ymax, length(BEB$Site)), col=pal.trans[3], lty=1, lwd=1)
  text(length(BEB$Site),ymax-0.30, substitute(E(omega) - sqrt(Var(omega)) > 1.25), pos=4)
  text(length(BEB$Site),ymax-0.15, substitute(P(omega > 1) > list(x), list(x=limit1) ), pos=4)
  text(length(BEB$Site),ymax,      substitute(P(omega > 1) > list(x), list(x=limit2) ), pos=4)
  
  
  if (outfile != '') {
    dev.copy(pdf, outfile, height=10, width=20)
    #dev.copy(png, outfile, height=1000, width=2000)    
    dev.off()    
  }
  
  par(old.par)
  
}



plotCladePoints <- function(Clade, M8, BS, scaling=1, palidx=1, labelfun=NULL) {
  
  #print(Clade)
  if (nrow(Clade) > 0) {
      
      Clade_M8 <- Clade[intersect(rownames(Clade), rownames(M8)),,drop=FALSE]
      Clade_BS <- Clade[intersect(rownames(Clade), rownames(BS)),,drop=FALSE]
              
      for (i in 1:nrow(Clade)) {
        draw.circle(Clade[i,"w2"], Clade[i,"w3"], scaling*Clade[i,"p2"]/4,
                    border=pal.dark[palidx], col=pal.trans[palidx], lty=1)
      }  
      points(Clade[,"w2"], Clade[,"w3"], pch=19, col=pal.dark[palidx], cex=0.25)
      points(Clade_M8[,"w2"], Clade_M8[,"w3"], pch=3, col=pal.dark[5], lwd=1.5, cex=1.5)
      points(Clade_BS[,"w2"], Clade_BS[,"w3"], pch=4, col=pal.dark[3], lwd=1.5, cex=1.5)  
      
      if (!is.null(labelfun)) {
        labels <- labelfun(rownames(Clade))
        
        text(Clade[,"w2"], Clade[,"w3"], labels, offset=0.25, cex=0.5, pos=4)   
      }
  }
}



# Plot differences in selection between clades
#   pattern highlights all rows containing the pattern in rowname
#   split is where the zoomed part cuts off
#   Clade 1 - w2
#   Clade 2 - w3
plotClade <- function(Clade, M8, BS, M3 = NULL, clades=c("Clade 1", "Clade 2"), pattern='', outfile='', split=2, palidx=1, size=10, labelfun=NULL) {  
  
  if (!is.null(M3)) {
    M3 <- cbind(M3, M3[,"w2"])
    colnames(M3)[length(colnames(M3))] <- "w3"
  }
  
  patsearch <- pattern == ''
  
  #old.par <- par(mar = c(0, 0, 0, 0))
  old.par <- par()
  par(mai = c(1,1,2,2))
  plot(1, type='n', xlim=c(0,split),ylim=c(0,split), asp=1, las=1, 
       xlab=paste("dN/dS on the ",clades[1],"clade"), ylab=paste("dN/dS on the",clades[2], 'clade'))
  if (!is.null(M3)) {
    plotCladePoints(M3[grep(pattern, rownames(M3), invert=!patsearch),,drop=FALSE], M8, BS, palidx=6, labelfun=labelfun)  
    plotCladePoints(M3[grep(pattern, rownames(M3), invert= patsearch),,drop=FALSE], M8, BS, palidx=palidx[2], labelfun=labelfun)  
  }
  
  #plotCladePoints(Clade[which(Clade[,"w2"] <= split | Clade[,"w3"] <= split),], M8, BS, palidx=palidx)
  plotCladePoints(Clade[grep(pattern, rownames(Clade), invert=!patsearch),,drop=FALSE], M8, BS, palidx=palidx[1], labelfun=labelfun)
  plotCladePoints(Clade[grep(pattern, rownames(Clade), invert= patsearch),,drop=FALSE], M8, BS, palidx=palidx[2], labelfun=labelfun)              
  abline(a=0,b=1,col=pal.dark[3],lty=2, lwd=1.5)  
  
  par(xpd=TRUE)
  legend(x = 'bottomright', inset=c(-0.2,0), pch=c(3,4), pt.lwd=2, lty=0, pt.cex=1.5, 
         col=c(pal.dark[5], pal.dark[3]), bty='n', title.adj=0, 
         title="       Evidence\n       for positive\n       selection:", legend=c("On all\nbranches\n", paste("Between\n",clades[1],"\nand ",clades[2], sep='')))
  if (pattern == '') {
    legend(x = 'topleft', inset=c(0,-0.15), pch=21, pt.cex=2, pt.bg=c(pal.trans[palidx[1]], pal.trans[6]),
           col=c(pal.dark[palidx[1]], pal.dark[6]), horiz=TRUE, bty='n', title.adj=0,
           title=paste("Difference in evolutionary pressures\nbetween",clades[1],"and",clades[2],"clades:"), 
           legend = c("Significant", "Not significant"))
  } else {
    legend(x = 'topleft', inset=c(0,-0.15), pch=21, pt.cex=2, pt.bg=c(pal.trans[palidx[1]], pal.trans[6], pal.trans[palidx[2]]),
           col=c(pal.dark[palidx[1]], pal.dark[6], pal.dark[palidx[2]]), horiz=FALSE, bty='n', title.adj=0,
           title=paste("Difference in evolutionary pressures\nbetween",clades[1],"and",clades[2],"clades:"), 
           legend = c("Significant", "Not significant", pattern))    
  }
  par(xpd=FALSE)
  
  # Set up inset plot and clear area
  par(mai = c(size/3,size/3,2,2), new=TRUE)
  plot(1, type='n', xlim=c(0,12),ylim=c(0,12), xlab='', ylab='',
       xaxt="n", yaxt="n", asp=1)
  axis(1,at=c(0:10,12), labels=c(0:10, expression(infinity)), las=1)
  axis(2,at=c(0:10,12), labels=c(0:10, expression(infinity)), las=1)
  polygon(c(-1,-1, 13, 13), c(-1, 13, 13, -1), col='white', lty=0, border=NA)
  par(new = TRUE)
  plot(1, type='n', xlim=c(0,12),ylim=c(0,12), xlab='', ylab='',
       xaxt="n", yaxt="n", asp=1)
  
  if (!is.null(M3)) {
    M3[which(M3[,"w2"] > 10),"w2"] <- 12
    M3[which(M3[,"w3"] > 10),"w3"] <- 12    
    plotCladePoints(M3[grep(pattern, rownames(M3), invert=!patsearch),,drop=FALSE], scaling=12.0/split, M8, BS, palidx=6, labelfun=labelfun)  
    plotCladePoints(M3[grep(pattern, rownames(M3), invert= patsearch),,drop=FALSE], scaling=12.0/split, M8, BS, palidx=palidx[2], labelfun=labelfun)  
  }
  
  Clade[which(Clade[,"w2"] > 10),"w2"] <- 12
  Clade[which(Clade[,"w3"] > 10),"w3"] <- 12    
  plotCladePoints(Clade[grep(pattern, rownames(Clade), invert=!patsearch),,drop=FALSE], scaling=12.0/split, M8, BS, palidx=palidx[1], labelfun=labelfun)
  plotCladePoints(Clade[grep(pattern, rownames(Clade), invert= patsearch),,drop=FALSE], M8, BS, scaling=12.0/split, palidx=palidx[2], labelfun=labelfun)    
  
  abline(a=0,b=1,col=pal.dark[3],lty=2, lwd=1.5)  
  #polygon(c(0,0,split,split), c(0,split,split,0), col=NA, border=pal.dark[3], lwd=1.5, lty=2)  
  lines(c(-1,split, split, split), c(split, split, split, -1), col=pal.dark[3], lwd=1.5, lty=2)
    
  if (outfile != '') {
    dev.copy(pdf, outfile, height=size, width=size)       
    dev.off()    
  }
  
  par(old.par)
}






# Deprecated
plotCladeInverted <- function(Clade, M8, BS, pattern='', outfile='', split=2) {  
  
  
  par(mai = c(1,1,1,1))
  plot(1, type='n', xlim=c(0,12),ylim=c(0,12), xlab="dN/dS (Bombus)", ylab="dN/dS (Apis)",
       xaxt="n", yaxt="n", asp=1)
  axis(1,at=c(0:10,12), labels=c(0:10, expression(infinity)))
  axis(2,at=c(0:10,12), labels=c(0:10, expression(infinity)))
  abline(a=0,b=1,col=pal.dark[3],lty=2, lwd=1.5)  
  #polygon(c(0,0,split,split), c(0,split,split,0), col=pal.light[3], border=pal.dark[3])  
  
  Clade[which(Clade[,"w2"] > 10),"w2"] <- 12
  Clade[which(Clade[,"w3"] > 10),"w3"] <- 12    
  plotCladePoints(Clade[which(Clade[,"w2"] > 0 | Clade[,"w3"] > 0),], M8, BS, scaling=1)  
  
  
  par(mai = c(4,4,1,1))
  par(new = TRUE)
  plot(1, type='n', xlim=c(0,split),ylim=c(0,split), xlab='', ylab='', asp=1)
  abline(a=0,b=1,col=pal.dark[3],lty=2, lwd=1.5)  
  plotCladePoints(Clade[which(Clade[,"w2"] <= split | Clade[,"w3"] <= split),], M8, BS)
  
  legend(x = 'topright', inset=c(0.05,0.1), pch=c(3,4), pt.lwd=1.5, lty=0, pt.cex=1.5, col=pal.dark[1],
         legend=c("Complete phylogeny", "Connecting branch"))    
  
  #idx1 <- grep(pattern, rownames(Clade))
  #idx2 <- grep(pattern, rownames(Clade), invert=TRUE)
  #points(Clade[idx1,"w2"], Clade[idx1,"w3"], pch=19, col=pal.trans[3])
  
  if (outfile != '') {
    dev.copy(pdf, outfile, height=10, width=10)    
    dev.off()    
  }
}


read4Species <- function(basepath, run, alpha) {
  
    resultpath <- paste(basepath,"PAML_",run,"_Results/",sep="")  
    
    #############################################################################################
    # Read data from PAML
      M0  <- read.table(paste(resultpath,'Sites-Best.M0.w',sep=''),header=TRUE)
      M1a <- read.table(paste(resultpath,'Sites-Best.M1.w',sep=''),header=TRUE)
      M2a <- read.table(paste(resultpath,'Sites-Best.M2.w',sep=''),header=TRUE)
      M3  <- read.table(paste(resultpath,'Sites-Best.M3.w',sep=''),header=TRUE)
      M7  <- read.table(paste(resultpath,'Sites-Best.M7.w',sep=''),header=TRUE)
      M8  <- read.table(paste(resultpath,'Sites-Best.M8.w',sep=''),header=TRUE)
      M8a <- read.table(paste(resultpath,'Sites.M8A-1.000.w',sep=''),header=TRUE)
      
      BSAlt  <- read.table(paste(resultpath,'BranchSiteAlt-Best.Tree1.w',sep=''),header=TRUE)
      BSNull <- read.table(paste(resultpath,'BranchSiteNull-1.000.Tree1.w',sep=''),header=TRUE)
      CladeC <- read.table(paste(resultpath,'CladeC-Best.Tree1.w',sep=''),header=TRUE)
      CladeD <- read.table(paste(resultpath,'CladeD-Best.Tree1.w',sep=''),header=TRUE)
    
    # Read alignment data
      alignments <- read.table(paste(resultpath,"/Alignments.csv",sep=''),header=TRUE)
      groups     <- read.table(paste(basepath, 'Groups_Final.csv',sep=''), header=TRUE, sep=',', as.is=TRUE)

    
    
    #############################################################################################
    # Perform tests
    
    # Site models
      cat("M3 vs M0: Look for site-to-site variation in omega\n")
      M3vsM0 <- doLRT(M3, M0, 4, alpha)
      
      cat("M2a vs M1a: Conservative test for positive selection\n")
      M2avsM1a <- doLRT(M2a, M1a, 2, alpha)
      
      cat("M8 vs M7: Less conservative test for positive selection\n")
      M8vsM7 <- doLRT(M8, M7, 2, alpha)
      
      cat("M8 vs M8a: Test for positive selection (check fit of beta distribution)\n")
      M8vsM8a <- doLRT(M8, M8a, c(0,2), alpha)
  
  # Branch-site model
      cat("Branch-site test for positive selection (branch between Bombus and Apis under selection)\n")
      BStest <- doLRT(BSAlt, BSNull, c(0,1), alpha)
  
  # Clade model
      cat("Clade model C (Bombus and Apis in different clades)\n")
      CladetestC <- doLRT(CladeC,M1a, 3, alpha)
  
      cat("Clade model D (Bombus and Apis in different clades)\n")
      CladetestD <- doLRT(CladeD,M3, 1, alpha)

  #############################################################################################
  # Build data structure
  
  models <- list(M0=M0, M1a=M1a, M2a=M2a, M3=M3, M7=M7, M8=M8, M8a=M8a, 
                 BSAlt=BSAlt, BSNull=BSNull, CladeC=CladeC, CladeD=CladeD)
  tests  <- list(M3vsM0=M3vsM0, M2avsM1a=M2avsM1a, M8vsM7=M8vsM7, M8vsM8a=M8vsM8a,
                 BStest=BStest, CladetestC=CladetestC, CladetestD=CladetestD)
    
  return(list(models=models, tests=tests, alignments=alignments, groups=groups))
}



read5Species <- function(basepath, run, alpha) {
  
  resultpath <- paste(basepath,"PAML_",run,"_Results/",sep="")
  
  #############################################################################################
  # Read data from PAML
      M0  <- read.table(paste(resultpath,'Sites-Best.M0.w',sep=''),header=TRUE)
      M1a <- read.table(paste(resultpath,'Sites-Best.M1.w',sep=''),header=TRUE)
      M2a <- read.table(paste(resultpath,'Sites-Best.M2.w',sep=''),header=TRUE)
      M3  <- read.table(paste(resultpath,'Sites-Best.M3.w',sep=''),header=TRUE)
      M7  <- read.table(paste(resultpath,'Sites-Best.M7.w',sep=''),header=TRUE)
      M8  <- read.table(paste(resultpath,'Sites-Best.M8.w',sep=''),header=TRUE)
      M8a <- read.table(paste(resultpath,'Sites.M8A-1.000.w',sep=''),header=TRUE)
      
      BSAltTree1  <- read.table(paste(resultpath,'BranchSiteAlt-Best.Tree1.w',sep=''),header=TRUE)
      BSNullTree1 <- read.table(paste(resultpath,'BranchSiteNull-1.000.Tree1.w',sep=''),header=TRUE)
      BSAltTree2  <- read.table(paste(resultpath,'BranchSiteAlt-Best.Tree2.w',sep=''),header=TRUE)
      BSNullTree2 <- read.table(paste(resultpath,'BranchSiteNull-1.000.Tree2.w',sep=''),header=TRUE)
      BSAltTree3  <- read.table(paste(resultpath,'BranchSiteAlt-Best.Tree3.w',sep=''),header=TRUE)
      BSNullTree3 <- read.table(paste(resultpath,'BranchSiteNull-1.000.Tree3.w',sep=''),header=TRUE)
      BSAltTree4  <- read.table(paste(resultpath,'BranchSiteAlt-Best.Tree4.w',sep=''),header=TRUE)
      BSNullTree4 <- read.table(paste(resultpath,'BranchSiteNull-1.000.Tree4.w',sep=''),header=TRUE)
      BSAltTree5  <- read.table(paste(resultpath,'BranchSiteAlt-Best.Tree5.w',sep=''),header=TRUE)
      BSNullTree5 <- read.table(paste(resultpath,'BranchSiteNull-1.000.Tree5.w',sep=''),header=TRUE)
  
      CladeCTree1  <- read.table(paste(resultpath,'CladeC-Best.Tree1.w',sep=''),header=TRUE)
      CladeCTree2  <- read.table(paste(resultpath,'CladeC-Best.Tree2.w',sep=''),header=TRUE)
      CladeDTree1  <- read.table(paste(resultpath,'CladeD-Best.Tree1.w',sep=''),header=TRUE)
      CladeDTree2  <- read.table(paste(resultpath,'CladeD-Best.Tree2.w',sep=''),header=TRUE)
  
  # Read alignment data
      alignments <- read.table(paste(resultpath,"/Alignments.csv",sep=''),header=TRUE)
      groups     <- read.table(paste(basepath, 'Groups_Final.csv',sep=''),header=TRUE, sep=',', as.is=TRUE)
  
  
  #############################################################################################
  # Perform tests
  
  # Site models
      cat("M3 vs M0: Look for site-to-site variation in omega\n")
      M3vsM0 <- doLRT(M3, M0, 4, alpha)
      
      cat("M2a vs M1a: Conservative test for positive selection\n")
      M2avsM1a <- doLRT(M2a, M1a, 2, alpha)
      
      cat("M8 vs M7: Less conservative test for positive selection\n")
      M8vsM7 <- doLRT(M8, M7, 2, alpha)
      
      cat("M8 vs M8a: Test for positive selection (check fit of beta distribution)\n")
      M8vsM8a <- doLRT(M8, M8a, c(0,2), alpha)
  
  # Branch-site model
      cat("Branch-site test for positive selection (Tree 1: All 3 internal branches under positive selection)\n")
      BStest1 <- doLRT(BSAltTree1, BSNullTree1, c(0,1), alpha)
      
      cat("Branch-site test for positive selection (Tree 2: 2 branches between Bombus and Apis under positive selection)\n")
      BStest2 <- doLRT(BSAltTree2, BSNullTree2, c(0,1), alpha)
      
      cat("Branch-site test for positive selection (Tree 3: Branch to Megachile under positive selection)\n")
      BStest3 <- doLRT(BSAltTree3, BSNullTree3, c(0,1), alpha)
      
      cat("Branch-site test for positive selection (Tree 4: Branch to Bombus under positive selection)\n")
      BStest4 <- doLRT(BSAltTree4, BSNullTree4, c(0,1), alpha)
  
      cat("Branch-site test for positive selection (Tree 5: Branch to Apis under positive selection)\n")
      BStest5 <- doLRT(BSAltTree5, BSNullTree5, c(0,1), alpha)
  
  
  # Clade model C
      cat("Clade model C (2 clades compared to no clades - social/solitary)\n")
      CladeCtest1 <- doLRT(CladeCTree1,M1a, 3, alpha)
      
      cat("Clade model C (3 clades compared to no clades)\n")
      CladeCtest2 <- doLRT(CladeCTree2,M1a, 4, alpha)
      
      cat("Clade model C (3 clades compared to 2 clades)\n")
      CladeCtest3 <- doLRT(CladeCTree2, CladeCTree1, 1, alpha)
  
  
  # Clade model D
      cat("Clade model D (2 clades compared to no clades - social/solitary)\n")
      CladeDtest1 <- doLRT(CladeDTree1,M3, 1, alpha)
      
      cat("Clade model D (3 clades compared to no clades)\n")
      CladeDtest2 <- doLRT(CladeDTree2,M3, 2, alpha)
      
      cat("Clade model D (3 clades compared to 2 clades)\n")
      CladeDtest3 <- doLRT(CladeDTree2, CladeDTree1, 1, alpha)
  
  
  #############################################################################################
  # Build data structure
  
  models <- list(M0=M0, M1a=M1a, M2a=M2a, M3=M3, M7=M7, M8=M8, M8a=M8a, 
                 BSAltTree1=BSAltTree1, BSAltTree2=BSAltTree2, BSAltTree3=BSAltTree3, BSAltTree4=BSAltTree4, BSAltTree5=BSAltTree5, 
                 BSNullTree1=BSNullTree1, BSNullTree2=BSNullTree2, BSNullTree3=BSNullTree3, BSNullTree4=BSNullTree4, BSNullTree5=BSNullTree5,
                 CladeCTree1=CladeCTree1, CladeCTree2=CladeCTree2, CladeDTree1=CladeDTree1, CladeDTree2=CladeDTree2)
  tests  <- list(M3vsM0=M3vsM0, M2avsM1a=M2avsM1a, M8vsM7=M8vsM7, M8vsM8a=M8vsM8a,
                 BStest1=BStest1, BStest2=BStest2, BStest3=BStest3, BStest4=BStest4, BStest5=BStest5,
                 CladeCtest1=CladeCtest1, CladeCtest2=CladeCtest2, CladeCtest3=CladeCtest3,
                 CladeDtest1=CladeDtest1, CladeDtest2=CladeDtest2, CladeDtest3=CladeDtest3)
    
  return(list(models=models, tests=tests, alignments=alignments, groups=groups))
}

