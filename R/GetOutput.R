#############################################################################################
# Read data and do tests

cat('Reading 4 species case\n')
data.4species <- read4Species(path4, run, alpha)

cat('\nReading 5 species case\n')
data.5species <- read5Species(path5, run, alpha)

dir.create(outpath, showWarnings=FALSE)
dir.create(paste(outpath,'Venn/'  ,sep=''), showWarnings=FALSE)
dir.create(paste(outpath,'M8_4taxa/',sep=''), showWarnings=FALSE)
dir.create(paste(outpath,'M8_5taxa/',sep=''), showWarnings=FALSE)
dir.create(paste(outpath,'CladeD_Bombus_Apis/',sep=''), showWarnings=FALSE)
dir.create(paste(outpath,'CladeD_Social_Non/' ,sep=''), showWarnings=FALSE)
#dir.create(paste(outpath,'CladeC_Bombus_Apis/',sep=''), showWarnings=FALSE)
#dir.create(paste(outpath,'CladeC_Social_Non/' ,sep=''), showWarnings=FALSE)
dir.create(paste(outpath,'Tables_4taxa/', sep=''), showWarnings=FALSE)
dir.create(paste(outpath,'Tables_5taxa/', sep=''), showWarnings=FALSE)



#############################################################################################
# Patterns across the whole phylogeny

# Summary statistics of w of M0 model across all sites and branches
  M0.w.summary = matrix(nrow=5, ncol=2)
  colnames(M0.w.summary) <- c("4 species Global $\\omega$", "5 species Global $\\omega$")
  rownames(M0.w.summary) <- c("Mean", "Median","Variance", "Standard deviation", "Standard error")
  M0.w.summary[1,] <- c(mean(data.4species$models$M0$w0),      mean(data.5species$models$M0$w0))
  M0.w.summary[2,] <- c(median(data.4species$models$M0$w0),    median(data.5species$models$M0$w0))
  M0.w.summary[3,] <- c(var(data.4species$models$M0$w0),       var(data.5species$models$M0$w0))
  M0.w.summary[4,] <- c(sqrt(var(data.4species$models$M0$w0)), sqrt(var(data.5species$models$M0$w0)))
  M0.w.summary[5,] <- c(sqrt(var(data.4species$models$M0$w0)/length(data.4species$models$M0$w0)), 
                        sqrt(var(data.5species$models$M0$w0)/length(data.5species$models$M0$w0)))

  M0.dNlen.summary <- matrix(nrow=5, ncol=2)
  colnames(M0.dNlen.summary) <- c("4 species Tree Length (dN)", "5 species Tree Length (dN)")
  rownames(M0.dNlen.summary) <- c("Mean", "Median","Variance", "Standard deviation", "Standard error")
  M0.dNlen.summary[1,] <- c(mean(data.4species$models$M0$dNlen),      mean(data.5species$models$M0$dNlen))
  M0.dNlen.summary[2,] <- c(median(data.4species$models$M0$dNlen),    median(data.5species$models$M0$dNlen))
  M0.dNlen.summary[3,] <- c(var(data.4species$models$M0$dNlen),       var(data.5species$models$M0$dNlen))
  M0.dNlen.summary[4,] <- c(sqrt(var(data.4species$models$M0$dNlen)), sqrt(var(data.5species$models$M0$dNlen)))
  M0.dNlen.summary[5,] <- c(sqrt(var(data.4species$models$M0$dNlen)/length(data.4species$models$M0$dNlen)),  
                            sqrt(var(data.5species$models$M0$dNlen)/length(data.5species$models$M0$dNlen)))

  M0.dSlen.summary <- matrix(nrow=5, ncol=2)
  colnames(M0.dSlen.summary) <- c("4 species Tree Length (dS)", "5 species Tree Length (dS)")
  rownames(M0.dSlen.summary) <- c("Mean", "Median","Variance", "Standard deviation", "Standard error")
  M0.dSlen.summary[1,] <- c(mean(data.4species$models$M0$dSlen),      mean(data.5species$models$M0$dSlen))
  M0.dSlen.summary[2,] <- c(median(data.4species$models$M0$dSlen),    median(data.5species$models$M0$dSlen))
  M0.dSlen.summary[3,] <- c(var(data.4species$models$M0$dSlen),       var(data.5species$models$M0$dSlen))
  M0.dSlen.summary[4,] <- c(sqrt(var(data.4species$models$M0$dSlen)), sqrt(var(data.5species$models$M0$dSlen)))
  M0.dSlen.summary[5,] <- c(sqrt(var(data.4species$models$M0$dSlen)/length(data.4species$models$M0$dSlen)), 
                            sqrt(var(data.5species$models$M0$dSlen)/length(data.5species$models$M0$dSlen)))

  M0.title <- expression("Global"~omega, "Tree length (dN)", "Tree length (dS)")
  M0.4species.table <-  makeTable(cbind(M0.w.summary[,1], M0.dNlen.summary[,1], M0.dSlen.summary[,1]), 
                                  caption=expression("Statistics for the global"~omega~"ratio obtained by the M0 model (4 taxa tree)."), footnotes=c(" "),
                                  title=M0.title, col=1, rownames=TRUE, equalwidth=TRUE, floatcol=c(1,2,3), draw=FALSE)

  M0.title <- expression("Global"~omega, "Tree length (dN)", "Tree length (dS)")
  M0.5species.table <-  makeTable(cbind(M0.w.summary[,2], M0.dNlen.summary[,2], M0.dSlen.summary[,2]), 
                                  caption=expression("Statistics for the global"~omega~"ratio obtained by the M0 model (5 taxa tree)."), footnotes=c(" "),
                                  title=M0.title, col=1, rownames=TRUE, equalwidth=TRUE, floatcol=c(1,2,3), draw=FALSE)


# Slowest and Fastest 10 genes (M0 model)
  M0.slowest.4species <- data.4species$models$M0[order(data.4species$models$M0$w0)[1:30],c(1,3,6)]
  M0.slowest.5species <- data.5species$models$M0[order(data.5species$models$M0$w0)[1:30],c(1,3,6)]
  M0.fastest.4species <- data.4species$models$M0[order(data.4species$models$M0$w0, decreasing=TRUE)[1:30],c(1,3,6)]
  M0.fastest.5species <- data.5species$models$M0[order(data.5species$models$M0$w0, decreasing=TRUE)[1:30],c(1,3,6)]

  slowestfastest.title  <- expression("OrthoDB group"^a, "Gene"^b~"                              ", "Classification                                                     ", "Global"~omega^c,"Tree Length"^d~" ")
  slowestfastest.foot   <- expression(""^a~"Group identifiers are from OrthoDB 6 (http://cegg.unige.ch/orthodb6).",
                                      ""^b~"Unless otherwise specified, gene names are taken from the"~italic("A. mellifera")~"or"~italic("D. melanogaster")~"orthologs.", 
                                      ""^c~"Maximum likelihod estimate across all sites and branches.", 
                                      ""^d~"Tree length in synonymous substitutions per synonymous sites.")

  M0.slowest.4species.table <- makeTable(expandLabels(M0.slowest.4species, data.4species$groups),             
                                                      caption = "The 30 slowest evolving genes on the 4 taxa tree as determined by the M0 model",
                                                      title = slowestfastest.title, footnotes = slowestfastest.foot,  col=1, draw=FALSE)  
  
  M0.fastest.4species.table <- makeTable(expandLabels(M0.fastest.4species, data.4species$groups),             
                                                      caption = "The 30 fastest evolving genes on the 4 taxa tree as determined by the M0 model",
                                                      title = slowestfastest.title, footnotes = slowestfastest.foot,  col=4, draw=FALSE)
  
  M0.slowest.5species.table <- makeTable(expandLabels(M0.slowest.5species, data.5species$groups),             
                                                      caption = "The 30 slowest evolving genes on the 5 taxa tree as determined by the M0 model",
                                                      title = slowestfastest.title, footnotes = slowestfastest.foot,  col=1, draw=FALSE)
  
  M0.fastest.5species.table <- makeTable(expandLabels(M0.fastest.5species, data.5species$groups),             
                                                      caption = "The 30 fastest evolving genes on the 5 taxa tree as determined by the M0 model",
                                                      title = slowestfastest.title, footnotes = slowestfastest.foot,  col=4, draw=FALSE)  
  

# Positive selection on whole gene
  pos.4species <- getSignificant(data.4species$tests$M2avsM1a) | 
                  getSignificant(data.4species$tests$M8vsM7)   | 
                  getSignificant(data.4species$tests$M8vsM8a)

  if (!identical(pos.4species, getSignificant(data.4species$tests$M8vsM7))) 
    cat("M8 vs M7 does not encompass other tests for positive selection (4 species)\n")
  
  pos.5species <- getSignificant(data.5species$tests$M2avsM1a) | 
                  getSignificant(data.5species$tests$M8vsM7)   | 
                  getSignificant(data.5species$tests$M8vsM8a)
  
  if (!identical(pos.5species, getSignificant(data.5species$tests$M8vsM7))) 
    cat("M8 vs M7 does not encompass other tests for positive selection (5 species)\n")
  
  M8.4species.sig <- extractM8(data.4species$models$M8, data.4species$alignments, data.4species$tests$M8vsM7)
  M8.5species.sig <- extractM8(data.5species$models$M8, data.5species$alignments, data.5species$tests$M8vsM7)

  possel.title <- expression("OrthoDB Group"^a~" ", "Gene"^b~"                              ", "Classification                                                     ", 
                             "Sites"^c, italic(p)*"-value"^d, italic(q)*"-value"^e, "Positively selected sites"^f~"                                                           ")
  possel.foot  <- expression(""^a~"Group identifiers are from OrthoDB 6 (http://cegg.unige.ch/orthodb6).",
                             ""^b~"Unless otherwise specified, gene names are taken from the"~italic("A. mellifera")~"or"~italic("D. melanogaster")~"orthologs.", 
                             ""^c~"Total number of codons in the alignment after trimming with Gblocks.", 
                             ""^d~"From LRT comparing model M7 to M8.",
                             ""^e~"Multiple test correction by the method of Benjamini and Hochberg to control the false discovery rate.",
                             ""^f~"Sites are classified as under positive selection if the Bayesian posterior probability > 0.75 (> 0.95 in bold). Sites where"~E[omega]-sqrt(Var(omega)) > 1.25~"are in orange.",
                             "  Reference sequence taken from"~italic("A. mellifera"))

  M8.4species.sig.table <- makeTable(getPositiveSelectionTable(path4, run, 'M8', M8.4species.sig, data.4species$groups, 0.75, 0.95), 
                                                              caption = "Genes under positive selection (using FDR < 0.05) across the whole phylogeny (4 taxa tree)", 
                                                              title   = possel.title, footnotes = possel.foot,            
                                                              col=5, floatcols=c(5,6), justcol=7, maxwidth=50, draw=FALSE)

  M8.5species.sig.table <- makeTable(getPositiveSelectionTable(path5, run, 'M8', M8.5species.sig, data.5species$groups, 0.75, 0.95), 
                                                              caption = "Genes under positive selection (using FDR < 0.05) across the whole phylogeny (5 taxa tree)", 
                                                              title   = possel.title, footnotes = possel.foot,            
                                                              col=5, floatcols=c(5,6), justcol=7, maxwidth=50, draw=FALSE)


  # Overlap between genes under positive selection in 4 and 5 taxa cases
  venn.diagram(list("5 Taxa"=rownames(M8.5species.sig), "4 Taxa"=rownames(M8.4species.sig)),
               filename=paste(outpath,'Venn/PositiveSelection_sig.tiff',sep='/'), col="black", lty=1, lwd=4, margin=0.15,
               fill=pal.dark[c(1,3)], alpha=0.50, cex=2.5, fontfamily="sans", fontface="bold", 
               cat.cex=2, cat.fontfamily="sans", cat.fontface="bold", cat.pos=c(315,45), cat.dist=c(0.12,0.12))


  for (gene in rownames(M8.4species.sig))
      plotM8BEB(path4, run, gene, data.4species$groups, 0.75, 0.95, paste(outpath,'M8_4taxa/M8.4species.sig.',gene, '.pdf', sep=''))

  for (gene in rownames(M8.5species.sig))
      plotM8BEB(path5, run, gene, data.5species$groups, 0.75, 0.95, paste(outpath,'M8_5taxa/M8.5species.sig.',gene, '.pdf', sep=''))  




#############################################################################################
# Positive selection on internal branches

  BS.bombus_apis4.sig <- extractBSAlt(data.4species$models$BSAlt,      data.4species$alignments, data.4species$tests$BStest)
  BS.allinternal.sig  <- extractBSAlt(data.5species$models$BSAltTree1, data.5species$alignments, data.5species$tests$BStest1)
  BS.bombus_apis5.sig <- extractBSAlt(data.5species$models$BSAltTree2, data.5species$alignments, data.5species$tests$BStest2)
  BS.megachile.sig    <- extractBSAlt(data.5species$models$BSAltTree3, data.5species$alignments, data.5species$tests$BStest3)
  BS.bombus.sig       <- extractBSAlt(data.5species$models$BSAltTree4, data.5species$alignments, data.5species$tests$BStest4)
  BS.apis.sig         <- extractBSAlt(data.5species$models$BSAltTree5, data.5species$alignments, data.5species$tests$BStest5)
  
  BS.foot  <- expression(""^a~"Group identifiers are from OrthoDB 6 (http://cegg.unige.ch/orthodb6).",
                         ""^b~"Unless otherwise specified, gene names are taken from the"~italic("A. mellifera")~"or"~italic("D. melanogaster")~"orthologs.", 
                         ""^c~"Total number of codons in the alignment after trimming with Gblocks.", 
                         ""^d~"From LRT comparing Branch-site model A to a constrained version with"~omega[2]~"=1.",
                         ""^e~"Multiple test correction by the method of Benjamini and Hochberg to control the false discovery rate.",
                         ""^f~"Sites are classified as under positive selection if the Bayesian posterior probability > 0.75 (> 0.95 in bold).",
                         "  Reference sequence taken from"~italic("A. mellifera"))

  BS.bombus_apis4.sig.table <- makeTable(getPositiveSelectionTable(path4, run, 'BS.Tree1', BS.bombus_apis4.sig, data.4species$groups, 0.75, 0.95), 
                                         caption = expression("Genes under positive selection (using FDR < 0.05) on the branch between"~italic(Bombus)~"and"~italic(Apis)~"(4 taxa tree)"), 
                                         title   = possel.title, footnotes = BS.foot,            
                                         col=3, floatcols=c(5,6), justcol=7, maxwidth=50, draw=FALSE)
  # All internal table doesn't work because there is only one group
  BS.bombus_apis5.sig.table <- makeTable(getPositiveSelectionTable(path5, run, 'BS.Tree2',BS.bombus_apis5.sig, data.5species$groups, 0.75, 0.95), 
                                         caption = expression("Genes under positive selection (using FDR < 0.05) on the branch between"~italic(Bombus)~"and"~italic(Apis)~"(5 taxa tree)"), 
                                         title   = possel.title, footnotes = BS.foot,            
                                         col=3, floatcols=c(5,6), justcol=7, maxwidth=50, draw=FALSE)
  BS.megachile.sig.table    <- makeTable(getPositiveSelectionTable(path5, run, 'BS.Tree3',BS.megachile.sig, data.5species$groups, 0.75, 0.95), 
                                         caption = expression("Genes under positive selection (using FDR < 0.05) on the branch to"~italic(Megachile)~"(5 taxa tree)"), 
                                         title   = possel.title, footnotes = BS.foot,            
                                         col=3, floatcols=c(5,6), justcol=7, maxwidth=50, draw=FALSE)
  BS.bombus.sig.table       <- makeTable(getPositiveSelectionTable(path5, run, 'BS.Tree4',BS.bombus.sig,data.5species$groups, 0.75, 0.95), 
                                         caption = expression("Genes under positive selection (using FDR < 0.05) on the branch to"~italic(Bombus)~"(5 taxa tree)"),
                                         title   = possel.title, footnotes = BS.foot,            
                                         col=3, floatcols=c(5,6), justcol=7, maxwidth=50, draw=FALSE)
  BS.apis.sig.table         <- makeTable(getPositiveSelectionTable(path5, run, 'BS.Tree5',BS.apis.sig,data.5species$groups, 0.75, 0.95), 
                                         caption = expression("Genes under positive selection (using FDR < 0.05) on the branch to"~italic(Apis)~"(5 taxa tree)"),
                                         title   = possel.title, footnotes = BS.foot,            
                                         col=3, floatcols=c(5,6), justcol=7, maxwidth=50, draw=FALSE)


  # Genes positively selected on the different branches for the 5 species case
  venn.diagram(list("Apis"=rownames(BS.apis.sig), "Bombus"=rownames(BS.bombus.sig), "Megachile"=rownames(BS.megachile.sig)), 
               filename=paste(outpath,'Venn/BS_5species_sig.tiff',sep='/'), col="black", lty=1, lwd=4, margin=0.1,
               fill=pal.trans[c(4,1,5)], alpha=1, cex=2.5, fontfamily="sans", fontface="bold", 
               cat.cex=2, cat.fontfamily="sans", cat.fontface="bold", cat.pos=c(315,45,180), cat.dist=c(0.09,0.09,0.05))

  # Genes positvely selected on the complete phylogeny and only on the branch connecting Bombus and Apis on the 4 species case
  venn.diagram(list("Complete\nphylogeny"=rownames(M8.4species.sig), "Connecting\nbranch"=rownames(BS.bombus_apis4.sig)),
               filename=paste(outpath,'Venn/PositiveSelection_4species_sig.tiff',sep='/'), col="black", lty=1, lwd=4, margin=0.15,
               fill=pal.dark[c(1,3)], alpha=0.50, cex=2.5, fontfamily="sans", fontface="bold", 
               cat.cex=2, cat.fontfamily="sans", cat.fontface="bold", cat.pos=c(330,30), cat.dist=c(0.12,0.12))               

  # Genes positvely selected on the complete phylogeny and only on the branch connecting Bombus and Apis on the 5 species case
  venn.diagram(list("Complete\nphylogeny"=rownames(M8.5species.sig), "Connecting\nbranch"=rownames(BS.bombus_apis5.sig)),
               filename=paste(outpath,'Venn/PositiveSelection_5species_sig.tiff',sep='/'), col="black", lty=1, lwd=4, margin=0.15,
               fill=pal.dark[c(1,3)], alpha=0.50, cex=2.5, fontfamily="sans", fontface="bold", 
               cat.cex=2, cat.fontfamily="sans", cat.fontface="bold", cat.pos=c(330,30), cat.dist=c(0.12,0.12))


#############################################################################################
# Differences in selection between clades  (Clade model D)


# Differences between Bombus and Apis clades
  CladeD.Bombus_Apis4.sig    <- extract2Clade(data.4species$models$CladeD, data.4species$alignments, data.4species$tests$CladetestD)
  CladeD.Bombus_Apis4.notsig <- summarizeM3(data.4species$models$M3, data.4species$alignments)[getNotSignificantNames(data.4species$tests$CladetestD),]
  
  # Bombus vs Apis plots
  plotClade(CladeD.Bombus_Apis4.sig, M8.4species.sig, BS.bombus_apis4.sig, CladeD.Bombus_Apis4.notsig, clades=c('Bombus','Apis'),
            outfile=paste(outpath,'/CladeD_Bombus_Apis/CladeD.Bombus_Apis4.sig.pdf', sep=''), split=2.5, palidx=c(2,1), size=10)
  plotClade(CladeD.Bombus_Apis4.sig, M8.4species.sig, BS.bombus_apis4.sig, CladeD.Bombus_Apis4.notsig, clades=c('Bombus','Apis'),
            outfile=paste(outpath,'/CladeD_Bombus_Apis/CladeD.Bombus_Apis4.sig.orthodblabels.pdf', sep=''), split=2.5, palidx=c(2,1), size=10, labelfun=getOrthoDB)  
  plotClade(CladeD.Bombus_Apis4.sig, M8.4species.sig, BS.bombus_apis4.sig, CladeD.Bombus_Apis4.notsig, clades=c('Bombus','Apis'),
            outfile=paste(outpath,'/CladeD_Bombus_Apis/CladeD.Bombus_Apis4.sig.labels.pdf', sep=''), split=2.5, palidx=c(2,1), size=10, labelfun=function(x){getGroupInfo(x,data.4species$groups, "name")})  
  plotClade(CladeD.Bombus_Apis4.sig, M8.4species.sig, BS.bombus_apis4.sig, clades=c('Bombus','Apis'), 
            outfile=paste(outpath,'/CladeD_Bombus_Apis/CladeD.Bombus_Apis4.sig.siglabels.pdf', sep=''), split=2.5, palidx=c(2,1), size=10, labelfun=function(x){getGroupInfo(x,data.4species$groups, "name")})
  
  # Divide by classification, no apparent pattern
  #for (path in unique(data.4species$groups[,2])) {
  #  plotClade(CladeD.Bombus_Apis4.sig, M8.4species.sig, BS.bombus_apis4.sig, CladeD.Bombus_Apis4.notsig,
  #            outfile=paste(outpath,'/CladeD_Bombus_Apis/CladeD.Bombus_Apis4.sig.',path,'.pdf', sep=''), split=2.5, pattern=path, palidx=c(2,1))  
  #}

  # Bombus vs Apis table
  perm   <- order(CladeD.Bombus_Apis4.sig[,'p2'], decreasing=TRUE)
  
  clade.title <- expression("OrthoDB Group"^a~" ", "Gene"^b~"                              ", "Classification                                                     ", 
                            "Sites"^c, italic(p)*"-value"^d, italic(q)*"-value"^e, "Proportion"^f~" ")
  clade.foot  <- expression(""^a~"Group identifiers are from OrthoDB 6 (http://cegg.unige.ch/orthodb6).",
                            ""^b~"Unless otherwise specified, gene names are taken from the"~italic("A. mellifera")~"or"~italic("D. melanogaster")~"orthologs.", 
                            ""^c~"Total number of codons in the alignment after trimming with Gblocks.", 
                            ""^d~"From LRT comparing Clade model D to M3 with 3 rate categories.",
                            ""^e~"Multiple test correction by the method of Benjamini and Hochberg to control the false discovery rate.",
                            ""^f~"Maximum likelihood estimate for the proportion of sites evolving differently between clades.")  
  
  clade.4.table <- cbind(expandLabels(rownames(CladeD.Bombus_Apis4.sig), data.4species$groups), CladeD.Bombus_Apis4.sig[,c(2,7,8,5)])

  CladeD.Bombus_Apis4.sig.table <- makeTable(clade.4.table[perm,], 
                                             caption=expression("Genes that tested significant for evolving under different selective pressures (FDR < 0.05) between"~italic("Bombus")~"and"~italic("Apis")*", according to Clade model D (4 taxa tree)"),
                                             title=clade.title, footnotes=clade.foot, col=2, floatcols=c(5,6,7), draw=FALSE)



  M8idxs <- rownames(CladeD.Bombus_Apis4.sig) %in% intersect(rownames(M8.4species.sig),     rownames(CladeD.Bombus_Apis4.sig))
  BSidxs <- rownames(CladeD.Bombus_Apis4.sig) %in% intersect(rownames(BS.bombus_apis4.sig), rownames(CladeD.Bombus_Apis4.sig))
  clade.4.table[M8idxs,1] <- paste(clade.4.table[M8idxs,1],'A')
  clade.4.table[BSidxs,1] <- paste(clade.4.table[BSidxs,1],'B')



  CladeD.Bombus_Apis4.sig.atable <- makeTable(clade.4.table[perm,], 
                                              caption=expression("Genes that tested significant for evolving under different selective pressures (FDR < 0.05) between"~italic("Bombus")~"and"~italic("Apis")*", according to Clade model D (4 taxa tree)",
                                                                 "In addition, red rows tested significant for positive selection on the whole phylogeny, while purple rows tested significant for positive selection on the branch",
                                                                 "connecting"~italic("Bombus")~"and"~italic("Apis")*"."),
                                              title=clade.title, footnotes=clade.foot, col=2, floatcols=c(5,6,7), draw=FALSE)


# Differences between Social and Non-social clades 
  CladeD.Social_Non.sig    <- extract2Clade(data.5species$models$CladeDTree1, data.5species$alignments, data.5species$tests$CladeDtest1)
  CladeD.Social_Non.notsig <- summarizeM3(data.5species$models$M3, data.5species$alignments)[getNotSignificantNames(data.5species$tests$CladeDtest1),]

  
  


  # Social vs Non-social plots
  plotClade(CladeD.Social_Non.sig, M8.5species.sig, BS.megachile.sig, CladeD.Social_Non.notsig, clades=c('asocial','social'),
          outfile=paste(outpath,'/CladeD_Social_Non/CladeD.Social_Non.sig.pdf', sep=''), split=2.5, palidx=c(2,1), size=10)  
  plotClade(CladeD.Social_Non.sig, M8.5species.sig, BS.megachile.sig, CladeD.Social_Non.notsig, clades=c('asocial','social'),
          outfile=paste(outpath,'/CladeD_Social_Non/CladeD.Social_Non.sig.orthodblabels.pdf', sep=''), split=2.5, palidx=c(2,1), size=10, labelfun=getOrthoDB)  
  plotClade(CladeD.Social_Non.sig, M8.5species.sig, BS.megachile.sig, CladeD.Social_Non.notsig, clades=c('asocial','social'),
            outfile=paste(outpath,'/CladeD_Social_Non/CladeD.Social_Non.sig.labels.pdf', sep=''), split=2.5, palidx=c(2,1), size=10, labelfun=function(x){getGroupInfo(x,data.5species$groups, "name")})  
  plotClade(CladeD.Social_Non.sig, M8.5species.sig, BS.megachile.sig, clades=c('asocial','social'),
            outfile=paste(outpath,'/CladeD_Social_Non/CladeD.Social_Non.sig.siglabels.pdf', sep=''), split=2.5, palidx=c(2,1), size=10, labelfun=function(x){getGroupInfo(x,data.5species$groups, "name")})  
  # Divide by classification, no apparent pattern
  #for (path in unique(data.5species$groups[,2])) {  
  #  plotClade(CladeD.Social_Non.sig, M8.5species.sig, BS.megachile.sig, CladeD.Social_Non.notsig, clades=c('asocial','social'),
  #            outfile=paste(outpath,'/CladeD_Social_Non/CladeD.Social_Non.sig.', path, '.pdf', sep=''), split=2.5, pattern=path, palidx=c(2,1), size=10)  
  #}

  # Social vs Non-social table
  perm             <- order(CladeD.Social_Non.sig[,'p2'], decreasing=TRUE)
  
  clade.5.table <- cbind(expandLabels(rownames(CladeD.Social_Non.sig), data.5species$groups), CladeD.Social_Non.sig[,c(2,7,8,5)])
  
  CladeD.Social_Non.sig.table <- makeTable(clade.5.table[perm,], 
                                           caption="Genes that tested significant for evolving under different selective pressures (FDR < 0.05) between the social and solitary clades, according to Clade model D (5 taxa tree)",
                                           title=clade.title, footnotes=clade.foot, col=2, floatcols=c(5,6,7), draw=FALSE)


  M8idxs           <- rownames(CladeD.Social_Non.sig) %in% intersect(rownames(M8.5species.sig),  rownames(CladeD.Bombus_Apis4.sig))
  BSidxs.megachile <- rownames(CladeD.Social_Non.sig) %in% intersect(rownames(BS.megachile.sig), rownames(CladeD.Social_Non.sig))
  BSidxs.bombus    <- rownames(CladeD.Social_Non.sig) %in% intersect(rownames(BS.bombus.sig),    rownames(CladeD.Social_Non.sig))
  BSidxs.apis      <- rownames(CladeD.Social_Non.sig) %in% intersect(rownames(BS.apis.sig),      rownames(CladeD.Social_Non.sig))
  clade.5.table[M8idxs,1]           <- paste(clade.5.table[M8idxs,1],'a')
  clade.5.table[BSidxs.megachile,1] <- paste(clade.5.table[BSidxs.megachile,1],'M')
  clade.5.table[BSidxs.bombus,1]    <- paste(clade.5.table[BSidxs.bombus,1],'B')
  clade.5.table[BSidxs.apis,1]      <- paste(clade.5.table[BSidxs.apis,1],'A')

  CladeD.Social_Non.sig.atable <- makeTable(clade.5.table[perm,], 
                                           caption=expression("Genes that tested significant for evolving under different selective pressures (FDR < 0.05) between the social and solitary clades, according to Clade model D (5 taxa tree)",
                                                   "In addition, red rows tested significant for positive selection on the whole phylogeny, while purple rows tested significant for positive selection on the branch to"~italic("Megachile"),
                                                   "orange for the branch to"~italic("Apis")~"and blue for the branch to"~italic("Bombus")*"."),
                                           title=clade.title, footnotes=clade.foot, col=2, floatcols=c(5,6,7), draw=FALSE)

# Difference between all 3 clades

  # Genes that are significant for 3 different clades on 5 species case (forward selection)
  CladeD.3clade.sig <- extract2Clade(data.5species$models$CladeDTree2, data.5species$alignments, data.5species$tests$CladeDtest3)[intersect(getSignificantNames(data.5species$tests$CladeDtest3),
                                                                                                                                            getSignificantNames(data.5species$tests$CladeDtest1)),]

  perm              <- order(CladeD.3clade.sig[,'p2'], decreasing=TRUE)

  CladeD.3clade.sig.table <- makeTable(cbind(expandLabels(as.matrix(rownames(CladeD.3clade.sig)), data.5species$groups)[perm,], CladeD.3clade.sig[perm,c(2,7,8,5)]), 
                                           caption=expression("Genes that in addition to testing significant for evolving under different evolutionary pressures between the social and solitary clades, also tested ",
                                                              "significant for evolving under different pressures between all 3 clades ("*italic("Bombus")*","~italic("Apis")*","~italic("Megachile")*") using Clade model D (FDR < 0.05, 5 taxa tree)"),
                                           title=clade.title, footnotes=clade.foot, col=2, floatcols=c(5,6,7), draw=FALSE)


  
#############################################################################################
# Differences in selection between clades  (Clade model C)
#
# Compare Clade model C to Clade model D (Don't do complete analysis for Clade model C as well)


  CladeC.Bombus_Apis4.sig  <- extract2Clade(data.4species$models$CladeC, data.4species$alignments, data.4species$tests$CladetestC)
  CladeC.Social_Non.sig    <- extract2Clade(data.5species$models$CladeCTree1, data.5species$alignments, data.5species$tests$CladeCtest1)
  CladeC.3clade.sig <- extract2Clade(data.5species$models$CladeCTree2, data.5species$alignments, data.5species$tests$CladeCtest3)[intersect(getSignificantNames(data.5species$tests$CladeCtest3),
                                                                                                                                            getSignificantNames(data.5species$tests$CladeCtest1)),]
  
  
  


# Overlap between Clade model C and model D
  venn.diagram(list("Clade C"=rownames(CladeC.Bombus_Apis4.sig), "Clade D"=rownames(CladeD.Bombus_Apis4.sig)),
               filename=paste(outpath,'Venn/CladeDiff_Bombus_Apis4.tiff',sep='/'))

  venn.diagram(list("Clade C"=rownames(CladeC.Social_Non.sig), "Clade D"=rownames(CladeD.Social_Non.sig)),  
               filename=paste(outpath,'Venn/CladeDiff_Social_Non.tiff',sep='/'))
               
  venn.diagram(list("Clade C"=CladeC.3clade.sig, "Clade D"=CladeD.3clade.sig),
               filename=paste(outpath,'Venn/CladeDiff_3Clade.tiff',sep='/'))               




#############################################################################################
# Summary of everything

#plotGradient(data.4species$models$M0, data.4species$groups, 100, M8.4species.sig, BS.bombus_apis4.sig, BS.apis.sig, BS.bombus.sig, BS.megachile.sig, CladeD.Bombus_Apis4.sig, CladeD.Social_Non.sig)


  names <- as.character(data.4species$alignments[,1])
  M8.4species.col <- M8.5species.col <- BS.bombus_apis4.col <- BS.bombus_apis5.col <- BS.bombus.col <- BS.apis.col <- BS.megachile.col <- CladeD.bombus_apis.col <- CladeD.social_non.col <- CladeD.3clade.col <- 
  matrix(nrow=length(names), ncol=1, " ")
  M0.4species <- cbind(data.4species$alignments[,2], data.4species$models$M0$w)
  M0.5species <- matrix(nrow=length(names), ncol=2, " ")
  notin5 <- setdiff(names, data.5species$alignments[,1])
  for (i in 1:length(names)) {
      if (length(grep(names[i], rownames(M8.4species.sig))) > 0) M8.4species.col[i]         <- "   X"
      if (length(grep(names[i], rownames(M8.5species.sig))) > 0) M8.5species.col[i]         <- "   X"
      if (length(grep(names[i], rownames(BS.bombus_apis4.sig))) > 0) BS.bombus_apis4.col[i] <- "   X"
      if (length(grep(names[i], rownames(BS.bombus_apis5.sig))) > 0) BS.bombus_apis5.col[i] <- "   X"
      if (length(grep(names[i], rownames(BS.bombus.sig))) > 0)           BS.bombus.col[i]           <- "   X"
      if (length(grep(names[i], rownames(BS.apis.sig))) > 0)             BS.apis.col[i]             <- "   X"
      if (length(grep(names[i], rownames(BS.megachile.sig))) > 0)        BS.megachile.col[i]        <- "   X"
      if (length(grep(names[i], rownames(CladeD.Bombus_Apis4.sig))) > 0) CladeD.bombus_apis.col[i]  <- "   X"
      if (length(grep(names[i], rownames(CladeD.Social_Non.sig))) > 0)   CladeD.social_non.col[i]   <- "   X"
      if (length(grep(names[i], rownames(CladeD.3clade.sig))) > 0)       CladeD.3clade.col[i]       <- "   X"      
      
      if (names[i] %in% notin5) 
          M0.5species[i,] <- c("N/A", "N/A")
      else 
          M0.5species[i,] <- c(data.5species$alignments[which(data.5species$alignments[,1] == names[i]),2], data.5species$models$M0$w[which(data.5species$models$M0[,1] == names[i])])
  }

  

  Summary.4species <- cbind(expandLabels(as.matrix(names), data.4species$groups), M0.4species, M8.4species.col, BS.bombus_apis4.col, CladeD.bombus_apis.col)
  Summary.5species <- cbind(expandLabels(as.matrix(names), data.4species$groups), M0.5species, M8.5species.col, BS.bombus_apis5.col, BS.bombus.col, BS.apis.col, BS.megachile.col, CladeD.social_non.col)

  perm             <- order(Summary.4species[,3], Summary.4species[,2])


  Summary.4species.table <- makeTable(Summary.4species[perm,], caption="Summary of models on the 4 taxa tree.", col=1, draw=FALSE, 
                                      title=expression("OrthoDB Group"^a~" ", "Gene"^b~"                              ", "Classification                                                     ", 
                                                       "Sites"^c, "Global"~omega^1, "M8 vs M7"^2, "Branch-site"^3, "Clade D"^4), 
                                      footnotes=expression(""^a~"Group identifiers are from OrthoDB 6 (http://cegg.unige.ch/orthodb6).",
                                                           ""^b~"Unless otherwise specified, gene names are taken from the"~italic("A. mellifera")~"or"~italic("D. melanogaster")~"orthologs.", 
                                                           ""^c~"Total number of codons in the alignment after trimming with Gblocks.", 
                                                           " ",
                                                           ""^1~"Across the whole phylogeny using M0 model.",
                                                           ""^2~"Positive selection across the whole phylogeny.",
                                                           ""^3~"Positive selection on the branch between"~italic("Bombus")~"and"~italic("Apis")*".",
                                                           ""^4~"Different selective pressures between"~italic("Bombus")~"and"~italic("Apis")*".",
                                                           " ", 
                                                           " Using FDR < 0.05 on all tests."))


  Summary.5species.table <- makeTable(Summary.5species[perm,], caption="Summary of models on the 5 taxa tree.", col=4, draw=FALSE, 
                                      title=expression("OrthoDB Group"^a~" ", "Gene"^b~"                              ", "Classification                                                     ", 
                                                       "Sites"^c, "Global"~omega^1, "M8 vs M7"^2, "Branch-site"^3, "Branch-site"^4, "Branch-site"^5, "Branch-site"^6, "Clade D"^7), 
                                    footnotes=expression(""^a~"Group identifiers are from OrthoDB 6 (http://cegg.unige.ch/orthodb6).",
                                                         ""^b~"Unless otherwise specified, gene names are taken from the"~italic("A. mellifera")~"or"~italic("D. melanogaster")~"orthologs.", 
                                                         ""^c~"Total number of codons in the alignment after trimming with Gblocks.", 
                                                         " ",
                                                         ""^1~"Across the whole phylogeny using M0 model.",
                                                         ""^2~"Positive selection across the whole phylogeny.",
                                                         ""^3~"Positive selection on the branch between"~italic("Bombus")~"and"~italic("Apis")*".",
                                                         ""^4~"Positive selection on the branch to"~italic("Bombus")*".",
                                                         ""^5~"Positive selection on the branch to"~italic("Apis")*".",
                                                         ""^6~"Positive selection on the branch to"~italic("Megachile")*".",
                                                         ""^7~"Different selective pressures between"~italic("social")~"and"~italic("solitary")*".",
                                                         " ", 
                                                         " Using FDR < 0.05 on all tests."))


#############################################################################################
# Save Tables
dev.off()

# 4 taxa 
pdf(paste(outpath,"Tables_4taxa/M0.4species.table.pdf", sep=''), height=4, width=8)
grid.draw(M0.4species.table)
dev.off()
    
pdf(paste(outpath,"Tables_4taxa/M0.fastest.4species.table.pdf", sep=''), height=11, width=10)
grid.draw(M0.fastest.4species.table)
dev.off()

pdf(paste(outpath,"Tables_4taxa/M0.slowest.4species.table.pdf", sep=''), height=11, width=10)
grid.draw(M0.slowest.4species.table)
dev.off()
            
pdf(paste(outpath,"Tables_4taxa/M8.4species.sig.table.pdf", sep=''), height=18, width=20) 
grid.draw(M8.4species.sig.table)
dev.off()
    
pdf(paste(outpath,"Tables_4taxa/BS.Bombus_Apis4.sig.table.pdf", sep=''), height=15, width=20) 
grid.draw(BS.bombus_apis4.sig.table)
dev.off()

pdf(paste(outpath,"Tables_4taxa/CladeD.Bombus_Apis4.sig.table.pdf", sep=''), height=30, width=15) 
grid.draw(CladeD.Bombus_Apis4.sig.table)
dev.off()

pdf(paste(outpath,"Tables_4taxa/CladeD.Bombus_Apis4.sig.table.markup.pdf", sep=''), height=30, width=15) 
grid.draw(CladeD.Bombus_Apis4.sig.atable)
dev.off()

pdf(paste(outpath,"Tables_4taxa/Summary.4species.table.pdf", sep=''), height=50, width=15) 
grid.draw(Summary.4species.table)
dev.off()

# 5 taxa 
pdf(paste(outpath,"Tables_5taxa/M0.5species.table.pdf", sep=''), height=4, width=8)
grid.draw(M0.5species.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/M0.fastest.5species.table.pdf", sep=''), height=11, width=10)
grid.draw(M0.fastest.5species.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/M0.slowest.5species.table.pdf", sep=''), height=11, width=10)
grid.draw(M0.slowest.5species.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/M8.5species.sig.table.pdf", sep=''), height=8, width=20) 
grid.draw(M8.5species.sig.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/BS.Bombus_Apis5.sig.table.pdf", sep=''), height=10, width=20) 
grid.draw(BS.bombus_apis5.sig.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/BS.Bombus.sig.table.pdf", sep=''), height=10, width=20) 
grid.draw(BS.bombus.sig.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/BS.Apis.sig.table.pdf", sep=''), height=10, width=20) 
grid.draw(BS.apis.sig.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/BS.Megachile.sig.table.pdf", sep=''), height=30, width=20) 
grid.draw(BS.megachile.sig.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/CladeD.Social_Non.sig.table.pdf", sep=''), height=30, width=15) 
grid.draw(CladeD.Social_Non.sig.table)
dev.off()

pdf(paste(outpath,"Tables_5taxa/CladeD.Social_Non.sig.table.markup.pdf", sep=''), height=30, width=15) 
grid.draw(CladeD.Social_Non.sig.atable)
dev.off()

#pdf(paste(outpath,"Tables_5taxa/CladeD.3clade.sig.table.pdf", sep=''), height=10, width=15) 
#grid.draw(CladeD.3clade.sig.table)
#dev.off()

pdf(paste(outpath,"Tables_5taxa/Summary.5species.table.pdf", sep=''), height=50, width=15) 
grid.draw(Summary.5species.table)
dev.off()
