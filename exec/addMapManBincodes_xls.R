require(GeneFamilies)
library(gdata)

message("USAGE: Rscript path/2/GeneFamilies/exec/addMapManBincodes.R path/2/DGE_Results/xls path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

# Declare funtions to use

#' Function to add MapMan-Bins where possible:
addMapManBins <- function(expr.tbl) {
    cbind(expr.tbl, Reduce(rbind, mclapply(expr.tbl$Name, function(gene) {
        gene.san <- toLowerCutTail(gene)
        if (gene.san %in% maize.mapMan$IDENTIFIER.san) {
            maize.mapMan.i <- maize.mapMan[which(maize.mapMan$IDENTIFIER.san ==
                gene.san), ]
            data.frame(MapMan.BINCODE = paste(sort(unique(maize.mapMan.i$BINCODE)),
                collapse = ";"), MapMan.NAME = paste(sort(unique(maize.mapMan.i$NAME)),
                collapse = ";"), MapMan.Gene.IDs = paste(sort(unique(maize.mapMan.i$IDENTIFIER)),
                collapse = ";"), stringsAsFactors = FALSE)
        } else {
            data.frame(MapMan.BINCODE = NA, MapMan.NAME = NA, MapMan.Gene.IDs = NA,
                stringsAsFactors = FALSE)
        }
    })))
}

#' Function to execute Fischer test:
maizeFischerTest <- function(genes.de, genes.not.de, univ.genes = union(genes.de,
    genes.not.de)) {
    res.df <- Reduce(rbind, mclapply(mmBins, function(mmBin) {
        mmBin.expr.genes <- intersect(maize.mapMan[which(maize.mapMan$BINCODE ==
            mmBin), "IDENTIFIER.san"], maize.expr.w.mapMan)
        if (length(mmBin.expr.genes) > 0) {
            cont.tbl <- generateContingencyTable(genes.de, genes.not.de,
                mmBin.expr.genes, setdiff(univ.genes, mmBin.expr.genes),
                "DE", "MapManBin")
            p.val <- fisher.test(cont.tbl, alternative = "greater")$p.value
            data.frame(BINCODE = mmBin, p.value = p.val, stringsAsFactors = FALSE)
        } else NULL
    }))
    if (!is.null(res.df)) {
        # Correct for multiple hypothesis testing:
        res.df$p.adjusted <- p.adjust(res.df$p.value, method = "fdr")
        # Retain significant ones only:
        res.df.sign <- res.df[which(res.df$p.adjusted <= 0.05), ]
        if (nrow(res.df.sign) > 0) {
            # Add MapMan-Bin Names and Descriptions:
            res.df.sign$NAME <- sapply(res.df.sign$BINCODE, function(mmBin) maize.mapMan[which(maize.mapMan$BINCODE ==
                mmBin), "NAME"][[1]])
            res.df.sign$DESCRIPTION <- sapply(res.df.sign$BINCODE, function(mmBin) maize.mapMan[which(maize.mapMan$BINCODE ==
                mmBin), "DESCRIPTION"][[1]])
            # Done:
            res.df.sign
        } else {
            warning("No significant P-Values found.")
            res.df
        }
    } else NULL
}

# Here we begin the analysis

maize.mapMan$IDENTIFIER.san <- toLowerCutTail(maize.mapMan$IDENTIFIER)
maize.w.mapMan <- intersect(toLowerCutTail(names(maize.aas)), toLowerCutTail(maize.mapMan$IDENTIFIER))
maize.genes.expr <- c()
baseDGE <- c()

#' Read list with corrected names
Name_fixed <- read.table("inst/list_corrected_names.txt", header=T)
Name_fixed <- unlist(Name_fixed)

#` Path to folder that holds the DGE results
dir <- input.args[[1]]

#' Create list of all downregulated or upregulated files in folder
file_list <- list.files(path=dir, pattern=".xls") 
basename <- gsub(".xls", "", file_list, perl = TRUE)


for (i in 1:length(file_list)){

  #' Read all expression tables to analyze:
  assign(basename[i], 
  read.xls(paste(dir, file_list[i], sep=''))
  ) 
 
  # Substitute Name column with the correpondign Ensemble nomenclature.
  renameCommand <- paste(basename[i], "$Name <- Name_fixed", sep = "")
  eval(parse(text = renameCommand ))

  #` Sanitize Gene Identifiers
	# B73_BIODYN_Root_vs_B73_CONMIN_Root_down$IDENTIFIER.san <- toLowerCutTail(B73_BIODYN_Root_vs_B73_CONMIN_Root_down$X)

  sanitizeCommand <- paste(basename[i], "$IDENTIFIER.san <- toLowerCutTail(", basename[i], "$Name)", sep = "")
  eval(parse(text = sanitizeCommand ))


  #` Identify differentially expressed genes and with a log-fold change greater or equal than 2
	#   a <- B73_BIODYN_Leaf_12_vs_B73_BIODYN_Leaf_90[ which(B73_BIODYN_Leaf_12_vs_B73_BIODYN_Leaf_90$FDR.p.value <= 0.05 & B73_BIODYN_Leaf_12_vs_B73_BIODYN_Leaf_90$Log..fold.change > 1), ]
	# b <- B73_BIODYN_Leaf_12_vs_B73_BIODYN_Leaf_90[ which(B73_BIODYN_Leaf_12_vs_B73_BIODYN_Leaf_90$FDR.p.value <= 0.05 & B73_BIODYN_Leaf_12_vs_B73_BIODYN_Leaf_90$Log..fold.change < 1), ]
	
	DgeDownCommand <- paste("a <- ", basename[i], "[ which(", basename[i], "$FDR.p.value <= 0.05 & ", basename[i],"$Log..fold.change < -1), ]", sep = "")
	DgeUpCommand <- paste("b <- ", basename[i], "[ which(", basename[i], "$FDR.p.value <= 0.05 & ", basename[i],"$Log..fold.change > 1), ]", sep = "")
	eval(parse(text = DgeDownCommand ))
	eval(parse(text = DgeUpCommand ))
	
	DgeMergeCommand <- paste(basename[i], "_DE <- rbind(a,b)", sep = "")
	eval(parse(text = DgeMergeCommand ))




  #' Add MapMan Bincodes
	# B73_BIODYN_Root_vs_B73_CONMIN_Root_down <- addMapManBins(B73_BIODYN_Root_vs_B73_CONMIN_Root_down)
  addMapManCommand <- paste(basename[i], "_DE <- addMapManBins(", basename[i], "_DE)", sep = "")
  eval(parse(text = addMapManCommand ))

  #' Vector for saving DGE
  baseDGE[i] <- paste(basename[i], "_DE", sep = "" )

  #' Prepare exact Fischer tests:   
  extractGenesCommand <- paste( "a <- as.character(", basename[i], "_DE$Name)", sep = "")
  eval(parse(text = extractGenesCommand ))
  maize.genes.expr <- unique(append(maize.genes.expr, a))

  #' Write tables with MapMan Info:
	#` write.table(B73_BIODYN_Root_vs_B73_CONMIN_Root_down, file.path(input.args[[2]], "ReconstructXls", "B73_BIODYN_Root_vs_B73_CONMIN_Root_down_wMapMan.txt"), row.names = FALSE, sep = "\t")
  writeMapManCommand <- paste('write.table(', baseDGE[i], ', file.path(input.args[[2]], "ReconstructXls", "', baseDGE[i], '_wMapMan.txt"), row.names = FALSE, sep = "\t")', sep = '')
  eval(parse(text = writeMapManCommand ))

  
}

save(list = baseDGE, file = file.path(input.args[[2]], "ReconstructXls", "DGE_LeafwMapManBins.RData"))

maize.expr.w.mapMan <- intersect(maize.w.mapMan, toLowerCutTail(maize.genes.expr))
mmBins <- unique(maize.mapMan$BINCODE)

#' Fisher tests:
  #` B73_BIODYN_Root_vs_B73_CONMIN_Root_down_fish <- maizeFischerTest(B73_BIODYN_Root_vs_B73_CONMIN_Root_down$IDENTIFIER.san, setdiff(maize.expr.w.mapMan, B73_BIODYN_Root_vs_B73_CONMIN_Root_down$IDENTIFIER.san))

message("Running Fisher test")

for (i in 1:length(basename)){
  FisherCommand <- paste(basename[i], "_fish <- maizeFischerTest(unique(", basename[i], "_DE$IDENTIFIER.san), setdiff(maize.expr.w.mapMan, unique(", basename[i], "_DE$IDENTIFIER.san)))", sep = "")
  eval(parse(text = FisherCommand ))

	#' Write tables with Fisher results:
	  # write.table(B73_BIODYN_Root_vs_B73_CONMIN_Root_down_fish, file.path(input.args[[2]], "ReconstructXls", "B73_BIODYN_Root_vs_B73_CONMIN_Root_down_Fisher.txt"), row.names = FALSE, sep = "")

	getlength <- paste('len <- ncol(', basename[i], '_fish)', sep = '')
	eval(parse(text = getlength))

	if (len > 3 ) {
		writeFisherCommand <- paste('write.table(', basename[i], '_fish, file.path(input.args[[2]], "ReconstructXls", "', basename[i], '_Fisher.txt"), row.names = FALSE, sep = "\t")', sep = '')
		eval(parse(text = writeFisherCommand ))
	} else NULL
}


#' Save results:

message("Saving results")
fisher_results <- ls(pattern = "*_fish")
save(list = fisher_results, file = file.path(input.args[[2]], "ReconstructXls", "DGE_LeafFisher.RData"))
