require(GeneFamilies)

message("USAGE: Rscript path/2/GeneFamilies/exec/addMapManBincodes.R path/2/DGE_Results path/2/GeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)

maize.mapMan$IDENTIFIER.san <- toLowerCutTail(maize.mapMan$IDENTIFIER)

#' Read all expression tables to analyze:

#` Path to folder that holds the DGE results
			# c("/mnt/data/yaxal/RECONSTRUCT/Analyses/SalmonQuantsAllSamples/DGE_Root/")
dir <- input.args[[1]]

#' Create list of all downregulated or upregulated files in folder
file_list <- list.files(path=dir, pattern="regulated.txt") 
basename <- gsub("regulated.txt", "", file_list, perl = TRUE)

message("Reading tables")

for (i in 1:length(file_list)){
  assign(basename[i], 
  read.table(paste(dir, file_list[i], sep=''), sep = "\t", stringsAsFactors = FALSE, header = TRUE)
)}

#` Sanitize Gene Identifiers
# Generates the following command for each comparison:
# B73_BIODYN_Root_vs_B73_CONMIN_Root_down$IDENTIFIER.san <- toLowerCutTail(B73_BIODYN_Root_vs_B73_CONMIN_Root_down$X)

message("Sanitize identifiers")

for (i in 1:length(basename)){
  sanitizeCommand <- paste(basename[i], "$IDENTIFIER.san <- toLowerCutTail(", basename[i], "$X)", sep = "")
  eval(parse(text = sanitizeCommand ))
}

#' Function to add MapMan-Bins where possible:
addMapManBins <- function(expr.tbl) {
    cbind(expr.tbl, Reduce(rbind, mclapply(expr.tbl$X, function(gene) {
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


#' Add MapMan Bincodes

# Generates the following command for each comparison:
# B73_BIODYN_Root_vs_B73_CONMIN_Root_down <- addMapManBins(B73_BIODYN_Root_vs_B73_CONMIN_Root_down)

message("Add MapMan Bincodes")

for (i in 1:length(basename)){
  addMapManCommand <- paste(basename[i], " <- addMapManBins(", basename[i], ")", sep = "")
  eval(parse(text = addMapManCommand ))
}


#' Prepare exact Fischer tests:   
maize.w.mapMan <- intersect(toLowerCutTail(names(maize.aas)), toLowerCutTail(maize.mapMan$IDENTIFIER))

maize.genes.expr <- c()
for (i in 1:length(basename)){
  extractGenesCommand <- paste( "a <- ", basename[i], "$X", sep = "")
  eval(parse(text = extractGenesCommand ))
  maize.genes.expr <- unique(append(maize.genes.expr, a))
}

maize.expr.w.mapMan <- intersect(maize.w.mapMan, toLowerCutTail(maize.genes.expr))
mmBins <- unique(maize.mapMan$BINCODE)

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

#' Fisher tests:
#` B73_BIODYN_Root_vs_B73_CONMIN_Root_down_fish <- maizeFischerTest(B73_BIODYN_Root_vs_B73_CONMIN_Root_down$IDENTIFIER.san, setdiff(maize.expr.w.mapMan, B73_BIODYN_Root_vs_B73_CONMIN_Root_down$IDENTIFIER.san))

message("Running Fisher test")

for (i in 1:length(basename)){
  FisherCommand <- paste(basename[i], "_fish <- maizeFischerTest(", basename[i], "$IDENTIFIER.san, setdiff(maize.expr.w.mapMan, ", basename[i], "$IDENTIFIER.san))", sep = "")
  eval(parse(text = FisherCommand ))
}

#' Write tables with MapMan Info:
#` write.table(B73_BIODYN_Root_vs_B73_CONMIN_Root_down, file.path(input.args[[2]], "ReconstructResults", "B73_BIODYN_Root_vs_B73_CONMIN_Root_down_wMapMan.txt"), row.names = FALSE, sep = "\t")

message("Writing files")

for (i in 1:length(basename)){
  writeMapManCommand <- paste('write.table(', basename[i], ', file.path(input.args[[2]], "ReconstructResults", "', basename[i], '_wMapMan.txt"), row.names = FALSE, sep = "\t")', sep = '')
  eval(parse(text = writeMapManCommand ))
}

#' Write tables with Fisher results:
# write.table(B73_BIODYN_Root_vs_B73_CONMIN_Root_down_fish, file.path(input.args[[2]], "ReconstructResults", "B73_BIODYN_Root_vs_B73_CONMIN_Root_down_Fisher.txt"), row.names = FALSE, sep = "")

for (i in 1:length(basename)){

	getlength <- paste('len <- ncol(', basename[i], '_fish)', sep = '')
	eval(parse(text = getlength))

	if (len > 3 ) {
		writeFisherCommand <- paste('write.table(', basename[i], '_fish, file.path(input.args[[2]], "ReconstructResults", "', basename[i], '_Fisher.txt"), row.names = FALSE, sep = "\t")', sep = '')
		eval(parse(text = writeFisherCommand ))
	} else NULL
}

#' Save results:

message("Saving results")
fisher_results <- ls(pattern = "*_fish")
save(list = basename, file = file.path(input.args[[2]], "ReconstructResults", "DGE__RootwMapManBins.RData")
save(list = fisher_results, file = file.path(input.args[[2]], "ReconstructResults", "DGE_RootFisher.RData")
