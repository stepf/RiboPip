#!/usr/bin/env Rscript

# Wrapping around DESeq1 and DESeq2 for ribopip pipeline orchestration
# Expects tab-delimited input file containing read counts for each experimental
# replicate. For example:

# Gene_name	A_replicate1	A_replicate2	A_replicate3	B_replicate1	B_replicate2
# GeneX	12	28	5	80	90
# GeneY	5	80	100	2	3
# GeneZ	1390	700	480	13	9

args <- commandArgs(trailingOnly=TRUE)
exp_args <- 9

if (is.na(args[exp_args])) {
  stop(
    "Usage: $0 <deseq_version> <input_counts_table> <min_count_threshold> " +
    "<experiment1_name> <experiment2_name> " +
    "<experiment1_number> <experiment2_number> " +
    "<output_all_results> <output_padj_filtered_results>\n" +
    "e.g. $0 1 counts.tsv 128 treatment control results.all.csv results.pad.csv"
  )
}

deseq_ver <- args[1]
infname <- args[2]
minCount <- args[3]
treatment <- args[4]
control <- args[5]
treatment_num <- args[6]
control_num <- args[7]
outallname <- args[8]
outpadjname <- args[9]

cat("counts table: ",infname,"\n")
countTable <- read.table(infname, header=TRUE, row.names=1 )
cat("features:",nrow(countTable),"\n")

countTable <- countTable[rowSums(countTable) >= minCount,]
cat("features with sum(counts) >=",minCount,":",nrow(countTable),"\n")

if(deseq_ver == 1) {
  suppressMessages(require("DESeq"))
  total_num <- as.integer(treatment_num) + as.integer(control_num)
  design <- data.frame(
    row.names <- colnames(countTable),
    condition <- c(
      rep('treatment', treatment_num),
      rep('control', control_num)
    ),
    libType <- c(rep("single-end", total_num))
  )

  condition <- design$condition
  cds <- newCountDataSet(countTable, condition)
  cds <- estimateSizeFactors(cds)
  cat("Size factors:\n")
  sizeFactors(cds)

  cds <- estimateDispersions(
    cds,
    method="blind",
    sharingMode="fit-only",
    fitType="local"
  )
  res <- nbinomTest(cds, 'treatment', 'control')
} else if(deseq_ver == 2) {
  suppressMessages(require("DESeq2"))
  colData <- data.frame(
    row.names = colnames(countTable),
    condition = c(rep('treatment', treatment_num), rep('control', control_num)),
    type = c(rep("single-end", treatment_num + control_num))
  )

  colData$condition <- relevel(colData$condition, 'control')

  dataset <- DESeqDataSetFromMatrix(
    countData=countTable,
    colData=colData,
    design=~condition
  )
  analysis <- DESeq(dataset)
  res <- results(analysis)
}

res_filtered <- res[res$padj < 0.05, ]
write.table(res[order(res$log2FoldChange),], file=outallname, sep="\t")
write.table(res_filtered[order(res$log2FoldChange),], file=outpadjname, sep="\t")
cat("significant changes:",nrow(res),"\n")
cat(warnings())
