library("DESeq2")

counts = read.csv(file=snakemake@input$counts, sep='\t')
annots = read.csv(file=snakemake@input$annotations, sep='\t')

organism = snakemake@wildcards$organism
experiment = snakemake@wildcards$experiment
comparison = snakemake@wildcards$comparison

exp_design = as.formula(snakemake@config[["Comparisons"]][[organism]][[experiment]]$Design)

sample_use_df = read.csv(file=snakemake@config$input[['sample table']], sep='\t')

if ('Included Samples' %in% names(snakemake@config[["Comparisons"]][[organism]][[experiment]])) {
    included_samples = snakemake@config[["Comparisons"]][[organism]][[experiment]][['Included Samples']]
    sample_use_df = subset(sample_use_df, sample %in% included_samples)
} else if ('Excluded Samples' %in% names(snakemake@config[["Comparisons"]][[organism]][[experiment]])) {
    excluded_samples = snakemake@config[["Comparisons"]][[organism]][[experiment]][['Excluded Samples']]
    sample_use_df = subset(sample_use_df, !sample %in% excluded_samples)
}

row.names(sample_use_df) <- sample_use_df$sample

counts_subset = subset(counts, organism==snakemake@wildcards$organism)
rownames(counts_subset) = counts_subset$ID
counts_subset = subset(counts_subset, select=rownames(sample_use_df))

annots_subset = subset(annots, organism==snakemake@wildcards$organism)
rownames(annots_subset) = annots_subset$ID

# organism must be present in all samples. Good test to tell if it is from experience...
if (min(colMeans(counts_subset[sapply(counts_subset, is.numeric)],na.rm=TRUE)) < 1){
    stop('organism must be present in all samples')
}

dds = DESeqDataSetFromMatrix(countData = counts_subset, colData = sample_use_df, design = exp_design)
dds_processed = DESeq(dds)

vst_counts = assay(varianceStabilizingTransformation(dds_processed))
vst_counts = cbind(ID=rownames(vst_counts), vst_counts)
write.table(vst_counts, file=snakemake@output$vst, quote=FALSE, sep='\t', row.names=FALSE)
print("vst file exists:")
print(snakemake@output$vst)
print(file.exists(snakemake@output$vst))

rlog_counts = assay(rlog(dds_processed))
rlog_counts = cbind(ID=rownames(rlog_counts), rlog_counts)
write.table(rlog_counts, file=snakemake@output$rlog, quote=FALSE, sep='\t', row.names=FALSE)
print("rlog file exists:")
print(snakemake@output$rlog)
print(file.exists(snakemake@output$rlog))

comparison_factor = snakemake@config[["Comparisons"]][[organism]][[experiment]][['Comparisons']][[comparison]][['Factor']]
comparison_treatment = snakemake@config[["Comparisons"]][[organism]][[experiment]][['Comparisons']][[comparison]][['Treatment']]
comparison_control = snakemake@config[["Comparisons"]][[organism]][[experiment]][['Comparisons']][[comparison]][['Control']]
comparison_results = as.data.frame(results(dds_processed, contrast=c(comparison_factor, comparison_treatment, comparison_control)))
comparison_results = cbind(ID=rownames(comparison_results), symlog10baseMean = with(comparison_results, log10(baseMean + 1)), comparison_results)
write.table(comparison_results, file=snakemake@output$results, quote=FALSE, sep='\t', row.names=FALSE)
print("results file exists:")
print(snakemake@output$results)
print(file.exists(snakemake@output$results))

results_w_annot = merge(comparison_results, annots_subset, by=0)
write.table(results_w_annot, file=snakemake@output$results_w_annot, quote=FALSE, sep='\t', row.names=FALSE)

print("results_w_annot file exists:")
print(snakemake@output$results_w_annot)
print(file.exists(snakemake@output$results_w_annot))
