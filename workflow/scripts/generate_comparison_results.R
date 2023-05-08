library("DESeq2")

counts = read.csv(file=snakemake@input$counts, sep='\t')
annots = read.csv(file=snakemake@input$annotations, sep='\t')

organism = snakemake@wildcards$organism
experiment = snakemake@wildcards$experiment
comparison = snakemake@wildcards$comparison

saveRDS(snakemake, file = glue("{organism}.{experiment}.{comparison}.rds"))

design_string = snakemake@config$comparisons[[organism]][[experiment]]$Design


sample_use_df = read.csv(file=snakemake@config$input['sample table'], sep='\t')


# TODO: Filter for only included/excluded samples

counts_subset = subset(counts, organism==snakemake@wildcards$organism)
rownames(counts_subset) = counts_subset$ID
counts_subset = subset(counts_subset, select=rownames(sample_use_df))

annots_subset = subset(annots, organism==snakemake@wildcards$organism)
rownames(annots_subset) = annots_subset$ID

# organism must be present in all samples. Good test to tell if it is from experience...
if (min(colMeans(counts_subset[sapply(counts_subset, is.numeric)],na.rm=TRUE)) < 1){
    stop('organism must be present in all samples')
}

dds = DESeqDataSetFromMatrix(countData = counts_subset, colData = sample_use_df, design = ~ condition)
dds_processed = DESeq(dds)

vst_counts = assay(varianceStabilizingTransformation(dds_processed))
vst_counts = cbind(ID=rownames(vst_counts), vst_counts)
write.table(vst_counts, file=snakemake@output$vst, quote=FALSE, sep='\t', row.names=FALSE)

rlog_counts = assay(rlog(dds_processed))
rlog_counts = cbind(ID=rownames(rlog_counts), rlog_counts)
write.table(rlog_counts, file=snakemake@output$rlog, quote=FALSE, sep='\t', row.names=FALSE)

if (single_treatment) {
    comparison_results = as.data.frame(results(dds_processed, contrast=c("condition","treatment","control")))
    comparison_results = cbind(ID=rownames(comparison_results), symlog10baseMean = with(comparison_results, log10(baseMean + 1)), comparison_results)
    write.table(comparison_results, file=snakemake@output$results, quote=FALSE, sep='\t', row.names=FALSE)

    results_w_annot = merge(comparison_results, annots_subset, by=0)
    write.table(results_w_annot, file=snakemake@output$results_w_annot, quote=FALSE, sep='\t', row.names=FALSE)
} else {
    file.create(snakemake@output$results)
    file.create(snakemake@output$results_w_annot)
}
