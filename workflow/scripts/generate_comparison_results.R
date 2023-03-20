library("DESeq2")

counts = read.csv(file=snakemake@input$counts, sep='\t')
annots = read.csv(file=snakemake@input$annotations, sep='\t')

organism = snakemake@wildcards$organism
comparison = snakemake@wildcards$comparison

conditions = names(snakemake@config$comparisons[[organism]][[comparison]])
single_treatment = length(conditions) == 2 & 'control' %in% conditions & 'treatment' %in% conditions

sample_use_df = stack(data.frame(snakemake@config$comparisons[[organism]][[comparison]]))

rownames(sample_use_df) = sample_use_df$values
sample_use_df = subset(sample_use_df, select=-c(values))

counts_subset = subset(counts, organism==snakemake@wildcards$organism)
rownames(counts_subset) = counts_subset$ID
counts_subset = subset(counts_subset, select=rownames(sample_use_df))

annots_subset = subset(annots, organism==snakemake@wildcards$organism)
rownames(annots_subset) = annots_subset$ID

# organism must be present in all samples. Good test to tell if it is from experience...
if (min(colMeans(counts_subset[sapply(counts_subset, is.numeric)],na.rm=TRUE)) < 1){
    stop('organism must be present in all samples')
}

dds = DESeqDataSetFromMatrix(countData = counts_subset, colData = sample_use_df, design = ~ ind)
dds_processed = DESeq(dds)

vst_counts = assay(vst(dds_processed))
vst_counts = cbind(ID=rownames(vst_counts), vst_counts)
write.table(vst_counts, file=snakemake@output$vst, quote=FALSE, sep='\t', row.names=FALSE)

rlog_counts = assay(rlog(dds_processed))
rlog_counts = cbind(ID=rownames(rlog_counts), rlog_counts)
write.table(rlog_counts, file=snakemake@output$rlog, quote=FALSE, sep='\t', row.names=FALSE)

if (single_treatment) {
    dds$ind = relevel(dds$ind, ref = "control")
    comparison_results = as.data.frame(results(dds_processed))
    comparison_results = cbind(ID=rownames(comparison_results), symlog10baseMean = with(comparison_results, log10(baseMean + 1)), comparison_results)
    write.table(comparison_results, file=snakemake@output$results, quote=FALSE, sep='\t', row.names=FALSE)

    results_w_annot = merge(comparison_results, annots_subset, by=0)
    write.table(results_w_annot, file=snakemake@output$results_w_annot, quote=FALSE, sep='\t', row.names=FALSE)
} else {
    file.create(snakemake@output$results)
    file.create(snakemake@output$results_w_annot)
}
