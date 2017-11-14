

get_ipython().run_line_magic('load_ext', 'rpy2.ipython')

# Deseq analysis
get_ipython().run_cell_magic('R', '', "library('DESeq2')\nrpy2_deseq2 <- function(df,sample_types,types) {\n    sampleCondition <- unlist(sample_type)\n    df <- data.frame(experiment_df)\n    ddsHTSeq <- DESeqDataSetFromMatrix(df, DataFrame(sampleCondition), ~ sampleCondition)\n    colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$sampleCondition, levels=ctrl_treat)\n    dds<-DESeq(ddsHTSeq)\n    res<-results(dds)\n    res<-res[order(res$padj),]\n    deseq_df <- data.frame(res)\n    return(deseq_df)\n}")


def DESeq2(experiment_df, sample_type, types):
    get_ipython().run_line_magic('R', '-i experiment_df,sample_type,types')
    get_ipython().run_line_magic('R', 'deseq_df = rpy2_deseq2(experiment_df, sample_type, types)')
    get_ipython().run_line_magic('R', '-o deseq_df')
    return deseq_df


