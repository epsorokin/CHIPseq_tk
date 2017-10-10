require(edgeR)
#cds <- estimateDispersions( cds , method = "per-condition", sharingMode="fit-only")
#cds <- estimateDispersions( cds , method = "per-condition", fitType="local")
cds <- estimateDispersions( cds , method = "pooled", fitType="local")
res <- nbinomTest (cds,"exp","con")
resOrdered <- res [order (res$padj),]
res <- na.exclude(res)
resSig <- res[ res$padj < 0.05, ]
dim(resSig)
resSig
cds <- estimateDispersions( cds , method = "per-condition", fitType="local")
res <- nbinomTest (cds,"exp","con")
res[ res$padj < 0.05, ]
cds <- estimateDispersions( cds , method = "pooled")
cds <- estimateDispersions( cds , method = "pooled",fitType="local",sharingMode="maximum")
res <- nbinomTest (cds,"exp","con")
res[ res$padj < 0.05, ]

sig <- res[ res$padj < 0.05, ]

write.table(res,"deseq_out_nicd.txt",sep="\t")
