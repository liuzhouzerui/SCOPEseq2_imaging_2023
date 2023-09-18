library("scran");
output_dir = 'D:/apoptosis/TS543_48h/analysis/mw/'
prefix = 'myeloid_30'
r = '2407868'

matrixfile = paste0(output_dir,'/',prefix,'.',r,'.ds.matrix.txt')
samplefile = paste0(output_dir,'/',prefix,'.',r,'.ds.samples.txt')
normfile = paste0(output_dir,'/',prefix,'.',r,'.ds.norm.txt')

df <- read.table(matrixfile, sep="\t");
data <- as.matrix(df[,3:ncol(df)]);
clusters <- as.matrix(read.table(samplefile));
sfs <- as.vector(sizeFactors(computeSumFactors(SingleCellExperiment(list(counts=data)),clusters=clusters)));
write.table(sfs,normfile,sep="\t",row.names=FALSE,col.names=FALSE);
