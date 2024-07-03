library(dplyr)


# Set the working directory to the script's directory
#script_path <- normalizePath(sys.frame(1)$ofile)
setwd('/Users/rodolfotoscan/Documents/estudos/UJ/phd/projects/camda_2023/publication/0_writing/data/normalization')

ssu_std = read.csv("table_std.tsv",sep="\t")
ssu_cov50 = read.csv("table_mincov_50.tsv",sep="\t")
assemblies = read.csv("assemblies_stats_full-output.tsv",sep="\t")
reads =  distinct(read.csv("camda_reads_bases.tsv",sep="\t"))
df <- merge(merge(merge(ssu_std,ssu_cov50,by='sample'),assemblies,by='sample'),reads,by='sample')
#write.csv(df,file="metrics_us_samples.csv")
df <- subset(df,select = -c(A...,T...,G...,C...,N50,Max,Min,N...,N....1,Gap...,Uncertain.bp.,Ncount))

# Calculate the correlation matrix
#


#plot(df$count_std, df$count_cov50, type = "point", pch = 19,      col = "red", xlab = "x", ylab = "y")
norm_metrics <- subset(df, select = -c(Contigs, bases) )


rownames(df) <- df$sample
df <- subset(df,select=-sample)


#cor_matrix <- cor(df)

rm(assemblies,reads,ssu_cov50,ssu_std,df)
