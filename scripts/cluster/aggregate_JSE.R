
inputs= c('amrfinder_binary', 'amrpp_binary', 'bowtie_binary', 'rgi_binary')
fpat<- '%s/similarity_stats_results_rep%d.csv'	
n_reps=500
reps<-list()
for (in_name in inputs){
	  reps[[in_name]] = list()
	for (J in 1:n_reps) {
	    reps[[in_name]][[J]]= read.csv( sprintf(fpat, in_name, J) )	[,3]
	
	}
}

fpat<-'%ssimilarity_stats_results.csv'
point_estimates<- list()
for (in_name in inputs)
	point_estimates[[in_name]]=  read.csv( sprintf(fpat,in_name))

LEGEND= point_estimates[[1]][,2]
n_samples= do.call(cbind,lapply( point_estimates, function(x) x[,4] ))
for (in_name in inputs)
	point_estimates[[in_name]]=  point_estimates[[in_name]][,3]

jackknife_SE<- function(replicates,r){
    
    replicate_mean<- Reduce('+', replicates)/length(replicates)
    sq_devs<- lapply(replicates, function(X) (X- replicate_mean)^2 )
    sq_devs_sum<- Reduce('+', sq_devs)
    constant_factor=  r/length(replicates)
    SE<- sqrt(constant_factor*sq_devs_sum)  
    return(SE)
}

all_datasets_point_estimates<- do.call(cbind, point_estimates)
colnames(all_datasets_point_estimates) = inputs
all_datasets_reps= list()
for (J in 1:n_reps)
	all_datasets_reps[[J]] = do.call(cbind,  lapply(inputs, function(in_name) reps[[in_name]][[J]])
					)

SE_all_datasets<- jackknife_SE( all_datasets_reps, r=4) # 4 because we leave 0.25 samples out

for (J in 1:n_reps)
	all_datasets_reps[[J]][ n_samples <60 ]  = NA
all_datasets_point_estimates[ n_samples< 60 ] = NA
SE_all_datasets[ n_samples < 60 ] = NA #when the number of samples with non zero norm is too small we leave variant out
				  #in particular such variants are amrfider based with unsignificant marker filtering

plot_order=c(5,1,6,2,7,3,8,4)
SE_all_datasets<- SE_all_datasets[plot_order,]
all_datasets_point_estimates<- all_datasets_point_estimates[plot_order,]
for (J in 1:n_reps)
	all_datasets_reps[[J]]=all_datasets_reps[[J]][plot_order,]
LEGEND=LEGEND[plot_order]

##PLOTS
titles=c('AMRFinderPlus','AMR++','Bowtie','RGI')
jpeg('SE.jpg', width=900, height= 500)
par(mfrow=c(2,4))
for (ROW in 1:2) for (COL in 1:4){
plot(rnorm(n_reps,0,0.05), sapply(all_datasets_reps,
			      function(x) x[(ROW-1)*4 + 1,COL] ),
	type="point", main= paste0(titles[[COL]],'\n', ifelse(ROW==1,'all variables', 'significant MI variables')),
	pch=16, col="darkgoldenrod2",
	ylim=c(-0.1,0.8), xlim=c(-1,4),
	xaxt="n", ylab=NA, xlab=NA)
for (k in 2:4)
points(rnorm(n_reps,k-1,0.05), sapply(all_datasets_reps,
			      function(x) x[(ROW-1)*4 + k,COL] ),
	     pch=16, col=ifelse(k<3,'darkgoldenrod2','cornflowerblue'))
points(0:3, all_datasets_point_estimates[ ((ROW-1)*4 + 1):((ROW-1)*4 + 4),COL ], pch=4, cex=3)
arrows(x0=0:3, 
y0=all_datasets_point_estimates[ ((ROW-1)*4 + 1):((ROW-1)*4 + 4),COL ]-SE_all_datasets[ ((ROW-1)*4 + 1):((ROW-1)*4 + 4),COL ],
       x1=0:3, 
y1=all_datasets_point_estimates[ ((ROW-1)*4 + 1):((ROW-1)*4 + 4),COL ]+SE_all_datasets[ ((ROW-1)*4 + 1):((ROW-1)*4 + 4),COL ],
code=3, col="black", lwd=2,angle=90, length=0.1)
if (!((ROW==2)&(COL==1)))
legend("topright", legend=c("SVD","plain"), col=c("darkgoldenrod2", "cornflowerblue"), pch=16 )
if (((ROW==2)&(COL==1)))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
}
dev.off()

