
inputs= c('amrfinder_binary', 'amrpp_binary', 'bowtie_binary', 'rgi_binary')
fpat<-'%s/similarity_stats_results_rearg%d.csv'
n_rearr=5000
reargs<-list()
for (in_name in inputs){
	  reargs[[in_name]] = list()
	for (R in 1:n_rearr) {
	    reargs[[in_name]][[R]]= read.csv( sprintf(fpat, in_name, R) )[,3]
	}
}

read.csv( sprintf(fpat, in_name, R) )[,2] -> LEGEND_r
fpat<-'%ssimilarity_stats_results.csv'
point_estimates<- list()
for (in_name in inputs)
	point_estimates[[in_name]]=  read.csv( sprintf(fpat,in_name))

LEGEND= point_estimates[[1]][,2]
n_samples= do.call(cbind,lapply( point_estimates, function(x) x[,4] ))
for (in_name in inputs)
	point_estimates[[in_name]]=  point_estimates[[in_name]][,3]

all_datasets_point_estimates<- do.call(cbind, point_estimates)
colnames(all_datasets_point_estimates) = inputs
all_datasets_reps= list()
for (R in 1:n_rearr)
	all_datasets_reps[[R]] = do.call(cbind,  lapply(inputs, function(in_name) reargs[[in_name]][[R]])
					)
##only take those with all variables
all_datasets_point_estimates<- all_datasets_point_estimates[c(1,2,5,6), ]
pv<- all_datasets_point_estimates
pv[,]=0
for (R in 1:n_rearr)
	pv = pv + (all_datasets_reps[[R]] >= all_datasets_point_estimates)
pv<- pv/n_rearr
titles=c('AMRFinderPlus','AMR++','Bowtie','RGI')
colnames(pv)=titles
pv[,]= p.adjust(pv, method="holm")
rownames(pv)=LEGEND_r
write.csv(pv,'H_corrected_pvs.csv')
