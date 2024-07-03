source('helpers.R')

prepare_rearrangements<- function(classVector, n_rearrangements){
cl_vals<- sort(unique(classVector))
n_cl<-  unlapply(cl_vals, function(clv) sum(classVector==clv))
L_samples_cl<- lapply(cl_vals, function(clv) which(classVector==clv))
samples=1:length(classVector)
rearrangements<- list()
for (J in 1:n_rearrangements) {
	print(J)
	# get a permutation of integers 1:N
	# assign label A to first N_A, label B to next N_B, C to next N_C... etc
	sample(samples, length(samples), replace=FALSE)-> permut
	Jth_potential <- classVector[ permut ]
	if (length(rearrangements)==0){ rearrangements[[J]]<- Jth_potential}else {
	 already_present= TRUE
	while(already_present) {
		 already_present= any( unlist(lapply(rearrangements, function(rng) all(rng==Jth_potential)) ))
		 if (already_present)
			{
			sample(samples, length(samples), replace=FALSE)-> permut
			Jth_potential <- classVector[ permut ]
			}
		}
	Jth_potential-> rearrangements[[J]]
									 }
	
	}
return(rearrangements)
}

