source('helpers.R')
unlapply<-function(...) unlist(lapply(...))
leaveFracOutMulticlassAllSubsets<- function( classVector, frac, n_replicates=200){
# function should not be used in practice, since it generates all possible subsets and samples from them
# number of all possible subsets to generate is probably unfeasibly large
# specify fraction of samples to leave out from each class
# and then generate every possible configuration of samples where each sample has 'frac' samples of each class not included
cl_vals<- sort(unique(classVector))
n_cl<-  unlapply(cl_vals, function(clv) sum(classVector==clv))
samples<- seq_along(classVector)
n_cl_2remove<- floor(frac*n_cl)
stopifnot(all(n_cl_2remove>0))
n_cl_reduced<- n_cl - n_cl_2remove

L_samples_cl<- lapply(cl_vals, function(clv) which(classVector==clv))
L_samples2remove_cl<- lapply(seq_along(L_samples_cl), function(cl)
					combn(L_samples_cl[[cl]], n_cl_2remove[[cl]],simplify=FALSE )
			    )
L_nSubsets_cl<- lapply(L_samples2remove_cl, function(L_subsets) length(L_subsets) )

do.call(expand.grid, lapply(L_nSubsets_cl, function(nSubsets_cl) 1:nSubsets_cl ) ) -> subsetIDXcombinations

combinations_chosen<- sample(1:nrow(subsetIDXcombinations), n_replicates, replace=FALSE)

subsetIDXcombinations<- subsetIDXcombinations[combinations_chosen,]
replicates<-list()
for (i in 1:nrow(subsetIDXcombinations)){
	subsetIDXcombinations[i,]-> which_subset_from_each_cl
	do.call(c,lapply(seq_along(which_subset_from_each_cl), function(j)
			L_samples2remove_cl[[j]][[ which_subset_from_each_cl [[j]] ]]))-> ith_2remove
	samples[-ith_2remove] -> replicates[[ length(replicates) + 1 ]] 
					 }
return(replicates)
}

leaveFracOutMulticlassSampleEachSubset<- function( classVector, frac, n_replicates=200){
# this function shall be used in practice, since it generates subsets element by element
# instead of sampling from preprepared set of all possible subsets
# specify fraction of samples to leave out from each class
# and then generate every possible configuration of samples where each sample has 'frac' samples of each class not included
cl_vals<- sort(unique(classVector))
n_cl<-  unlapply(cl_vals, function(clv) sum(classVector==clv))
samples<- seq_along(classVector)
n_cl_2remove<- floor(frac*n_cl)
stopifnot(all(n_cl_2remove>0))
n_cl_reduced<- n_cl - n_cl_2remove

L_samples_cl<- lapply(cl_vals, function(clv) which(classVector==clv))

replicates<- list()
for (J in 1:n_replicates) {

	Jth_potential<- unlist(lapply(seq_along(L_samples_cl), function(cl_nr)
					sample( L_samples_cl[[cl_nr]], n_cl_2remove[[cl_nr]], replace=FALSE))
		    )
	if (length(replicates)==0){ replicates[[J]]<- Jth_potential}else {
	 already_present= TRUE
	while(already_present) {
		 already_present= any( unlist(lapply(replicates, function(rpl) all(rpl==Jth_potential)) ))
		 if (already_present)
			{
			Jth_potential<- unlist(lapply(seq_along(L_samples_cl), function(cl_nr)
					sample( L_samples_cl[[cl_nr]], n_cl_2remove[[cl_nr]], replace=FALSE))
					      )
			}
		}
	Jth_potential-> replicates[[J]]
									 }
	
	}
return(replicates)
}
