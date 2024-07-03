# we resume where SE_jacknife_prep.R ends, reload the result, and re-run the determinisitc preprocessing up to 
# subset generation 
args=commandArgs(trailingOnly=TRUE)
input_data= commandArgs(trailingOnly=TRUE)[[1]] # first arg determines variant
source('entropies.R')
source('helpers.R')
source('click-que.R')
filter_nonsignif=TRUE
f_name=paste0(input_data,'.csv')
read.table(f_name, sep = ',', header = TRUE)-> X
setwd(input_data)
# n samples per city
table(X$city)-> city_hist
# city of each sample
classVector<- X$city

#remove the class with 1 sample inside from first tool
if (any(city_hist==1)){
names(city_hist)[city_hist==1]-> cities2remove
removeIDX<- which(classVector %in% cities2remove)
classVector<- classVector[-removeIDX]
X <-X[-removeIDX, ]
}

#encode cities as integers
	#string names
	unique(classVector)->uq_y_v; 
	#number labels
	nuq_y_v<- seq_along(uq_y_v); 
	#placeholder
	ny<- rep(NA, length(classVector)); 
	for(u in seq_along(uq_y_v)) 
	{
	 ny[ classVector== uq_y_v[[u]] ] = nuq_y_v[[u]]
	}
	# table: city name, code
	classCityLookup<- cbind(uq_y_v, nuq_y_v)
	colnames(classCityLookup)<- c("city name", "code")
	classVector<-ny


#prepare matrix of markers-- main source for computatation
Ma <- as.matrix(X[,3:ncol(X)]) #Ma for markers

#if data is binary, remove markers with only 0 or only 1 values
unlapply<- function(...) unlist(lapply(...))

if (length(unique(Ma))==2) {
	unlapply(as.data.frame(Ma), 
	       function(x) length(unique(x))>1 
	       )-> okCols; Ma[,okCols]
}

J<-  as.integer(args[[2]]) # read the index of the J-th replicate from command line
#IMPORTANT: replace  actual sets by leave-0.25.samples-from-each-city-subsets

readRDS("LeaveQrtReps.rds") ->LeaveQrtReps
Ma <- Ma [ -LeaveQrtReps[[J]], ]
classVector <- classVector [  -LeaveQrtReps[[J]] ]
#Afterwards: processing is exactly the same as in city_similarity.R
#            but saved files differ

euc_cosine_mat<- function(mtrx){ #rows are vectors between which similarity is to be computed

prod_mtr<- mtrx %*% t(mtrx)
if (any(is.na(prod_mtr))) print("NA in prod_mtr")
self_prods<- diag(prod_mtr)
(self_prods==0) -> norm0mask
self_prods<- self_prods[!norm0mask]
prod_mtr<- prod_mtr[!norm0mask,
		    !norm0mask]
list(
S=prod_mtr/sqrt(self_prods %o% self_prods), #cosine similarity matrix
norm0mask=norm0mask) #  which rows (columns (samples)) were excluded
}

# cosine on SVD embedding of matrix 'Ma', using k first principal components
# k cannot be >= nrow(Ma)
SVDcos<- function(Ma,k= floor(0.3*min(nrow(Ma),ncol(Ma)))){
svd(Ma)-> UDV
UD<- UDV$u %*%  diag(UDV$d)
#watch out for this
if (any(is.na(UD))) print("NA in UD")
euc_cosine_mat(UD[,1:k])
}

# custom MDFS like feature selection but for decision having more than 2 classes
# works only for discrete X (tested for binary)
# uses entropies.R (provided by Krzysztof)
featureImportance<- function(y,X, getSignif=TRUE, method="holm", lvl=0.05){
n=ncol(X)
stopifnot(nrow(X)==length(y))
yX<- cbind(y,X)
indices<-compute.indices(2)
signifMask<-MI<-df.MI<-pv.MI<-matrix(NA,1,n)
PV<- vector(length=n)
for (j in 1:n) {
 counter<-compute.counter(as.data.frame(yX[, c(1,j) ]))
 entropies.cells<-compute.entropies.cells(counter,indices)
 entropies<-entropies.cells$entropies
 cells<-entropies.cells$cells
 MI[1,j]<-mutual.information(entropies,indices,c(0,1),c(1,0)) #mutual information between 1 and j variable
if (getSignif){
 df.MI[1,j]<-df.mutual.information(cells,indices,c(0,1),c(1,0)) 
 pv.MI[1,j]<-theoretical.pvalue(MI[1,j],df.MI[1,j])
	      }
} 
if (getSignif) {
PV<- pv.MI[1,]
pv.MI[1,] <- p.adjust(pv.MI[1,],method)
signifMask[1,] <- pv.MI[1,] < lvl
		}
list(MI=MI[1,],
     pv=pv.MI[1,],
     PV=PV,
     sgnf=signifMask[1,])
}


# feature importance comp
if (filter_nonsignif)
	Ma_imps <- featureImportance(classVector, Ma)


list(
if_svd=c(TRUE,FALSE)
    )-> variant_factors

if (filter_nonsignif)
	variant_factors$Ma_subset=c("all", "significant")
print("# signif")
print( length(Ma_imps$sgnf))
variant_factors$sparse=c(TRUE,FALSE)

#combination df describing each similarity variant
expand.grid(variant_factors) -> similVariants
# list to contain results of all variants
variousSimilarities<- as.list( 1: nrow(similVariants))
apply(apply(similVariants,1, function(ROW) paste0( names(variant_factors),'.', ROW) ),2,paste0, collapse='|')->vS_names
names(variousSimilarities)<- vS_names
for (v in 1:nrow(similVariants)) { #loop thru variants
if (!similVariants[v, 'sparse']){ #first pass: only non-sparse matrices( sparse are calc. based on non-sparse)
	Ma_v<- Ma
	if (filter_nonsignif)
		if (similVariants[v, 'Ma_subset']=='significant' ){
			Ma_v[,Ma_imps$sgnf] -> Ma_v   #choose only significant columns for decision
		}
	if (similVariants[v,'if_svd']){   #compute cos on  UD from SVD of samples or just on samples
		cosRes_v<- SVDcos(Ma_v)
		} else { cosRes_v<- euc_cosine_mat(Ma_v)    }
	if (any(is.na(cosRes_v$S))) print(paste("NA in cosine similarity! At:\n", vS_names[[v]] ))
	variousSimilarities[[v]]<-cosRes_v   
	variousSimilarities[[v]]$S<-cosRes_v$S +1    # [-1,1] -> [0,2]
	}
}

count_cliques<- function(W, tholds, minClqSize=3, ...) {
clq_counts<-rep(0, length(tholds))
diag(W)<-0
W.t<-W
for (tr in seq_along(tholds)) {
	TR<- tholds[[tr]]
	W.t[W.t < TR] = 0
	tmp_clq<- greedyCliquesWeighted(W.t,...)
	tmp_clq<-zeroOutSmall(tmp_clq, minClqSize)
	tmp_tab<-fastTable(tmp_clq)
	clq_counts[[tr]]<-sum(tmp_tab$value!=0)
			      }
which_maxCount<-which.max(clq_counts)
maxCount<- clq_counts[[which_maxCount]]
finThr<- tholds[[which_maxCount]]
W.t<-W
W.t[W.t < finThr ] = 0
prop_clqs<- greedyCliquesWeighted( W.t,...)
list(maxCount=maxCount,
     used_thr=finThr,
     clq_mem=prop_clqs)
}


lapply(seq_along(variousSimilarities), function(v) 
			{ 
			if (!similVariants[v, 'sparse']){
			x<- variousSimilarities[[v]]
			seq(from=min(x$S),
			    to= max(x$S),
			    length.out=50)  #look for thr in 50 steps
				} else {NA}
			}
	)-> thr_ranges;
 
for (v in 1:nrow(similVariants)) { #loop thru variants, add sparsified versions
if (similVariants[v, 'sparse']){
	cols2match <- which(colnames(similVariants)!='sparse')
	vals2match <- similVariants[v, ]
	setdiff(which(
	Reduce('&',
	lapply(cols2match, function(j) similVariants[,j]==vals2match[[j]])
	     )
	),v)-> v2 #get the row number of the non-sparse variant	
		   #and feed it to the function below
	count_res<-count_cliques(variousSimilarities[[v2]]$S, thr_ranges[[v2]], extended_mode=FALSE, mode="average")
	S_thr<- variousSimilarities[[v2]]$S #zero out weak edges
	S_thr[ S_thr < count_res$used_thr ]= 0
	variousSimilarities[[v]] <- variousSimilarities[[v2]]
	variousSimilarities[[v]]$S <- S_thr #substitute thresholded
	}
}

###Statistics
### avg(S_in) - avg(S_out)
###function to compute the statistic 
###note: zero-norm samples should be removed both from classVector and simil_mat
S_io_diff<- function( simil_mat, classVector) {
	stopifnot( nrow(classVector)==nrow(simil_mat) )
	classes<- unique(classVector)
	print(table(classVector))
	class_places<- lapply( classes, function(cl)
				 classVector==cl
	      		     )
	avg_S_in<- mean(do.call(c, lapply(class_places, 
			function(cp) { block<- simil_mat[cp,cp, drop=FALSE]
				       if (nrow(block)==1) block[[1]] else mean(block[lower.tri(block)])
				     } )
			)
		    )
	off_diag_simils<- list()
	for (i in seq_along(class_places)) 
	for (j in seq_along(class_places))
		if (i < j){ cp_i<- class_places[[i]]; cp_j<- class_places[[j]]
	   block<-simil_mat[ cp_i, cp_j, drop=FALSE ]
	   off_diag_simils[[ 
	     length(off_diag_simils) + 1 
			  ]] = mean(block)
		       	  }
	avg_S_out<- mean(do.call(c,off_diag_simils))
	avg_S_in - avg_S_out
}

###Computing the statistic for all variants of similarity matrix
vS_stats<-vector()
for (v in seq_along(variousSimilarities))
{
variousSimilarities[[v]]$norm0mask -> n0m
vS_stats[[v]]= S_io_diff( simil_mat= variousSimilarities[[v]]$S,
			  classVector= classVector[ !n0m]) #exclude samples 
							   #not in S
}


cbind(names(variousSimilarities), vS_stats) -> S_io_df
colnames(S_io_df)<- c('variant','mean(S_in) - mean(S_out)')
S_io_df<- as.data.frame(S_io_df)
S_io_df[,'n_Samples'] = unlapply( 
			variousSimilarities, function(x) sum(!x$norm0mask)
				   )


write.csv( S_io_df, paste0('similarity_stats_results_rep',
			    J,
			    '.csv'
			   )
	  )

