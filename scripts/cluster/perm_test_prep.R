
source('perm_test.R')
#table (sample, makrer + additional info_
input_data= commandArgs(trailingOnly=TRUE)[[1]]
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

set.seed(123)
rearrgs_1to4<- list()   # 4 different rearrangements since for each dataset 4 similarity variants to be tested.
for (j in 1:4) rearrgs_1to4[[j]]=prepare_rearrangements(classVector=classVector, n_rearrangements=5000)

saveRDS(rearrgs_1to4, "rearrangements_1to4.rds")
