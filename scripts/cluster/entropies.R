

#################################################################
#Counting objects in each combination of classes 

#a) the fast way using a table function

#x<-as.data.frame(x)
compute.counter<-function(x #x MUST BE a data frame
                          ) {
 return(table(x))
} 

#b) manually to be reproduced in c code; slower in R, NOT for current use

#for (i in 1:ncol(x)) x[,i]<-as.numeric(as.factor(x[,i]))-1
compute.counter.manually<-function(x #x MUST contain numbers from 0 without gaps
                          ) {
 dimension<-ncol(x)
 ncells<-colMaxs(x)+1
 counter<-array(0,dim=ncells)
 for (i in 1:nrow(x)){
  index<-x[i,dimension]
  for (j in (dimension-1):1) index<-index*ncells[j]+x[i,j]
  index<-index+1
  counter[index]<-counter[index]+1
 }
 return(counter)
} 

#################################################################
#Computing indices for all the combination of n variables
#To be run once at the beginning

#a) using expand.grid R function

compute.indices<-function(dimension) {
 
 return(t(!expand.grid(rep(list(0:1),dimension)))) #from 1,1,..1 to 0,0,..0
}

#b) manually to be reproduced in c code; slower in R, slightly different order, NOT for current use
compute.indices.manually<-function(dimension) {
 indices<-rep(1,dimension)
 index.all<-1:dimension
 for (d in 1:(dimension-1)) {
  indices.d<-combn(index.all,dimension-d)
  for (i in 1:ncol(indices.d)) indices<-cbind(indices,(1:dimension)%in%indices.d[,i])
 }
 return(cbind(indices,rep(0,dimension)))
} 


######################################################################
#Computing entropies from the counters

compute.entropies.cells<-function(counter,indices,pseudo.count=0.25) {
 dimension<-length(dim(counter))
 counter<-counter+pseudo.count

 entropies<-sum(counter*log(counter))
 cells<-length(counter)

 for (i in 2:(ncol(indices)-1)) {
  index<-c(which(!indices[,i]),which(indices[,i]))
  subcounter<-aperm(counter,index)         
  for (j in 1:sum(!indices[,i])) subcounter<-colSums(subcounter)  
  entropies[i]<-sum(subcounter*log(subcounter))
  cells[i]<-length(subcounter)  #so far, we don't take zeros into account      
 }
 
 entropies[length(entropies)+1]<-sum(counter)*log(sum(counter))
 cells[length(cells)+1]<-1

 return(list(entropies=entropies,cells=cells))
} 
 
#######################################################################
# Various entropy functions

#entropy of a subset S
total.entropy<-function(entropies,indices,S #S MUST BE given as bit mask
                        ) {
 return(entropies[length(entropies)]-
        entropies[colSums(indices!=S)==0])
}

#entropy of a subset S1 conditioned on |S2
conditional.entropy<-function(entropies,indices,S1,S2 #S1, S2 as bit masks
                              ) {
 return(entropies[colSums(indices!=S2)==0]-
        entropies[colSums(indices!=(S1|S2))==0])
}

#mutual information between S1 and S2
mutual.information<-function(entropies,indices,S1,S2 #S1, S2 as bit masks
                             ) {
 return(entropies[colSums(indices!=(S1|S2))==0]+
        entropies[length(entropies)]-
        entropies[colSums(indices!=S1)==0]-
        entropies[colSums(indices!=S2)==0]
 )
}

#mutual information between S1 and S2 conditioned on |S3
conditional.mutual.information<-function(entropies,indices,S1,S2,S3 #S1,S2,S3 as bit masks
                                         ) {
 return(entropies[colSums(indices!=(S1|S2|S3))==0]+
         entropies[colSums(indices!=S3)==0]-
         entropies[colSums(indices!=(S1|S3))==0]-
         entropies[colSums(indices!=(S2|S3))==0]
 )
}

#interaction information 
interaction.information<-function(entropies,indices) {
 return(sum((-1)^rowSums(indices)*entropies))     
}

###########################################################################
#number of degrees of freedom

df.mutual.information<-function(cells,indices,S1,S2) {
 return((cells[colSums(indices!=S1)==0]-1)*
        (cells[colSums(indices!=S2)==0]-1))
}

df.conditional.mutual.information<-function(cells,indices,S1,S2,S3) {
 return((cells[colSums(indices!=S1)==0]-1)*
        (cells[colSums(indices!=S2)==0]-1)*
        cells[colSums(indices!=S3)==0])
}

###########################################################################
#theoretical p-value 

theoretical.pvalue<-function(information,df) return(pchisq(2*information,df=df,lower.tail=F,log.p=F))

log.theoretical.pvalue<-function(information,df) return(pchisq(2*information,df=df,lower.tail=T,log.p=T))



############################################################################
#Example use cases


#n<-100 #number of objects
#p<-10 #number of variables
#
#method='holm' #correction for multiple tests
#threshold=0.05 #significance level
#
##example discrete data
#x<-matrix(floor(4*runif(n*p)),n,p)
#
##data preparation
#for (i in 1:ncol(x)) x[,i]<-as.numeric(as.factor(x[,i]))-1
#x<-as.data.frame(x)
#
##############################################################
##mutual information between pairs of x
#indices<-compute.indices(2)
#MI<-df.MI<-pv.MI<-matrix(NA,p,p)
#
#for (i in 2:p) for (j in 1:(i-1)) {
# counter<-compute.counter(x[,c(i,j)])
# 
# entropies.cells<-compute.entropies.cells(counter,indices)
# entropies<-entropies.cells$entropies
# cells<-entropies.cells$cells
# 
# MI[i,j]<-mutual.information(entropies,indices,c(0,1),c(1,0)) #mutual information between first and second variable
# df.MI[i,j]<-df.mutual.information(cells,indices,c(0,1),c(1,0)) 
# pv.MI[i,j]<-theoretical.pvalue(MI[i,j],df.MI[i,j])
#} 
#
##pairs significantly connected
#significant.vec<-p.adjust(pv[lower.tri(pv.MI,diag=F)],method)<threshold
#significant.MI<-pv.MI
#significant.MI[lower.tri(significant.MI,diag=F)]<-significant.vec
#
#
######################################################################
##conditional mutual information between pairs of x given x[,1]
#indices<-compute.indices(3)
#CMI<-df.CMI<-pv.CMI<-matrix(NA,p,p)
#
#for (i in 3:p) for (j in 2:(i-1)) {
# counter<-compute.counter(x[,c(i,j,1)])
# 
# entropies.cells<-compute.entropies.cells(counter,indices)
# entropies<-entropies.cells$entropies
# cells<-entropies.cells$cells
# 
# CMI[i,j]<-conditional.mutual.information(entropies,indices,c(1,0,0),c(0,1,0),c(0,0,1)) 
#                              #mutual information between first and second variable given third
# df.CMI[i,j]<-df.conditional.mutual.information(cells,indices,c(1,0,0),c(0,1,0),c(0,0,1)) 
# pv.CMI[i,j]<-theoretical.pvalue(CMI[i,j],df.CMI[i,j])
#} 
#
##pairs significantly connected
#significant.vec<-p.adjust(pv.CMI[lower.tri(pv.CMI,diag=F)],method)<threshold
#significant.CMI<-pv.CMI
#significant.CMI[lower.tri(significant.CMI,diag=F)]<-significant.vec
#
#
