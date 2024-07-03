args= commandArgs(trailingOnly=TRUE)
args[[1]]-> var_name
# cli arg can be either 'name.of.the.variant' or number of the variant
if( nchar(var_name)==1) var_name= as.integer(var_name)

readRDS('vS.rds')->vS #similarities
readRDS('cV.rds')->cV #class Vector
ccodes<-read.csv('CityCodes.csv')
cV<- cV[ !vS[[var_name]]$norm0mask ]
#perm<- sample(1:length(cV), length(cV), replace=FALSE)
#cV<- cV[ perm ]
W<- vS[[ var_name ]]$S
#W<- W[ perm, perm ]
W<- W[ order(cV), order(cV) ]
cV<- cV [ order(cV) ]
cN<- cV
for (uqv in unique(cV) )
	cN[ cV==uqv ] = ccodes[ ccodes$code==uqv,2]
colnames(W)<-rownames(W)<- cN
labelColors<-c("red","green","blue","yellow","magenta","black")
	 
jpeg( paste0('correct',var_name,'.jpg'), 550,550 )
heatmap(W, Rowv=NA, Colv=NA, revC=TRUE, scale="none",
	RowSideColors= labelColors[cV],
	ColSideColors= labelColors[cV]
	)
dev.off()

