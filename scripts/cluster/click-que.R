#i use this for more readable code
max2minOrder<- function(x)
{
  stopifnot(is.numeric(x))
  order(-x)
}


#additional stuff for correlation matrix:
#find correlation thresholding value based on the 
#giant component size of the induced network
phi.TH.search<- function(W, phi=(1+sqrt(5))/2, minW=0.024,
                         giantComponentBound=0.9*ncol(W),
                         maxTH= 0.9 )
{
  stopifnot(all(diag(W)==0) )
  W.values<- W [lower.tri(W)]
  W.values<- W.values[ W.values > minW]
  W.values<- W.values[order(W.values)]
  N<-length(W.values)
  n=as.integer(N/(1+phi))
  m= N - n
  thIDX<- n
  TH<- W.values[[thIDX]]
  while (TH <= maxTH)
  {
    A<- (W > TH)*1.0
    G<-graph_from_adjacency_matrix(A)
    G_components<-components(G)
    gcSize<- G_components$csize[[1]]
    if (gcSize <= giantComponentBound)
      return(
        list(A=A,
             comps=G_components,
             TH=TH)
      )
    else
    {
      N<- m
      n=as.integer(N/(1+phi))
      m= N - n
      thIDX<- thIDX+ n
      TH<- W.values[[thIDX]]
    }
    message(gcSize)
    message(TH)
  }
  warning("did not converge, reached maxTH cutoff threshold.")
  return(
    list(A=A,
         comps=G_components,
         TH=TH)
  )
  
}

#clique finding

fastTable<- function(values)
{
    uqs<-unique(values)
    histo<-unlist(lapply(uqs, function(x) sum(values==x)))
    list(value=uqs,
         count=histo)
}

greedyCliquesWeighted<- function(W,extended_mode=FALSE, mode="basic",
                                return_list=FALSE){
  if(is.null(dim)) return(1)
  stopifnot(all(diag(W)==0))
  N<-ncol(W)
  if(is.null(N)) return(1)
  if(N==1) return(1)
  if(N==0) return(1)
  degs<-colSums(W)
  cl_memship<- rep(0,N)
  if (return_list)
    if (extended_mode) {warning("for extended mode return_list forced to FALSE"); return_list=FALSE}
  if (return_list)
    cl_list<-list()
  clID<-0
  Zeros<- cl_memship==0
  colnames(W)<- paste0(1:N,"n")
  rownames(W)<- paste0(1:N,"n")
  Zero2Zero<- W[ Zeros, Zeros ] # edges joining nodes not assigned anywhere yet
  progress_counter=0
  while(any(Zero2Zero > 0)) #if yes, there is some edge that can be new cluster
  {
    progress_counter= progress_counter + 1
    w_max<- max(Zero2Zero)
    ij_max<-which(Zero2Zero==w_max,arr.ind=TRUE)[1,]  #  max_e1 -  max_e2
    max_e1<-which(rownames(W)==rownames(Zero2Zero)[ij_max[[1]]])
    max_e2<-which(rownames(W)==colnames(Zero2Zero)[ij_max[[2]]])
    clID=clID+1  #new cluster creation
    cl_memship[ c(max_e1, max_e2) ] = clID
    if (return_list) cl_list[[clID]]<- c(max_e1, max_e2)
    if ( degs[[max_e1]] >= degs[[max_e2]] )
      current.node=max_e1
    else
      current.node=max_e2 
      potential<- which( (W[,current.node]!=0) &
                             (cl_memship==0)&
                        (rowSums(W[,cl_memship==clID,drop=FALSE])>0)
                             ) #consider only:
    #- nodes connected to current
    #- which are not clustered.
    if (length(potential))
    {
      #this is basic mode
      #BASIC MODE START
      #"add a thing that boss likes the most"
      if (mode=='basic'){   #modes ifelse start
      neigh.W<- W[current.node,potential]
      neighDecrOrder<- max2minOrder(neigh.W)
      #print(neighDecrOrder)
      potential<- potential[neighDecrOrder]
      for (POTENTAT in potential)
      {
        #                       print(POTENTAT)
        new.cl<-c(POTENTAT,which(cl_memship==clID))
        W.new.cl<-W[new.cl,new.cl]
        W.new.cl.eIDX<- lower.tri(W.new.cl)
        W.new.cl.edges<- W.new.cl[ W.new.cl.eIDX]
        #                       print(any(W.new.cl.edges==0))
        if (any(W.new.cl.edges==0)) #dont join
          invisible() #do nothing
        else {                       #join
          cl_memship[[POTENTAT]]=clID
          if (return_list) cl_list[[clID]]<- c( cl_list[[clID]], POTENTAT)
            }
      }
      #BASIC MODE END
      } else if (mode %in% c('average', 'maxOut')){
      #average mode:
      #"add a thing that is on average the most likable in the clique"
      #each time before something is joined...
      #...compute average over connections' weights...
      #... to whole clique. Join the thing that has the largest avg.
      #maxOut mode:
      #"add a thing that has the highest potential of extending the clique"
      #each time before something is joined...
      #...compute number of connections to the nodes outside...
      #... of clique. Join the thing that has the largest such number.
      while(length(potential))
      {
      CL_SIZE<- sum(cl_memship==clID)
      cliquePerserversMask<- rowSums(W[, cl_memship==clID,drop=FALSE]!=0)==CL_SIZE
      if (any(cliquePerserversMask))
        potential <- which( cliquePerserversMask&(cl_memship==0))
      else
        potential<-NULL
      if (length(potential))
         {
         if (mode=='average')
            {
            POTENTATpos<- which.max(rowMeans(W[potential,cl_memship==clID,drop=FALSE]))
            if(length(potential)==0) print(potential)
            POTENTAT<- potential[[POTENTATpos]]
            cl_memship[[POTENTAT]]=clID
            if (return_list) cl_list[[clID]]<- c( cl_list[[clID]], POTENTAT)
            }
          if (mode=='maxOut')
            {
            POTENTATpos<- which.max(rowSums(W[potential,cl_memship==0,drop=FALSE]>0))
            POTENTAT<- potential[[POTENTATpos]]
            cl_memship[[POTENTAT]]=clID
            if (return_list) cl_list[[clID]]<- c( cl_list[[clID]], POTENTAT)
            }
         }
       }
      } else { stop("wrong mode argument")} #modes if else end
    } #if !is.null(potential) end
    Zeros<- cl_memship==0
    #extension
    if (extended_mode)
    {
    continue_flag=TRUE#set continue_flag as TRUE
    emeregency_stopper=0
    while (continue_flag){#while continue_flag
    clSIZE=sum(cl_memship==clID)
    PossibleSwaps= (rowSums(W[,cl_memship==clID]>0)==(clSIZE-1))&
                    (cl_memship==0)
    if (!any(PossibleSwaps)) break
    #this ugly indexing is for the purpose of having indexes relative to original network size and order
    SwapForWhat<-(1:nrow(W))[cl_memship==clID][max.col(-(W[PossibleSwaps,cl_memship==clID]>0))] 
    SwapForWhatFreq<- fastTable(SwapForWhat) 
    SwappabilityOrder<- max2minOrder(SwapForWhatFreq[['count']])
    SwappableSeq<- SwapForWhatFreq[['value']][SwappabilityOrder]
    #print(SwappableSeq)
    for_broken=FALSE
    if (length(SwappableSeq) > 1){
    for (potentialInClusterSwap in SwappableSeq) #note this loop 
        {                               # gets broken immidiately
                                        # when potentialInClusterSwap
                                        #is good enough.
                                        #but it continues if it was not 
    #    print(potentialInClusterSwap)
        PosSwpPos<- which(PossibleSwaps)
        whichPosSwp<- which(SwapForWhat==potentialInClusterSwap)
        pocr<-rep(FALSE,nrow(W))
        pocr[PosSwpPos[whichPosSwp]]<- TRUE
    #    print((W[pocr,pocr]>0)*1)
        PotentialOutClusterReplacement<-((W[potentialInClusterSwap,]==0)&pocr)
    #    print((W[(cl_memship==clID)|PotentialOutClusterReplacement,(cl_memship==clID)|PotentialOutClusterReplacement]>0)*1)
    #    print((W[(cl_memship==clID)|pocr,(cl_memship==clID)|pocr]>0)*1)
    #    print('aux mat and cliques')
    #    print(W[PotentialOutClusterReplacement,PotentialOutClusterReplacement]>0)
        aux_cliques<-greedyCliquesWeighted(W[PotentialOutClusterReplacement,PotentialOutClusterReplacement],
                                            mode=mode)
        aux_cliques<-uniqueSingletonLabels(aux_cliques)
        auxHist<-fastTable(aux_cliques)
        largestAuxClique<- auxHist$value[which.max(auxHist$count) ]
        cliquePosition<-(aux_cliques==largestAuxClique)
        PotentialOutClusterReplacement[PotentialOutClusterReplacement]=cliquePosition
    #    print(PotentialOutClusterReplacement[PotentialOutClusterReplacement])
        if (sum(PotentialOutClusterReplacement)>1)
            { #make a swap if there are at least 2 replacements
    #            print("before swap")
    #            print((W[cl_memship==clID,cl_memship==clID]>0)*1)
    #            print("2replace:")
    #            print(potentialInClusterSwap)
    #            print("cl and replacement")
    #            print((W[(cl_memship==clID)|PotentialOutClusterReplacement,(cl_memship==clID)|PotentialOutClusterReplacement]>0)*1)
                cl_memship[PotentialOutClusterReplacement]<-clID
                cl_memship[[potentialInClusterSwap]]<-0
            #update bookkeeping variables and quit for 
            #since operations inside depend on cl_memship 
                Zeros<- cl_memship==0
    #            message(sum(cl_memship==0))
    #            print("after swap")
    #            print((W[cl_memship==clID,cl_memship==clID]>0)*1)
                for_broken=TRUE
                break #quit for
            } 
        } #end for interating on in cluster things
        if (!for_broken) continue_flag=FALSE
            } else continue_flag=FALSE # this happens if no swap was possible
        if (all(!Zeros)) continue_flag=FALSE
    emeregency_stopper=emeregency_stopper+1
    if (emeregency_stopper> 1000)
        stop('reached max iter')
    if (emeregency_stopper%%10==0) message(sum(cl_memship==0))
    }  #end while
    } #end extended mode
    if (all(!Zeros)) break
    Zero2Zero<- W[ Zeros, Zeros ] # edges joining nodes not assigned anywhere yet
  }
  if (return_list) return(cl_list)
  return(cl_memship)
}


isClique<-function(WorA,subgr=NULL){
    if (!is.null(subgr))  WorA<- WorA[subgr,subgr,
                                        drop=FALSE]
    diag(WorA)<-0
    sum((WorA==0))==nrow(WorA)
}

areCliques<-function(WorA,cl_mem){
non0<- unique(cl_mem[cl_mem!=0])
all(
unlapply(non0, function(lbl) isClique(WorA, cl_mem==lbl))
    )
}

#to be applied after greedyCliquesWeighted,
#if one wants the singletons to be assigned unique labels

uniqueSingletonLabels<- function(cl_mem){
  
  singletonMask<- cl_mem==0
  n_clusters<- max(unique(cl_mem))
  cl_mem[singletonMask] = (n_clusters+1):(sum(singletonMask)+n_clusters)
  return(cl_mem)
}

# puts 0 at clusters/ cliques smaller than mcl
zeroOutSmall<- function(labelV, mcl) {
rment<- labelV
non0labs<- unique(labelV[labelV!=0])
non0cnts<- unlapply(non0labs, function(lab) sum(labelV==lab))
for (l in seq_along(non0labs)) 
    if (non0cnts[[l]] < mcl )
       rment[ labelV==non0labs[[l]] ]=0
return(rment)
}

# change non-zero labels to 1:N where N is number of non-zero labels (uniq)

tidyUpLabels<- function(labelV) {
uqn0<-unique(labelV[labelV!=0])
rmentL<- seq_along(uqn0)
rment<- rep(0, length(labelV))
for (l in seq_along(uqn0))
   rment[ labelV == uqn0[[l]] ]= rmentL[[l]]
return(rment)
}

#clique similarity weighted adjacency

cliqueSimilarity<- function(cl_mem, A){
  stopifnot(all(diag(A)==0))
  stopifnot(length(unique(cl_mem))<=ncol(A))
  stopifnot(all(A %in% c(0,1)))
  #all the tests below should pass
  #if they do not, use tidyUpLabels
  cl_mem_uq<-1:max(cl_mem)
  non0<- unique(cl_mem[cl_mem!=0])
  stopifnot(all(cl_mem_uq %in% non0))
  stopifnot(all(non0 %in% cl_mem_uq))
  #to optimize, cbind is slow
  componentIndicator<-do.call(cbind,
                              lapply(cl_mem_uq, function(k) cl_mem==k)) #nxk
  n_con<-A %*% componentIndicator  # n_con[i,k] number of connections of node i to cluster k
  numer<- t(componentIndicator) %*% n_con   # number of connections from cluster ki to cluster kj
  compSizes<- colSums(componentIndicator)
  denom<- compSizes %o% compSizes #for 100% correctness we would have to put..
  #... diag(denom)[[i]]<- compSizes[[i]]*(compSizes[[i]]-1)/2
  numer/denom
}




# joinCliques(threshold,minCoreSize, Similarity)
# 1. take the largest clique, say it is i, if size(i)>= minCoreSize then use it as core, mark it  as checked
# 2. for each of the non checked cliques j:
#    a. if Similarity(i,j)> threshold:
#        attach j to the neighbourhood if i
#        mark i as already checked
#3. if there is some clique k that is not checked yet, take largest, if size(k)>= minCoreSize then use it as core, 
#   mark as checked
#.4. Do the step 2 for the new core.
#5. Go back to 3 if there is still some non checked clique
#6. We are done.
##INITIALIZATION OF K-MEANS LIKE ALGO

joinCliques<- function(SimMat, cl_sizes=NULL, cl_mem=NULL,
                       threshold=0.5,
                       minCoreSize=3 
                       )
{
  if (is.null(cl_sizes) && is.null(cl_mem))
        stop(paste0("give either sizes of the cliques,\n",
                    "cl_sizes[[i]] being the size of clique with label i,\n",
                    "or clique membership labels with cl_mem"))
  if(is.null(cl_sizes)){  #label correctness tests & cl_sizes calc
    non0<-unique(cl_mem[cl_mem!=0])
    cl_mem_uq<-1:max(cl_mem)
    stopifnot(all(cl_mem_uq %in% unique(cl_mem)))
    stopifnot(all(unique(cl_mem) %in% cl_mem_uq))
    cl_sizes<-unlapply(cl_mem_uq, function(k) sum(cl_mem==k))
  }

  stopifnot(ncol(SimMat)==length(cl_sizes))
  S_th <- (SimMat > threshold)*1.
  diag(S_th)<-0
  if (all(lower.tri(S_th)==0)) warning("all values in SimMat are lower than threshold")
  if (all(lower.tri(S_th)==1)) warning("all values in SimMat are higher than threshold")

  bigMask<- cl_sizes >= minCoreSize 
  nodes_checked<- rep(FALSE, ncol(S_th))
  coreList<-vector(mode="list")
  potCoreMask<- (!nodes_checked)&bigMask 
  while (any(potCoreMask))
  {
    maxVal<-max(cl_sizes[potCoreMask])
    core<-which((cl_sizes==maxVal)&potCoreMask)[[1]]
    #message(core)
    nodes_checked[[core]]=TRUE
    leaves<-which((S_th[core,]==1) & !(nodes_checked))
    coreList[[length(coreList)+1]]<-list(
      core=core,
      leaves=leaves)
    nodes_checked[leaves]=TRUE
    potCoreMask<- (!nodes_checked)&bigMask
  }
  return(coreList)
}

#take coreList and number of cliques, return matrix L of FALSEs and TRUEs:
#L is n x C where C is nubmer of cores, n is number of cliques
#L[l,c]=T if clique l is a leaf of core c, FALSE otherwise
coreList2matrix<- function(coreList, n){
 C<- length(coreList)
 coress<- unlapply(coreList, function(x) x$core)
 leavess<-unlapply(coreList, function(x) x$leaves)
 stopifnot(length(unique(coress))==length(coress))
 stopifnot(all(!(coress %in% leavess)))
 L<-matrix(rep(FALSE, n*C), ncol=C)
 for (co in seq_along(coreList))
 {
    if (length(coreList[[co]]$leaves))
    L[ coreList[[co]]$leaves, co ] <- TRUE
 }
 attr(L, 'coreIDs')<- coress
 return(L)
}


#### Reassign leaves step:
### for leaf in all leaves:
###     for core in all cores excluding core(leaf):
###            let core_m = arg max_{cores} S_{core, leaf}
###            if S_{core_m, leaf} > S_{core(leaf), leaf):
###               let core_new(leaf)=core

reassignLeaves<- function(SimMat,threshold,L)
{
    stopifnot(ncol(SimMat)>=ncol(L))
    stopifnot(all(rowSums(L) <2))
    S_th <-SimMat
    S_th[(SimMat <= threshold)]=0
    diag(S_th)<-0
    coreIDs<- attr(L, 'coreIDs')
    stopifnot(length(unique(coreIDs))==length(coreIDs))
    LeafMask<-rowSums(L) > 0 
    stopifnot(all(!(coreIDs %in% which(LeafMask))))
    new_L<- matrix(rep(FALSE, length(L)),ncol=ncol(L))
    if (any(LeafMask)){
    #ATTENTION. this function max.col
    #IS ACTUALLY NONDETERMINISTIC!
    #if we use seed this will run deterministically
    #but not otherwise...
    #we can also add ties.method="first"
    maxCorePerLeaf<- max.col(S_th[LeafMask ,coreIDs, drop=FALSE])
    for (k in seq_along(coreIDs))
    {
        LeafMaskReplacement<- new_L[LeafMask,,drop=FALSE]
        LeafMaskReplacement[maxCorePerLeaf==k,k]<-TRUE
        new_L[LeafMask,]<-LeafMaskReplacement
    }  #@UP 1. rows Leaves & which have k as core
                       }   
    attr(new_L, 'coreIDs')<- coreIDs
    return(new_L)
}


#### Refine cores step:
### for clique cluster in all cliq-cls:
###     for leaf of current core:
####            check if all leaves and core are connected to leaf
####            calculate its degree (0 if not connected)
####    let core_m = arg_max_{leaves} deg_{leaf}
refineCores<- function(SimMat, L,threshold,SimMatColSums=NULL){
    stopifnot(all(rowSums(L) <2))
    new_L<- L
    A_th<- (SimMat > threshold)*1.
    if (is.null(SimMatColSums)) SimMatColSums<- colSums(SimMat)
    allCliqs<-c(1:length(SimMatColSums))
    coreIDs<- attr(L, 'coreIDs')
    stopifnot(length(unique(coreIDs))==length(coreIDs))
    LeafMask<-rowSums(L) > 0 
    stopifnot(all(!(coreIDs %in% which(LeafMask))))
    new_coreIDs<- coreIDs
    for (k in seq_along(coreIDs))
    {
        L[,k]-> kth_coreLeafMask
        kth_coreLeafMask[[coreIDs[[k]]]]=TRUE
        allConnectedLeaves<-(rowSums(A_th[kth_coreLeafMask,
                                kth_coreLeafMask,drop=FALSE])==(sum(kth_coreLeafMask)-1))
        if (any(allConnectedLeaves))
        {
        positionOfMaxInAllConnected<-which.max(SimMatColSums[kth_coreLeafMask][allConnectedLeaves])
        new_coreIDs[[k]] <-  allCliqs[kth_coreLeafMask][allConnectedLeaves][[positionOfMaxInAllConnected]]
        if (new_coreIDs[[k]]!=coreIDs[[k]]){
        new_L[new_coreIDs[[k]],k]<-FALSE
        new_L[coreIDs[[k]],k]<-TRUE}
        }
    }
    stopifnot(length(unique(new_coreIDs))==length(new_coreIDs))
    LeafMask<-rowSums(new_L) > 0 
    stopifnot(all(!(new_coreIDs %in% which(LeafMask))))
    attr(new_L,'coreIDs')<- new_coreIDs
    return(new_L)
}


RunCliqueJoin<- function(SimMat, cl_sizes=NULL, cl_mem=NULL,
                       threshold=0.5,
                       minCoreSize=3,maxIter=50 ,
                       monitorStateChangeN=FALSE
                       )
{
coreList<- joinCliques(SimMat, cl_sizes, cl_mem, threshold, minCoreSize)
if (!length(coreList)) {message('no cores available, i stopped'); return(NULL)}
L<- coreList2matrix(coreList, ncol(SimMat))
SimMatColSums<- colSums(SimMat)
i=0
stateChange=TRUE
if (monitorStateChangeN) N_CHANGES<-vector()

while ( (i<maxIter) &(stateChange))
{i<- i+ 1
    L_1<-reassignLeaves(SimMat,threshold,L)
    if (i==2) fL<-L_1
    L_2<-refineCores(SimMat,L_1,threshold, SimMatColSums)
    n_changesIter<- sum(L_2!=L)
    if (n_changesIter==0) stateChange=FALSE
    if (monitorStateChangeN) N_CHANGES[[i]]<- n_changesIter
    L  <- L_2
}
#if (i>= maxIter) message("reached max iterations.")
if (!monitorStateChangeN) return(L) else {
 return( list(n_changes=N_CHANGES,fL=fL, L=L))
                }
}

coreMatrix2coreList<-function(L)
{
    stopifnot(!is.null(attr(L,'coreIDs')))
    stopifnot(length(attr(L,'coreIDs'))==ncol(L))
    coreIDs<-attr(L,'coreIDs')
    coreList<-vector(mode="list")
    for (k in seq_along(coreIDs))
    {
        coreList[[length(coreList)+1]] = list(core=coreIDs[[k]],
                                              leaves=which(L[,k]))
    }
    return(coreList)
}

conDensity<- function(W){
  W<- W [ lower.tri(W) ]
 sum(W)/ length(W)
}

clusterGrowing<- function(cl_mem, SimMat,origW,threshold=0.5){
    non0<-unique(cl_mem[cl_mem!=0])
    cl_mem_uq<-1:max(cl_mem)
    stopifnot(all(cl_mem_uq %in% unique(cl_mem)))
    stopifnot(all(unique(cl_mem) %in% cl_mem_uq))
    cl_dens<- unlapply(cl_mem_uq, function(lab) conDensity(origW[ cl_mem==lab, cl_mem==lab ])
                      )
    dense2looseOrder<-max2minOrder(cl_dens)
    dense2looseLabs<- cl_mem_uq[dense2looseOrder]
    dense2looseDens<- cl_dens[dense2looseOrder]
    B<- (SimMat[dense2looseOrder,dense2looseOrder]>threshold)*.1
    checked<- rep(FALSE, ncol(B))
    clusters<- rep(-1, ncol(B))
    B_deg<- colSums(B)
    checked[B_deg==0]<-TRUE
    clusters[B_deg==0]<-0
    C=0
    while (any(!checked))
    {
        C=C+1
        clique_nr<- which(!checked)[[1]]
        checked[[clique_nr]]<-TRUE
        clusters[[clique_nr]]<-C
        current_cluster<- clusters==C
        neigh<- (colSums(B[current_cluster,,drop=FALSE])>0)
        neigh[current_cluster]<-FALSE
        while (any(neigh)){
        checked[neigh]<-TRUE
        clusters[neigh]<-C
        current_cluster<- clusters==C
        neigh<- (colSums(B[current_cluster,,drop=FALSE])>0)
        neigh[current_cluster]<-FALSE
                        }
    }
    rbind(dense2looseLabs,
          clusters)
}

cliqueChaining<- function(cl_mem, SimMat,origW,threshold=0.5){
    non0<-unique(cl_mem[cl_mem!=0])
    cl_mem_uq<-1:max(cl_mem)
    stopifnot(all(cl_mem_uq %in% unique(cl_mem)))
    stopifnot(all(unique(cl_mem) %in% cl_mem_uq))
    cl_dens<- unlapply(cl_mem_uq, function(lab) conDensity(origW[ cl_mem==lab, cl_mem==lab ])
                      )
    dense2looseOrder<-max2minOrder(cl_dens)
    dense2looseLabs<- cl_mem_uq[dense2looseOrder]
    dense2looseDens<- cl_dens[dense2looseOrder]
    B<- (SimMat[dense2looseOrder,dense2looseOrder]>threshold)*.1
    checked<- rep(FALSE, ncol(B))
    clusters<- rep(-1, ncol(B))
    B_deg<- colSums(B)
    checked[B_deg==0]<-TRUE
    clusters[B_deg==0]<-0
    C=0
    steps<-list()
    while (any(!checked))
    {
        C=C+1
        clique_nr<- which(!checked)[[1]]
        checked[[clique_nr]]<-TRUE
        clusters[[clique_nr]]<-C
        current_cluster<- clusters==C
        neigh<- (colSums(B[current_cluster,,drop=FALSE])>0)
        neigh[current_cluster]<-FALSE
        while (any(neigh)){
        checked[neigh]<-TRUE
        clusters[neigh]<-C
        current_cluster<- clusters==C
        neigh<- (colSums(B[current_cluster,,drop=FALSE])>0)
        neigh[current_cluster]<-FALSE
        steps<- c(steps, 
                  list(
                    rbind(dense2looseLabs,
                          clusters)
                  )
                  )
                        }
    }
    list(compar=
    rbind(dense2looseLabs,
          clusters),
         steps=steps)
}
cliqueChainingPreferential<- function(cl_mem, SimMat,origW,threshold=0.5){
max_slack=0
    non0<-unique(cl_mem[cl_mem!=0])
    cl_mem_uq<-1:max(cl_mem)
    stopifnot(all(cl_mem_uq %in% unique(cl_mem)))
    stopifnot(all(unique(cl_mem) %in% cl_mem_uq))
    cl_dens<- unlapply(cl_mem_uq, function(lab) conDensity(origW[ cl_mem==lab, cl_mem==lab ])
                      )
    dense2looseOrder<-max2minOrder(cl_dens)
    dense2looseLabs<- cl_mem_uq[dense2looseOrder]
    dense2looseDens<- cl_dens[dense2looseOrder]
    S<- SimMat[dense2looseOrder, dense2looseOrder]
    B<- (S>threshold)*.1
    checked<- rep(FALSE, ncol(B))
    clusters<- rep(-1, ncol(B))
    B_deg<- colSums(B)
    checked[B_deg==0]<-TRUE
    clusters[B_deg==0]<-0
    C=0
    steps<-list()
    while (any(!checked))
    {
        C=C+1
        clique_nr<- which(!checked)[[1]]
        checked[[clique_nr]]<-TRUE
        clusters[[clique_nr]]<-C
        current_cluster<- clusters==C
        neigh<- (colSums(B[current_cluster,,drop=FALSE])>0)
        neigh[current_cluster]<-FALSE
        while (any(neigh)){
        neighS<- max(apply(S[current_cluster,neigh,drop=FALSE],2,max))
        s_neigh<-which.max(neighS) 
        added<-which(neigh)[s_neigh]
        checked[added]<-TRUE
        clusters[added]<-C
        current_cluster<- clusters==C
        neigh<- (colSums(B[current_cluster,,drop=FALSE])>0)
        neigh[current_cluster]<-FALSE
        n_try=0
        sl_tr<-threshold
        while ((!any(neigh))&n_try<max_slack){
         sl_tr<- sl_tr-0.01
         neigh<- (colSums(S[current_cluster,,drop=FALSE]>sl_tr)>0)
         n_try<-n_try+1
        neigh[current_cluster]<-FALSE}
        steps<- c(steps, 
                  list(
                    rbind(dense2looseLabs,
                          clusters)
                  )
                  )
                        }
    }
    list(compar=
    rbind(dense2looseLabs,
          clusters),
         steps=steps)
}


cliques2clusters<- function(cliqueMembership,
                            cliqueClustering)
{
    if (is.matrix(cliqueClustering)) cliqueClustering<- coreMatrix2coreList(cliqueClustering)
    stopifnot(length(unlist(cliqueClustering))<=length(unique(cliqueMembership)))
    the_leaves<- unlist(lapply(cliqueClustering, function(x) x$leaves))
    the_cores<- unlist(lapply(cliqueClustering, function(x) rep(x[['core']],
                                                                length(x$leaves))                 ))

    coreMembership<-unlist(lapply(seq_along(cliqueMembership),
            function(i)
            if (cliqueMembership[[i]] %in% the_leaves) {#if the clique of ith node is a leaf, 
                                                        #return the label of the clique that is the core to which leaf belongs.
                the_cores[[which(the_leaves==cliqueMembership[[i]])]]
            }else{                                      #otherwise, the clique of i-th node is a singleton OR a core, so just return the clique label
                cliqueMembership[[i]]}
                ))
    return(coreMembership) 
}





betweenDensity<- function(cl_mem,W){
 free<- cl_mem==0
 freeW<- W[free,free]
 freeW<- freeW[ lower.tri(freeW) ]
 not_free<- !free
 allW<- W[!free,!free]
 allW<- allW[ lower.tri(allW) ]
 non0labs<- unique(cl_mem[!free])
 clW<- unlapply( non0labs, function(lab) { thisW<- W[cl_mem==lab, 
                                                     cl_mem==lab]
                                            thisW[lower.tri(thisW)]
                                         }
                )
 (sum(allW)-sum(clW)) /((length(allW)- length(clW)))
}

inClDensity<- function(cl_mem,W){
 free<- cl_mem==0
 non0labs<- unique(cl_mem[!free])
 clW<- unlapply( non0labs, function(lab) { thisW<- W[cl_mem==lab, 
                                                     cl_mem==lab]
                                            thisW[lower.tri(thisW)]
                                         }
                )
 (sum(clW)) /(length(clW))
}

clCoverage<- function(cl_mem, W){
 allW<- W[lower.tri(W)]
 free<- cl_mem==0
 non0labs<- unique(cl_mem[!free])
 clW<- unlapply( non0labs, function(lab) { thisW<- W[cl_mem==lab, 
                                                     cl_mem==lab]
                                            thisW[lower.tri(thisW)]
                                         }
                )
 (sum(clW)) /(sum(allW))
}

extendCliques<-function(old_labels, newW, #newcomer should be last
              mode='average')
              {
               stopifnot(all(diag(newW)==0))
               new_labels<- c(old_labels, max(old_labels)+1)
                newc<- ncol(newW)
                possible= newW[newc,][-newc]>0
                if (any(possible)){
                old_labelsUq<- unique(old_labels)
                potClqsStr<-unlapply(old_labelsUq,
                function(ol) 
               (all(newW[newc, new_labels==ol]>0))*mean(newW[newc,new_labels==ol])
                )
                if (any(potClqsStr>0))
                    {
                    newcomer_l<- old_labelsUq[which.max(potClqsStr)]
                    new_labels[newc]<- newcomer_l
                    }
                } 
                return(new_labels)
              }

extendCliqueJoin<-function(L_old, new_cl_mem,
                 newA, sc_thr=0.5,
                 maxA=1){
   #stopifnot(all(newA %in% c(0,1)) )
   oldCores<-attr(L_old, 'coreIDs')
   new_cl_mem_uq<-unique(new_cl_mem) 
   L_new<- rep(FALSE, ncol(L_old)*length(new_cl_mem_uq))
   L_new<- matrix(L_new, ncol=ncol(L_old),
                         nrow=length(new_cl_mem_uq))
   L_new[1:nrow(L_old),]<-L_old
   attr(L_new, 'coreIDs')<- oldCores
   oldClq<-unlist(coreMatrix2coreList(L_old))
   joined2old<- new_cl_mem_uq %in% oldClq
   new_cl_mem_uq<- new_cl_mem_uq[!joined2old]
   print(new_cl_mem_uq)
   for (new_mem_lab in new_cl_mem_uq)
   {    oldCoreSims<-unlapply(oldCores, function(old_lab)
        {
            oldCorePos<-new_cl_mem==old_lab
            newClqPos<- new_cl_mem== new_mem_lab
            N_con<-sum(newA[oldCorePos,newClqPos])
            N_possible<- maxA*sum(oldCorePos)*sum(newClqPos)
            SimScore<-N_con/N_possible
            SimScore*(SimScore>sc_thr)
        })
        if (max(oldCoreSims)>0){
        newMemCore<- oldCores[ which.max(oldCoreSims)]
        print((newMemCore))
        L_new[new_mem_lab ,oldCores==newMemCore]<-TRUE
        }
   }
    return(L_new)
                 }

cliqueDensitySimilarity<- function(cl_mem, W){
  stopifnot(all(diag(W)==0))
  stopifnot(length(unique(cl_mem))<=ncol(W))
  cl_mem_uq<-1:max(cl_mem)
  #to optimize, cbind is slow
  componentIndicator<-do.call(cbind,
                              lapply(cl_mem_uq, function(k) cl_mem==k)) #nxk
  compSizes<- colSums(componentIndicator) #1xk
  wSums<- t(componentIndicator) %*% (W %*% componentIndicator) #kxk
           # sum of the weights inside cliques (on the diagonal)
           # and between cliques (off-diagonal)
  diag(wSums)<- diag(wSums)/2 #diagonal was counted twice
  maxEdges<- compSizes %o% compSizes  #kxk, maximum possible number of edges
  diag(maxEdges)<- (compSizes *(compSizes-1))/2# diagonal was counted twice
  #S_ij:= U_dens_ij / max(U_dens_ij, I_dens_ij)
  #where U_dens_ij is density of union of clique i and j
  #I_dens_ij is mean of densities of clique i and j
  #density of clique k is (sum of weights in k)/(number of possible edges)
  I_maxEdg<- diag(maxEdges)
  I_wSum  <- diag(wSums)
  I_maxEdg[I_maxEdg==0]<-1 #to avoid NAs below
  I_dens<- outer(I_wSum/I_maxEdg,
                 I_wSum/I_maxEdg, 
                 function(x,y) (x+y)/2) #kxk
  ijs<-as.list(as.data.frame(t(expand.grid(1:ncol(wSums),1:ncol(wSums)))))
  U_wSum<- matrix(ncol=ncol(wSums), nrow=nrow(wSums))
  U_maxEdg<-matrix(ncol=ncol(wSums), nrow=nrow(wSums))
  for (ij in ijs) U_wSum[ij[[1]], ij[[2]] ] <- I_wSum[[ij[[1]]]] + I_wSum[[ij[[2]]]] + wSums[ ij[[1]], ij[[2]] ]
  for (ij in ijs) U_maxEdg[ij[[1]], ij[[2]] ] <- diag(maxEdges)[[ij[[1]]]] + diag(maxEdges)[[ij[[2]]]] + maxEdges[ ij[[1]], ij[[2]] ]
  U_maxEdg[ U_maxEdg==0 ] <- 1
  numer<- U_wSum/U_maxEdg 
  diag(maxEdges)<-I_maxEdg
  denom<- pmax( I_dens, U_wSum) + I_dens
  RES<-numer/denom
  RES<-wSums/maxEdges
  diag(RES)<-1
  RES
}
