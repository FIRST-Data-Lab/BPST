tdata <- function(V,Tr){
  m <- dim(Tr)[1]
  n <- dim(Tr)[2]
  E <- c()
  TE <- matrix(nrow=m,ncol=1)
  TE <- TE[,-1]
  numEdges <- 0
  for (i in 1:m){
    Ti <- Tr[i,]
    for (j in 1:3){
      edge <- c(min(Ti[j],Ti[j%%3+1]), max(Ti[j], Ti[j%%3+1]))
      #edge <- sort(c(Ti[j],Ti[j%%3+1]))
      if (length(E)>0){
        edgenum <- which((edge[1]==E[,1]) & (edge[2]==E[,2]))
      } else edgenum <- c()
      if (length(edgenum)==0){
        E <- rbind(E,edge)
        numEdges <- numEdges+1
        TE <- cbind(TE,sparseMatrix(c(),c(),dims=c(m,1)))
        edgenum <- numEdges
      }
      TE[i,edgenum] <- 1
    }
  }
  TE <- as.matrix(TE)
  numV <- dim(V)[1]
  TV <- Matrix(0,nrow=m,ncol=numV,sparse=T)
  for (i in 1:m){
    TV[i,Tr[i,]] <- matrix(1,ncol=3,nrow=1)
  }
  EV <- Matrix(0,nrow=numEdges,ncol=numV,sparse=T)
  for (i in 1:numEdges){
    EV[i,E[i,]] <- matrix(1,nrow=1,ncol=2)
  }
  Bdr <- findbdt(Tr,V,E,TE,EV)
  list(E=E,TE=TE,TV=TV,EV=EV,Bdr=Bdr)
}
