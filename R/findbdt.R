findbdt <- function(Tr,V,E,TE,EV){
	BE <- which(colSums(TE)==1)
	B <- c()
	# BE = a list of bdr edges, but not yet in cc order
	while(length(BE)>0){
		enum <- BE[1]
		BEE <- matrix(E[BE,],ncol=2)
		e <- E[enum,]
		tt <- which(TE[,enum]!=0)
		Tri <- Tr[tt,]
		V1 <- V[Tri[1],]
		V2 <- V[Tri[2],]
		V3 <- V[Tri[3],]
		P <- (V1+V2+V3)/3
		VL <- V[e[1],]
		VR <- V[e[2],]
		if (triarea(VL,VR,P)>0) {
			B1 <- e
		} else {B1 <- c(e[2],e[1])}
		I1 <- 1
		loop <- 0
		i <- 1
		while(loop==0){
			v1 <- B1[i]
			v2 <- B1[i+1]
		# [v1,v2] is the last bdr edge, with the cc direction from v1 to v2
			next_e <- which(((v2==BEE[,1]) & (v1!=BEE[,2])) | 
			                  ((v2==BEE[,2])&(v1!=BEE[,1])))
			if(identical(BEE[next_e,1],v2)){
				v3 <- BEE[next_e,2]
			} else {v3 <- BEE[next_e,1]}
			if(length(v3)>0 && length(intersect(v3,B1[1]))==0){
				B1 <- c(B1,v3)
				I1 <- c(I1,next_e)
				i <- i+1
			} else {I1 <- c(I1,next_e);loop <- 1}
		}
		BE <- setdiff(BE,BE[I1])
		B1 <- matrix(B1,ncol=1)
		B <- newcol(B,B1)
	}
	return(B)
}
