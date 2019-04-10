smoothness <- function(V,Tr,d,r){
	fit_tdata <- tdata(V,Tr)
	E <- fit_tdata$E
	TE <- fit_tdata$TE
	TV <- fit_tdata$TV
	EV <- fit_tdata$EV
	B <- fit_tdata$Bdr
	int <- which(colSums(TE)>1)
	n <- length(int)
	N <- 0
	Neq <- 0
	I_fit <- CrArrays(d,r)
	I1 <- I_fit$I1
	I2 <- I_fit$I2
	Is <- list()
	Js <- list()
	Ks <- list()
	for(j in 0:r){
		N <- N+((j+1)*(j+2)/2+1)*(d+1-j)
		Neq <- Neq+d+1-j
		ind_fit <- indices(j)
		LI <- ind_fit$I
		LJ <- ind_fit$J
		LK <- ind_fit$K
		Is[[j+1]] <- LI
		Js[[j+1]] <- LJ
		Ks[[j+1]] <- LK
	}
# Neq = the number of equations generated in each triangle
# N = the number of nonzero coefficients in the equations generated in
# each triangle.

	m <- (d+1)*(d+2)/2
	Index1 <- rep(0,N*n)
	Index2 <- rep(0,N*n)
	Values <- Index1
	A <- rbind(c(1,2),c(2,3),c(3,1),c(2,1),c(3,2),c(1,3))
	for(j in 1:n){
		k <- int[j]
		v1 <- E[k,1]
		v2 <- E[k,2]
		AdjT <- which(TE[,k]!=0)
		t1 <- AdjT[1]
		t2 <- AdjT[2]
		T1 <- Tr[t1,]
		T2 <- Tr[t2,]
		i1 <- which(T1==v1)
		i2 <- which(T1==v2)
		j1 <- which(T2==v1)
		j2 <- which(T2==v2)
		e1 <- which(A[,1]==i1 & A[,2]==i2)
		e2 <- which(A[,1]==j1 & A[,2]==j2)
		if(length(e1)>0 & e1>3){
			e1 <- e1-3
			Temp <- T1
			T1 <- T2
			T2 <- Temp
			Temp <- e1
			e1 <- e2
			e2 <- Temp
			Temp <- t1
			t1 <- t2
			t2 <- Temp
		} else {e2 <- e2-3}
		v4 <- setdiff(T2,c(v1,v2))
		V4 <- V[v4,]
		lam_fit <- bary(V[T1[1],],V[T1[2],],V[T1[3],],V4[1],V4[2])
		lam1 <- lam_fit[[1]]
		lam2 <- lam_fit[[2]]
		lam3 <- lam_fit[[3]]
		lambda <- rbind(lam1,lam2,lam3)
		if(e1==2){
			temp <- lambda[1]
			lambda[1] <- lambda[2]
			lambda[2] <- lambda[3]
			lambda[3] <- temp
		} else if(e1==3){
			temp <- lambda[2]
			lambda[2] <- lambda[3]
			lambda[3] <- temp
		}
		VarCT <- 0
		EqCt <- 0
		for(i in 0:r){
			Lambda <- factorial(i)/gamma(Is[[i+1]]+1)*gamma(Js[[i+1]]+1)*
			  gamma(Ks[[i+1]]+1)*lambda[1]^Is[[i+1]]*
			  lambda[2]^Js[[i+1]]*lambda[3]^Ks[[i+1]]
			T1mat <- I1[[i+1]][[e1]]
			T2vector <- I2[[i+1]][[e2]]
			numeq <- dim(T1mat)[1]
			numVar <- (dim(T1mat)[2]+1)*numeq

			T1Values <- matrix(1,nrow=dim(T1mat)[1],
				ncol=dim(T1mat)[2])%*%diag(c(Lambda))

			eqnums <- (j-1)*Neq+EqCt+diag(1:numeq)%*%
				matrix(1,nr=numeq,nc=dim(T1mat)[2]+1)

			Index1[((j-1)*N+VarCT+1):((j-1)*N+VarCT+numVar)] <- eqnums

			Index2[((j-1)*N+VarCT+1):((j-1)*N+VarCT+numVar)] <- 
			  cbind(((t1-1)*m+T1mat),((t2-1)*m+T2vector))

			Values[((j-1)*N+VarCT+1):((j-1)*N+VarCT+numVar)] <- 
				rbind(matrix(c(T1Values),ncol=1),matrix(-1,
				nrow=dim(T2vector)[1],ncol=dim(T2vector)[2]))

			VarCT <- VarCT+numVar
			EqCt <- EqCt+numeq
		}
	}
  H_sparse <- sparseMatrix(Index1,Index2,x=Values,dims=c(n*Neq,dim(Tr)[1]*m))
  H <- H_sparse
  # H=as.matrix(H_sparse)
	return(H)
}
