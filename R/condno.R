"condno" <-
function(W, type){

type<-type[1]	#only takes first letter as type

#	library(Matrix)
	W <- as.matrix(W)
	Winv <- Rmatsolve(W)
#	Winv <- solve.Matrix(W)
	if (type=="l1"){
		n1<-sum(abs(W))
		n2<-sum(abs(Winv))
		}
#else{
#	n1 <- norm.Matrix(W, type = type)
#	n2 <- norm.Matrix(Winv, type = type)
#		}

if(type=="1" | type=="o" | type=="O"){
	n1<-max(abs(apply(W,2,sum)))
	n2<-max(abs(apply(Winv,2,sum)))
		}
if(type=="i" | type=="I"){
	n1<-max(abs(apply(W,1,sum)))
	n2<-max(abs(apply(Winv,1,sum)))
		}
if(type=="f" | type=="F"){
	q<-as.vector(W)
	r<-as.vector(Winv)
	n1<-sum(abs(q)^2)
	n2<-sum(abs(r)^2)
		}
if(type=="m" | type=="M"){
	n1<-max(abs(W))
	n2<-max(abs(Winv))
		}

	condno <- n1 * n2
	condno
}
