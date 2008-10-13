"LinearPred" <-
function(pointsin,X,coeff,nbrs,remove,intercept,neighbours){

#does local linear prediction for the point remove based on N points (with 
#intercept as default);

#library(Matrix);

Xneighbours<-X[nbrs];
Xneighbours<-as.column(Xneighbours);
Xremove<-X[remove];

if (intercept){
	Xneighbours<-cbind(1,Xneighbours);
	Xremove<-as.row(c(1,Xremove));
		}

if (length(nbrs)>=2){
#	temp<-Matrix(t(Xneighbours)%*%Xneighbours);	
#	mm<-solve(temp,t(Xneighbours));

temp<-crossprod(Xneighbours)
mm<-Rmatsolve(temp)%*%t(Xneighbours)


#	mm<-solve(t(Xneighbours)%*%Xneighbours,t(Xneighbours));
	bhat<-mm%*%matrix(coeff[nbrs],ncol=1);
	
	pred<-Xremove%*%bhat;

	
	weights<-matrix(Xremove,nrow=1)%*%mm;  #works out neighbour weights

#	pred<-sum(as.column(weights)*coeff[nbrs])
	#weights1[1]<-(X[nbrs[2]]-X[remove])/2;
	#weights1[2]<-(X[remove]-X[nbrs[1]])/2;
	
		}
else{
    mm<-0;
    bhat<-1;
    weights<-1;
    pred<-coeff[nbrs];
    };

#coeff[remove]<-coeff[remove]-pred; #works out the detail coefficient for
				   # removed point


return(list(Xneighbours=Xneighbours,mm=mm,bhat=bhat,weights=weights,pred=pred,coeff=coeff));

}
