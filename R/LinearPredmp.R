"LinearPredmp" <-
function(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g){

#does local linear prediction for the point remove based on N points (with 
#intercept as default);
# with multiple point consideration

#library(Matrix);

if (length(pointsin)==2){
#if (all(nbrs==pointsin)){
Xneighbours<-X[nbrs];
}
else{
Xneighbours<-X[newnbrs];
}

Xneighbours<-as.column(Xneighbours);
Xremove<-X[remove];

if (intercept){
	Xneighbours<-cbind(1,Xneighbours);
	Xremove<-as.row(c(1,Xremove));
		}

if (length(nbrs)>=2){
#	temp<-as.matrix(t(Xneighbours)%*%Xneighbours);	
#	mm<-solve.Matrix(temp,t(Xneighbours));

temp<-crossprod(Xneighbours)
mm<-Rmatsolve(temp)%*%t(Xneighbours)

#	mm<-solve(t(Xneighbours)%*%Xneighbours,t(Xneighbours));
	coeff1<-NULL


for (i in 1:length(nbrs)){
	coeff1<-cbind(coeff1,as.row(coefflist[[nbrs[i]]]))


	}



	bhat<-mm%*%matrix(coeff1,ncol=1);

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
    pred<-coefflist[[nbrs]];
    };

#coeff[remove]<-coeff[remove]-pred; #works out the detail coefficient for
				   # removed point

coeff<-as.row(coeff)

#print("done linear")

return(list(Xneighbours=Xneighbours,mm=mm,bhat=bhat,weights=weights,pred=pred,coeff=coeff));

}
