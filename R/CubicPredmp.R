"CubicPredmp" <-
function(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g){

#does local cubic  prediction for the point remove based on N 
#points (with intercept as default);

#library(Matrix);

if (length(pointsin)==2){
#if (all(nbrs==pointsin)){
Xneighbours<-cbind(X[nbrs],X[nbrs]^2,X[nbrs]^3);
}
else{
Xneighbours<-cbind(X[newnbrs],X[newnbrs]^2,X[newnbrs]^3)
}
Xremove<-cbind(X[remove],X[remove]^2,X[remove]^3);

if (intercept){
	Xneighbours<-cbind(1,Xneighbours);
	Xremove<-as.row(c(1,Xremove));
		}

if (length(nbrs)>=4){       #possible to have curve with intercept
	
#	temp<-as.matrix(t(Xneighbours)%*%Xneighbours);	
#	mm<-solve.Matrix(temp,t(Xneighbours));

temp<-crossprod(Xneighbours)
mm<-Rmatsolve(temp)%*%t(Xneighbours)

#	mm<-solve(t(Xneighbours)%*%Xneighbours,t(Xneighbours));

coeff1<-NULL
#	if (all(nbrs==pointsin)){
for (i in 1:length(nbrs)){
	coeff1<-cbind(coeff1,as.row(coefflist[[nbrs[i]]]))
	}
#}
#else{
#	for (i in 1:length(nbrs)){
#	coeff1<-cbind(coeff1,coefflist[[nbrs[i]]])
#	}
#}
	
	bhat<-mm%*%matrix(coeff1,ncol=1);
	
	pred<-Xremove%*%bhat;
		
	weights<-matrix(Xremove,nrow=1)%*%mm;  #works out neighbour weights
	
		}
	
if(length(nbrs)==3){
#if we don't have enough (4) neighbours to define a cubic 
#automatically, we do quadratic prediction		

		
	Xneighbours<-Xneighbours[,1:(2+intercept)]
	Xremove<-as.row(Xremove[,1:(2+intercept)])	


#this does it^^^

#	temp<-as.matrix(t(Xneighbours)%*%Xneighbours);	
#	mm<-solve.Matrix(temp,t(Xneighbours));

temp<-crossprod(Xneighbours)
mm<-Rmatsolve(temp)%*%t(Xneighbours)

#	mm<-solve(t(Xneighbours)%*%Xneighbours,t(Xneighbours));

	coeff1<-NULL
#	if (all(nbrs==pointsin)){
for (i in 1:length(nbrs)){
	coeff1<-cbind(coeff1,as.row(coefflist[[nbrs[i]]]))
	}
#}
#else{
#	for (i in 1:length(nbrs)){
#	coeff1<-cbind(coeff1,coefflist[[nbrs[i]]])
#	}
#}

	bhat<-mm%*%matrix(coeff1,ncol=1);
	
	pred<-Xremove%*%bhat;
		
	weights<-matrix(Xremove,nrow=1)%*%mm;  #works out neighbour weights
	
		}
if(length(nbrs)==2){
	
#we dont have enough neighbours for a cubic, so we do linear
#prediction...
		 
	Xneighbours<-Xneighbours[,1:(1+intercept)]
	Xremove<-as.row(Xremove[,1:(1+intercept)])	

#this does it^^^

#	temp<-as.matrix(t(Xneighbours)%*%Xneighbours);	
#	mm<-solve.Matrix(temp,t(Xneighbours));

temp<-crossprod(Xneighbours)
mm<-Rmatsolve(temp)%*%t(Xneighbours)

#	mm<-solve(t(Xneighbours)%*%Xneighbours,t(Xneighbours));
	coeff1<-NULL
#	if (all(nbrs==pointsin)){
for (i in 1:length(nbrs)){
	coeff1<-cbind(coeff1,as.row(coefflist[[nbrs[i]]]))
	}
#}
#else{
#	for (i in 1:length(nbrs)){
#	coeff1<-cbind(coeff1,coefflist[[nbrs[i]]])
#	}
#}
	bhat<-mm%*%matrix(coeff1,ncol=1);

	pred<-Xremove%*%bhat;
		
	weights<-matrix(Xremove,nrow=1)%*%mm;  #works out neighbour weights

	}
		
if(length(nbrs)==1){
    mm<-0;
    bhat<-1;
    weights<-1;
    pred<-coefflist[[nbrs]];
    };

#coeff[remove]<-coeff[remove]-pred; #works out the detail coefficient for
				   # removed point


coeff<-as.row(coeff)
#print("done Cubic")

return(list(Xneighbours=Xneighbours,mm=mm,bhat=bhat,weights=weights,pred=pred,coeff=coeff));

}
