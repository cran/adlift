"CubicPred" <-
function(pointsin,X,coeff,nbrs,remove,intercept,neighbours){

#does local cubic  prediction for the point remove based on N 
#points (with intercept as default);

#library(Matrix);

Xneighbours<-cbind(X[nbrs],X[nbrs]^2,X[nbrs]^3);

Xremove<-cbind(X[remove],X[remove]^2,X[remove]^3);

if (intercept){
	Xneighbours<-cbind(1,Xneighbours);
	Xremove<-as.row(c(1,Xremove));
		}

if (length(nbrs)>=4){       #possible to have curve with intercept
	
#	temp<-Matrix(t(Xneighbours)%*%Xneighbours);	
#	mm<-solve(temp,t(Xneighbours));

temp<-crossprod(Xneighbours)
mm<-Rmatsolve(temp)%*%t(Xneighbours)

#	mm<-solve(t(Xneighbours)%*%Xneighbours,t(Xneighbours));
	bhat<-mm%*%matrix(coeff[nbrs],ncol=1);
	
	pred<-Xremove%*%bhat;
		
	weights<-matrix(Xremove,nrow=1)%*%mm;  #works out neighbour weights
	#weights1[1]<-(X[nbrs[2]]-X[remove])/2;
	#weights1[2]<-(X[remove]-X[nbrs[1]])/2;
	
		}
	
if(length(nbrs)==3){
#if we don't have enough (4) neighbours to define a cubic 
#automatically, we do quadratic prediction		

		
	Xneighbours<-Xneighbours[,1:(2+intercept)]
	Xremove<-as.row(Xremove[,1:(2+intercept)])	
	#print(Xremove)
	#print(Xneighbours)

#this does it^^^

#	temp<-Matrix(t(Xneighbours)%*%Xneighbours);	
#	mm<-solve(temp,t(Xneighbours));

temp<-crossprod(Xneighbours)
mm<-Rmatsolve(temp)%*%t(Xneighbours)


#	mm<-solve(t(Xneighbours)%*%Xneighbours,t(Xneighbours));
	bhat<-mm%*%matrix(coeff[nbrs],ncol=1);
	
	pred<-Xremove%*%bhat;
		
	weights<-matrix(Xremove,nrow=1)%*%mm;  #works out neighbour weights
	
		}
if(length(nbrs)==2){
	
#we dont have enough neighbours for a cubic, so we do linear
#prediction...
		 
	Xneighbours<-Xneighbours[,1:(1+intercept)]
	Xremove<-as.row(Xremove[,1:(1+intercept)])	

#this does it^^^

#	temp<-Matrix(t(Xneighbours)%*%Xneighbours);	
#	mm<-solve(temp,t(Xneighbours));

temp<-crossprod(Xneighbours)
mm<-Rmatsolve(temp)%*%t(Xneighbours)


#	mm<-solve(t(Xneighbours)%*%Xneighbours,t(Xneighbours));

	bhat<-mm%*%matrix(coeff[nbrs],ncol=1);

	pred<-Xremove%*%bhat;
		
	weights<-matrix(Xremove,nrow=1)%*%mm;  #works out neighbour weights
	}
		
if(length(nbrs)==1){
    mm<-0;
    bhat<-1;
    weights<-1;
    pred<-coeff[nbrs];
    };

#coeff[remove]<-coeff[remove]-pred; #works out the detail coefficient for
				   # removed point


return(list(Xneighbours=Xneighbours,mm=mm,bhat=bhat,weights=weights,pred=pred,coeff=coeff));

}

