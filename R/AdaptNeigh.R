"AdaptNeigh" <-
function(pointsin,X,coeff,nbrs,remove,intercept,neighbours){

#does local adaptive prediction for the point remove based on N 
#points (chooses method of prediction and intercept);

#library(Matrix);

mindetails<-NULL
minindices<-NULL
results<-list()
tempres<-list()
newinfo<-list()
nlist<-list()

N<-length(pointsin);
r<-which(pointsin==remove);

min1<-min(N-1,neighbours)
min2<-min(N-1,2*neighbours)

closest<-FALSE   #does prediction schemes with neighbours on each side

for (k in 1:min1){

out1<-getnbrs(X,remove,pointsin,k,closest)
nbrs<-out1$nbrs
index<-out1$index

out2<-AdaptPred(pointsin,X,coeff,nbrs,remove,intercept,neighbours)
nlist[[k]]<-out1
tempres[[k]]<-out2

}


closest<-TRUE

for (k in 1:min2){

out1<-getnbrs(X,remove,pointsin,k,closest)
nbrs<-out1$nbrs
index<-out1$index

out2<-AdaptPred(pointsin,X,coeff,nbrs,remove,intercept,neighbours)
nlist[[k+min1]]<-out1
tempres[[k+min1]]<-out2
}


for (i in 1:(min1+min2)){
minindices[i]<-tempres[[i]][[10]]
mindetails[i]<-tempres[[i]][[9]][minindices[i]]
}

totalminindex<-order(abs(mindetails))[1]
pred<-coeff[remove]-mindetails[totalminindex]
coeff[remove]<-mindetails[totalminindex]
clo<-NULL
int<-NULL
scheme<-NULL
if(totalminindex<=neighbours){clo<-FALSE}else{clo<-TRUE}

results<-tempres[[totalminindex]]
nbrs<-nlist[[totalminindex]]$nbrs
index<-nlist[[totalminindex]]$index

newinfo<-list(clo=clo,totalminindex=totalminindex,nbrs=nbrs,index=index,mindetails=mindetails,minindices=minindices)


return(list(results=results,newinfo=newinfo))
}