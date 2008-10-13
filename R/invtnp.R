"invtnp" <-
function(X,coeff,lengths,lengthsremove,pointsin,removelist,neighbrs,schemehist,interhist,nadd=length(X)-2,intercept=TRUE,neighbours=1,closest=FALSE,LocalPred=LinearPred) {

X<-as.row(X);
coeff<-as.row(coeff);
n<-length(X);
N<-length(pointsin);
m<-length(removelist);
d<-neighbours;

nadd<-min(nadd,m);

gamlist<-list()
predlist<-NULL
alphalist<-list()
lengthlist<-list()



for (j in 1:nadd) {

N<-length(pointsin);

remove<-removelist[m-j+1];
lengthrem<-lengthsremove[m-j+1];
nbrs<-neighbrs[[m-j+1]];

index<-NULL;
for (i in 1:length(nbrs)) {
	index[i]<-which(pointsin==nbrs[i]);
			}



B<-(X[remove]>X[nbrs])		#checks where to place removed point
nt<-sum(B)			#by counting number of X positions on its 
				#left
if (nt==0){
r<-which(pointsin==nbrs[1])}
if (nt==length(nbrs)){r<-which(pointsin==nbrs[length(nbrs)])+1;}

if ((nt>0)&(nt<length(nbrs))){			
	r<-which(pointsin==nbrs[nt+1])  
	} 

if (is.null(schemehist)==FALSE){
    if (schemehist[m-j+1]=="Linear"){
       res<-LinearPred(pointsin,X,coeff,nbrs,remove,intercept=interhist[m-j+1],neighbours)}
    if (schemehist[m-j+1]=="Quad"){
       res<-QuadPred(pointsin,X,coeff,nbrs,remove,intercept=interhist[m-j+1],neighbours)}
    if (schemehist[m-j+1]=="Cubic"){
       res<-CubicPred(pointsin,X,coeff,nbrs,remove,intercept=interhist[m-j+1],neighbours)}
			    }
else{ 
res<-LocalPred(pointsin,X,coeff,nbrs,remove,intercept,neighbours)}
if(length(res)==2){
	l<-res[[1]]}
else{l<-res}
gamweights<-l[[4]]
#pred<-l[[5]]
bhat<-l[[3]]

l1<-UndoPointsUpdate(X,coeff,nbrs,index,remove,r,N,pointsin,gamweights,lengths,lengthrem);
coeff<-l1$coeff;
lengths<-l1$lengths;
alpha<-l1$alpha;

pred<-sum(as.column(gamweights)*coeff[nbrs])
coeff[remove]<-coeff[remove]+pred

gamlist[[j]]<-gamweights
predlist[j]<-pred
alphalist[[j]]<-alpha


removelist<-setdiff(removelist,remove);


if (r==1){
	lengths<-c(lengthrem,lengths)
	pointsin<-c(remove,pointsin)
	}
if (r==(N+1)){
	lengths<-c(lengths,lengthrem)
	pointsin<-c(pointsin,remove)
	}
if ((r>1)&(r<(N+1))){
	lengths<-c(lengths[1:(r-1)],lengthrem,lengths[r:N])
	pointsin<-c(pointsin[1:(r-1)],remove,pointsin[r:N])
	}

#pointsin<-as.row(pointsin);
#lengths<-as.row(lengths);

lengthlist[[j]]<-lengths

		}#end for

return(list(X=X,coeff=coeff,lengths=lengths,lengthsremove=lengthsremove,pointsin=pointsin,removelist=removelist));
}

