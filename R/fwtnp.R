"fwtnp" <-
function(input,f,inputtype="points",nkeep=2,intercept=TRUE,initboundhandl="reflect",updateboundhandl="add",neighbours=1,closest=FALSE,LocalPred=LinearPred){

#does the 1D single coefficient lifting transform
#inputtype determines which procedure to follow in order to apply the transform
#input is either a vector of points, or a vector of interval lengths, together
#with a startpoint in its first entry. If startpoint is NA, default is used
#(input matches inputtype).
#nkeep is no. points to keep after transform is performed (at least 2).
#intercept is boolean, denoting whether regression with  intercept is required
#boundaryhandling determines which boundary policy to use.
#closest indicates the way to choose neighbours (boolean)
#neighbours is the number of neighbours on each side, or total number if 
#closest is chosen.



if (inputtype=="lengths"){
	output<-pts(input);
	X<-output$X;
	lengths<-output$lengths;
	origlengths<-lengths
			}
if (inputtype=="points"){
	X<-input;
	I<-intervals(X,initboundhandl)$i;  #creates interval endpoints with
				 	   #edge correction

	lengths<-lengthintervals(X,I,type="midpoints",neighbours,closest)$lengths;
				    #^^^^^^^^^^^^^^^^^ 
	origlengths<-lengths
			}

X<-as.row(X);
f<-as.row(f);
nkeep<-max(nkeep,2);  #ensures not too many steps are performed

n<-length(X);

removelist<-NULL;	#contains removed points
lengthsremove<-NULL;    #contains interval lengths of removed points
neighbrs<-list();       #list containing the neighbours of the removed 
			#point at each step

gamlist<-list()
predlist<-NULL
alphalist<-list()
lengthlist<-list()
schemehist<-NULL       #records the history of prediction method (linear,..)
interhist<-NULL    #records the history of intercept (T,F)
detailslist<-list()
clolist<-NULL
minindexlist<-NULL
mindetailslist<-list()
indiceslist<-list()
history<-list()

pointsin<-matrix(1:n,1,n);
pointsin<-pointsin[order(X)];

coeff<-f;

#if (nkeep==length(x)){
#print("this is a null decomposition...")
#return(list(X=X,coeff=coeff,origlengths=origlengths,pointsin=pointsin))
#}

#else{
for (j in 1:(n-nkeep)){
#print(j)
remove<-order(lengths)[1];   #finds (sorted)index of point to remove from list
remove<-pointsin[remove];    #finds initial index of point to remove

removelist[j]<-remove;       #updates list of removed points

#if (closest & (neighbours>=length(pointsin))){
#neighbours<-neighbours-1;
#			}; #decreases the number of closest neighbours in case
			   # pointsin contains too few points

out<-getnbrs(X,remove,pointsin,neighbours,closest);

nbrs<-out$n;
index<-out$index;


res<-LocalPred(pointsin,X,coeff,nbrs,remove,intercept,neighbours)

if(length(res)==2){    #(for AN) l is the corresponding AP info, 
l<-res[[1]]		#i.e. of length 10
newinfo<-res[[2]]
nbrs<-newinfo[[3]]
index<-newinfo[[4]]
clolist[j]<-newinfo[[1]]
minindexlist[j]<-newinfo[[2]]
mindetailslist[[j]]<-newinfo[[5]]
indiceslist[[j]]<-newinfo[[6]]
history[[j]]<-c(indiceslist[[j]][[minindexlist[j]]],minindexlist[j])
		}
else{l<-res}  #l has length either 6 (LP,QP,CP) or 10 (AP)


neighbrs[[j]]<-nbrs  #puts appropriate nbrs into neighbrs according to 
		     #LocalPred choice
#coeff<-l[[6]];
weights<-l[[4]];

pred<-l[[5]];
bhat<-l[[3]];
if(length(l)==6){scheme<-NULL;int<-NULL;details<-NULL}else{

scheme<-l[[8]];
int<-l[[7]];
details<-l[[9]]
}
coeff[remove]<-coeff[remove]-pred


if (inputtype=="points"){
	l1<-PointsUpdate(X,coeff,nbrs,index,remove,pointsin,weights,lengths,updateboundhandl);
	coeff<-l1$coeff;
	lengths<-l1$lengths;
	r<-l1$r;
	weights<-l1$weights;
	N<-l1$N;
	alpha<-l1$alpha;
		}


#if (inputtype=="lengths"){
#	l2<-LengthsUpdate(X,coeff,nbrs,index,remove,pointsin,weights,lengths);
#	coeff<-l2$coeff;
#	lengths<-l2$lengths;
#	r<-l2$r;
#	weights<-l2$weights;
#	N<-l2$N;
#	alpha<-l2$alpha;
#		}		

lengthsremove[j]<-lengths[r];
gamlist[[j]]<-weights
predlist[j]<-pred
alphalist[[j]]<-alpha
schemehist[j]<-scheme
interhist[j]<-int
detailslist[[j]]<-details

lengths<-lengths[setdiff(1:length(pointsin),r)];
pointsin<-setdiff(pointsin,remove);

lengthlist[[j]]<-lengths
}  #end j for loop

#} nkeep loop

N<-length(pointsin);

return(list(X=X,coeff=coeff,origlengths=origlengths,lengths=lengths,lengthsremove=lengthsremove,pointsin=pointsin,removelist=removelist,neighbrs=neighbrs,neighbours=neighbours,schemehist=schemehist,interhist=interhist,clolist=clolist,gamlist=gamlist,alphalist=alphalist));
}







