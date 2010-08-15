"fwtnpmp" <-
function(input,f,inputtype="points",nkeep=2,intercept=TRUE,initboundhandl="reflect",neighbours=1,closest=FALSE,LocalPredmp=LinearPredmp,mpdet="ave"){

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
#mpdet is either "ave" (averaged multiple details) or "min" for the minimum one


xold<-NULL
fold<-NULL
g<-list()
X<-NULL


xold<-input
fold<-f

isu<-adjustx(xold,fold,"mean")
X<-isu$sepx
f<-isu$sepf
g<-isu$groups


coefflist<-list()
mp<-matrix(0,1,length(g))

for (i in 1:length(g)){
coefflist[[i]]<-fold[g[[i]]]
mp[i]<-(length(g[[i]])>1)
}


#if (inputtype=="lengths"){
#	output<-pts(X);
#	X<-output$X;
#	lengths<-output$lengths;
#	origlengths<-lengths
#			}
if (inputtype=="points"){
	X<-X;
	I<-intervals(X,initboundhandl);  #creates interval endpoints with
				 	   #edge correction

	lengths<-lengthintervals(X,I,type="midpoints",neighbours,closest);
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
newneighbrs<-list()

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

for (j in 1:(n-nkeep)){
#print(j)
remove<-order(lengths)[1];   #finds (sorted)index of point to remove from list
remove<-pointsin[remove];    #finds initial index of point to remove

removelist[j]<-remove;       #updates list of removed points

#if (closest & (neighbours>=length(pointsin))){
#neighbours<-neighbours-1;
#			}; #decreases the number of closest neighbours in case
			   # pointsin contains too few points


d<-matrix(0,1,length(coefflist[[remove]]))

out<-getnbrs(X,remove,pointsin,neighbours,closest);

nbrs<-out$n;
index<-out$index;

newnbrs<-NULL
for (i in 1:length(nbrs)){
newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
}


res<-LocalPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)

if(length(res)==2){
l<-res[[1]]
newinfo<-res[[2]]
nbrs<-newinfo[[3]]
index<-newinfo[[4]]
clolist[j]<-newinfo[[1]]
minindexlist[j]<-newinfo[[2]]
mindetailslist[[j]]<-newinfo[[5]]
indiceslist[[j]]<-newinfo[[6]]
history[[j]]<-c(indiceslist[[j]][[minindexlist[j]]],minindexlist[j])
		}
else{l<-res}

if(length(res)==2){
newnbrs<-NULL
for (i in 1:length(nbrs)){
newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
}
}
newneighbrs[[j]]<-newnbrs
neighbrs[[j]]<-nbrs  #puts appropriate nbrs into neighbrs according to 
		     #LocalPred choice
coeff<-l[[6]];
weights<-l[[4]];

pred<-l[[5]];
bhat<-l[[3]];
if(length(l)==6){scheme<-NULL;int<-NULL;details<-NULL}else{

scheme<-l[[8]];
int<-l[[7]];
details<-l[[9]]
}

if (length(res)==6){# i.e. is lp,qp,cp
	if (length(pred)>1){
		md<-matrix(0,1,length(coefflist[[remove]]))
		pr<-matrix(0,1,length(coefflist[[remove]]))
		if (mpdet=="min"){
			for (i in 1:length(coefflist[[remove]])){
				pr[i]<-order(abs(coefflist[[remove]][i]-pred))[1]
				md[i]<-(coefflist[[remove]][i]-pred)[pr[i]]
			}#end for 
		}
		else{
			for (i in 1:length(coefflist[[remove]])){
				md[i]<-mean(coefflist[[remove]][i]-pred)
			}
		}

		if (mpdet=="min"){
			sel<-order(abs(md))[1]
			coefflist[[remove]]<-md[sel]
			pred<-pred[pr[sel]]
			coeff[remove]<-coefflist[[remove]]
		}
		else{
			coefflist[[remove]]<-mean(md)
			coeff[remove]<-coefflist[[remove]]
			pred<-mean(pred)
		}


	}#end if

	else{

		for (i in 1:length(coefflist[[remove]])){
			d[i]<-coefflist[[remove]][i]-pred
		}
		aved<-mean(d)
		mind<-min(d)

		if (mpdet=="min"){
			coefflist[[remove]]<-mind
		}
		else{
			coefflist[[remove]]<-aved
		}
		coeff[remove]<-coefflist[[remove]]

	} #end else (pred)
}
else{
	coefflist[[remove]]<-coeff[remove]
}


if (inputtype=="points"){
	l1<-PointsUpdatemp(X,coefflist,nbrs,newnbrs,index,remove,pointsin,weights,lengths);
	coefflist<-l1$coeff;
	lengths<-l1$lengths;
	r<-l1$r;
	weights<-l1$weights;
	N<-l1$N;
	alpha<-l1$alpha;
		}

coeff[nbrs]<-coeff[nbrs]+alpha*coeff[remove];
	

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

N<-length(pointsin);

for (i in 1:length(pointsin)){
coefflist[[pointsin[i]]]<-mean(coefflist[[pointsin[i]]])
coeff[pointsin[i]]<-coefflist[[pointsin[i]]]
}

return(list(X=X,coeff=coeff,coefflist=coefflist,origlengths=origlengths,lengths=lengths,lengthsremove=lengthsremove,pointsin=pointsin,removelist=removelist,neighbrs=neighbrs,newneighbrs=newneighbrs,neighbours=neighbours,schemehist=schemehist,interhist=interhist,clolist=clolist,gamlist=gamlist,alphalist=alphalist,g=g,mp=mp));
}
