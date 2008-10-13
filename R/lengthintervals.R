"lengthintervals" <-
function(X,I,type="midpoints",neighbours,closest){

#computes the lengths of the intervals, given an intervals vector I
#if using closest neighbour method, option: lengths as local neighbour distance
#average

o<-as.column(order(X));

lengths<-matrix(0,1,length(X));   #empty lengths vector

d<-neighbours;

if (closest){
initialnbrs<-matrix(0,length(X),neighbours);  #sets up empty matrices to hold
						# nbrs,
initialindex<-matrix(0,length(X),d);		# indices

if (type=="average"){
	for (i in 1:length(X)){
		out<-getnbrs(X,o[i],order(X),neighbours,closest);
							 #find nbrs of point i
		initialnbrs[i,]<-out$nbrs;
		initialindex[i,]<-out$index;
			}

    for (i in 1:length(X)){	
	 lengths[i]<-sum(abs(rep(X[o[i]],times=d)-X[initialnbrs[i,]]))/d;
			}
		};

if (type=="midpoints"){
	for (i in 1:(length(X))){
		lengths[i]<-I[i+1]-I[i];
				}
			    
	initialnbrs<-NULL;
	initialindex<-NULL;
		}

}
else{

	for (i in 1:(length(X))){
		lengths[i]<-I[i+1]-I[i];
				}
			    
	initialnbrs<-NULL;
	initialindex<-NULL;

}

initialnbrs<-cbind(o,initialnbrs);
return(list(lengths=lengths,initialnbrs=initialnbrs,initialindex=initialindex));

}
