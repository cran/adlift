"PointsUpdate" <-
function(X,coeff,nbrs,index,remove,pointsin,weights,lengths,updateboundhandl){

#does the update lifting step based on nbrs of remove

r<-which(pointsin==remove);

N<-length(pointsin);

###update the interval lengths (>=2 nbrs)###

if ((r>=2)&(r<=(N-1))){
	lengths[index]<-as.row(lengths[index]);
	weights<-as.row(weights);
	lengths[index]<-lengths[index]+lengths[r]*weights;

	}
else{
#	if(updateboundhandl=="reflect"){
#
#		if(r==1){
#			lengths[2]<-X[pointsin[3]]-X[pointsin[2]];
#			}
#		if(r==N){
#			lengths[N-1]<-X[pointsin[N-1]]-X[pointsin[N-2]];
#			}
#				};  #end reflect
#
#	if(updateboundhandl=="stop"){
#
#	  if(r==1){
#	      lengths[2]<-length[2]+length[1]-(X[pointsin[2]]-X[pointsin[1]]);
#		    }
#	  if(r==N){
#	    lengths[N-1]<-length[N]+length[N-1]-(X[pointsin[N]]-X[pointsin[N-1]])/2;
#			}
#				};  #end stop
#
	if(updateboundhandl=="add"){

		if(r==1){
			lengths[2]<-lengths[2]+lengths[1];
			}
		if(r==N){
			lengths[N-1]<-lengths[N-1]+lengths[N];
			}
				};  #end add

    }  #end else

###update the scaling function coefficients###

alpha<-matrix(0,1,length(nbrs));

if (length(nbrs)>=2){
alpha<-lengths[r]*lengths[index]/(sum(lengths[index]^2));

			#these are the update weights for scaling coefficients 



coeff[pointsin[index]]<-coeff[pointsin[index]]+alpha*coeff[remove];
		}
else{
	q<-which(pointsin==nbrs);
	alpha<-lengths[r]/lengths[q];
	coeff[pointsin[q]]<-coeff[pointsin[q]]+alpha*coeff[remove];
}

return(list(coeff=coeff,lengths=lengths,r=r,N=N,weights=weights,alpha=alpha));
}
