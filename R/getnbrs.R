"getnbrs" <-
function(X,remove,pointsin,neighbours,closest){

#finds the indices of the neighbours of X[remove]  

#I is (ordered) list of interval endpoints
#neighbours is no. neighbours on each side...
#the procedure reduces the number of neighbours on a respective side if they
#don't exist.

N<-length(pointsin);
d<-neighbours;
nbrs<-NULL;
r<-which(remove==pointsin);
range<-2:(N-1);

if (is.element(r,range)){     	#this loop checks whether d neighbours exist
	checkl<-matrix(0,1,d);  #on each side
	for (i in 1:d){
		checkl[i]<-(r-i<1); #checks the left...
		      };	              
				#... and counts how many places fail
	leftneigh<-d-sum(checkl);  
		               
	checkr<-matrix(0,1,d);
	for (i in 1:d){
		checkr[i]<-(r+i>N); #checks the right...
		      };
	rightneigh<-d-sum(checkr);     #counts no. places which fail...
	
	for (i in 1:leftneigh){nbrs[i]<-pointsin[r-leftneigh+i-1]};
	for (j in 1:rightneigh){nbrs[j+leftneigh]<-pointsin[r+j]};	

			};  #end if 


if (!closest){          #same as not in range?

	if (r==1){
		nbrs<-pointsin[2];
		leftneigh<-0;
		rightneigh<-1;
		};
	if (r==N){
		nbrs<-pointsin[N-1];
		leftneigh<-1;
		rightneigh<-0;
		};
	index<-setdiff((r-leftneigh):(r+rightneigh),r);
	}  #end !closest
else{
	if (r==1){
		index<-r+1
#		index<-(r+1):(r+d)
		nbrs<-pointsin[index];
		leftneigh<-0;
#		rightneigh<-d;
		rightneigh<-1
		}
	if (r==N){
		index<-(r-1)
#		index<-(r-d):(r-1)
		nbrs<-pointsin[index];
#		leftneigh<-d;
		leftneigh<-1
		rightneigh<-0;
		}
	if (is.element(r,range)){ 
		distances<-matrix(0,1,leftneigh+rightneigh);	
		for (i in 1: leftneigh){
		    distances[i]<-abs(X[remove]-X[pointsin[r-leftneigh+i-1]]);
					}
		for (j in 1:rightneigh){
		    distances[j+leftneigh]<-abs(X[remove]-X[pointsin[r+j]]);
					}
		d1<-min(d,N-1)		
		q<-order(distances)[1:d1];	
		nbrs<-nbrs[q];
		index<-setdiff((r-leftneigh):(r+rightneigh),r);
		index<-index[q];
		B<-matrix(0,1,d1);
		for (i in 1:d1){  
			B[i]<-(index[q[i]]<r);    
				};
		numleft<-sum(B);    #counts number of nbrs on left of remove
		leftneigh<-numleft;
		rightneigh<-d1-leftneigh;
		#index<-index[q];
				}; #end if is.element

}  #end else

index<-sort(index)
nbrs<-pointsin[index]
nbrs<-as.vector(nbrs);

return(list(nbrs=nbrs,index=index));



}



