"pts" <-
function(input){

#produces X vector according to lengths, assigning values to left endpoints of 
#interval lengths, according to a given startpoint.

startpoint<-input[1];
lengths<-input[2:length(input)];

n<-length(lengths);

if (is.na(startpoint)){
startpoint<-sum(lengths)/(n+1);
		}


X<-matrix(0,1,n+1);

X[1]<-startpoint;

for (i in 1:n){
	X[i+1]<-X[i]+lengths[i];
	}

lengths[n+1]<-sum(lengths)/(n+1);

return(list(lengths=lengths,X=X));

}
