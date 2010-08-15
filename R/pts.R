"pts" <-
function(input){

#produces X vector according to lengths, assigning values to left endpoints of 
#interval lengths, according to a given startpoint.

X<-NULL

startpoint<-input[1];
lengths<-input[2:length(input)];

n<-length(lengths);

if (is.na(startpoint)){
startpoint<-sum(lengths)/(n+1);
		}

X[1]<-startpoint;

X<-matrix(c(X[1],X[1]+cumsum(lengths)),nrow=1)

lengths[n+1]<-sum(lengths)/(n+1);

return(list(lengths=lengths,X=X));

}
