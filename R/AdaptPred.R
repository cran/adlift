"AdaptPred" <-
function(pointsin,X,coeff,nbrs,remove,intercept,neighbours){

#does local adaptive prediction for the point remove based on N 
#points (chooses method of prediction and intercept);

details<-NULL
results<-list()


intercept<-FALSE  #does prediction schemes with no intercept

out1<-LinearPred(pointsin,X,coeff,nbrs,remove,intercept)
pred1<-out1$pred
out2<-QuadPred(pointsin,X,coeff,nbrs,remove,intercept)
pred2<-out2$pred
out3<-CubicPred(pointsin,X,coeff,nbrs,remove,intercept)
pred3<-out3$pred

intercept<-TRUE
out4<-LinearPred(pointsin,X,coeff,nbrs,remove,intercept)
pred4<-out4$pred
out5<-QuadPred(pointsin,X,coeff,nbrs,remove,intercept)
pred5<-out5$pred
out6<-CubicPred(pointsin,X,coeff,nbrs,remove,intercept)
pred6<-out6$pred

details[1]<-coeff[remove]-pred1
details[2]<-coeff[remove]-pred2
details[3]<-coeff[remove]-pred3
details[4]<-coeff[remove]-pred4
details[5]<-coeff[remove]-pred5
details[6]<-coeff[remove]-pred6

minindex<-order(abs(details))[1]
pred<-coeff[remove]-details[minindex]
coeff[remove]<-details[minindex]
int<-TRUE
scheme<-NULL
if(minindex<=3){int<-FALSE}

if((minindex==1)|(minindex==4)){scheme<-"Linear"}
if((minindex==2)|(minindex==5)){scheme<-"Quad"}
if((minindex==3)|(minindex==6)){scheme<-"Cubic"}

wt<-paste("out",minindex,"$weights",sep="")


return(list(weights=eval(parse(text=wt)),pred=pred,coeff=coeff,int=int,scheme=scheme,details=details,minindex=minindex))

}
