"AdaptPredmp" <-
function(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g){

#does local adaptive prediction for the point remove based on N 
#points (chooses method of prediction and intercept);

#library(Matrix);

details<-NULL
results<-list()
d<-matrix(0,1,length(coefflist[[remove]]))

intercept<-FALSE  #does prediction schemes with no intercept

out1<<-LinearPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
pred1<-out1$pred
out2<-QuadPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
pred2<-out2$pred
out3<-CubicPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
pred3<-out3$pred

intercept<-TRUE
out4<-LinearPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
pred4<-out4$pred
out5<-QuadPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
pred5<-out5$pred
out6<-CubicPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
pred6<-out6$pred

p<-list()
p[[1]]<-pred1
p[[2]]<-pred2
p[[3]]<-pred3
p[[4]]<-pred4
p[[5]]<-pred5
p[[6]]<-pred6
#print(p)
pre<-matrix(0,1,6)
for (k in 1:6){
pred<-p[[k]]
#print(k)
#print("set")
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
		details[k]<-md[sel]
		pre[k]<-pred[pr[sel]]
		#coeff[remove]<-coefflist[[remove]]
	}
	else{
		details[k]<-mean(md)
		#coeff[remove]<-coefflist[[remove]]
		pre[k]<-mean(pred)
	}


}#end if

else{
#print(coefflist[[remove]])
#print(pred)
	for (i in 1:length(coefflist[[remove]])){
		d[i]<-coefflist[[remove]][i]-pred
	}
#print(d)
	aved<-mean(d)
	mind<-min(d)
	pre[k]<-pred	
#print(pre)
	if (mpdet=="min"){
		details[k]<-mind
	}
	else{
		details[k]<-aved
	}
	#coeff[remove]<-coefflist[[remove]]

} #end else (pred)
#print(details)

#details[1]<-coeff[remove]-pred1
#details[2]<-coeff[remove]-pred2
#details[3]<-coeff[remove]-pred3
#details[4]<-coeff[remove]-pred4
#details[5]<-coeff[remove]-pred5
#details[6]<-coeff[remove]-pred6

} #for k

minindex<-order(abs(details))[1]
pred<-pre[minindex]
coefflist[[remove]]<-details[minindex]
coeff[remove]<-coefflist[[remove]]
int<-NULL
scheme<-NULL
if(minindex<=3){int<-FALSE}
else{int<-TRUE}


if((minindex==1)|(minindex==4)){scheme<-"Linear"}
if((minindex==2)|(minindex==5)){scheme<-"Quad"}
if((minindex==3)|(minindex==6)){scheme<-"Cubic"}



if(minindex==1){
	for (i in 1:4){
results[[i]]<-out1[[i]]
	}
}
if(minindex==2){
	for (i in 1:4){
results[[i]]<-out2[[i]]
	}
}
if(minindex==3){
	for (i in 1:4){
results[[i]]<-out3[[i]]
	}
}
if(minindex==4){
	for (i in 1:4){
results[[i]]<-out4[[i]]
	}
}
if(minindex==5){
	for (i in 1:4){
results[[i]]<-out5[[i]]
	}
}
if(minindex==6){
	for (i in 1:4){
results[[i]]<-out6[[i]]
	}
}


results[[5]]<-pred
results[[6]]<-coeff
results[[7]]<-int
results[[8]]<-scheme
results[[9]]<-details
results[[10]]<-minindex


results
}
