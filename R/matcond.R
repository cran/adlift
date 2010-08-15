"matcond" <-
function(x,f,Pred,neigh,int,clo,keep){

 a<-transmatdual(x,f,Pred,neigh,int,clo,keep)
 W<-a$Wnew

 cno<-condno(W,type="F")
 
 v<-svd(W)$d[1]/(svd(W)$d[nrow(W)])


return(list(cno=cno,v=v,a=a))

}
