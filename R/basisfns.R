"basisfns" <-
function(x,f,pred,neigh,int,clo,keep,plot.f=FALSE,plot.bas=FALSE,separate=FALSE){

#produces plots of PRIMAL wavelet (psi) basis functions from a specified transform

#since multiplying Winv with a delta vector gives the columns of Winv, 
#the basis func. vectors are just the columns of winv.


tmd<-transmatdual(x,f,pred,neigh,int,clo,keep)

basmat<-matrix(0,length(f),length(f))
maxv<-NULL
out<-tmd$out
schemehist<-out$schemehist
interhist<-out$interhist
pointsin<-out$pointsin   #so that we know which basis functions are the scaling functions
removelist<-out$removelist
w<-tmd$Wnew

#winv<-Rmatsolve(w)
winv<-solve(w)

schhist<-matrix(0,1,length(f))
inthist<-matrix(0,1,length(f))
#reorders schemehist, interhist to be in order of x

if (plot.f==TRUE){
plot(sort(x),f[order(x)],type="l")
orig<-dev.cur()

}

for (k in 1:length(f)){
colour<-1

basmat[k,]<-winv[,k]   #each ROW of basmat is basis fun  

if(is.null(schemehist)){schhist<-NULL}
else{
q<-which(removelist==k)  #finds out the scheme of removelist in order of x
if (length(q)!=0){
schhist[k]<-schemehist[q]
inthist[k]<-interhist[q]
}


if(schhist[k]=="Linear"){colour<-5}
if(schhist[k]=="Quad"){colour<-12}
if(schhist[k]=="Cubic"){colour<-6}
if(schhist[k]=="0"){colour<-3}
}

maxv[k]<-order(abs(basmat[k,]))[length(f)]

if (plot.f==TRUE){
dev.set(orig)
#lines(rep(x[maxv[k]],times=11),-5:5,type="l",col=colour)  
					#plots line at peak of basis fn
lines(rep(x[k],times=11),-5:5,type="l",col=colour)
					#plots line at datapoint (x axis)
}

#basmat values will be in the "order" of x[1:length(f)] , so still need to sort


} #end for 

c<-setdiff(plot.bas,0)

if (length(c)!=0){

if (any(is.na(match(plot.bas,1:length(f))))){
cat("can't plot requested basis functions")
}
else{
if (separate==FALSE){

#graphsheet()		#for R, splus
X11()

plot(x,f,type="n")
newdev<-dev.cur()
}

for (i in plot.bas){

if(is.null(schhist[i])){colour<-i}
else{
if(schhist[i]=="Linear"){colour<-5}
if(schhist[i]=="Quad"){colour<-12}
if(schhist[i]=="Cubic"){colour<-6}
if(schhist[i]=="0"){colour<-3}
}

if (separate==TRUE){
X11()
plot(sort(x),basmat[i,][order(x)],type="l",col=colour)
}
else{
dev.set(newdev)
lines(sort(x),basmat[i,][order(x)],type="l",col=colour)
}
} #end for
} #end else
} #end if

return(list(out=out,maxv=maxv,schhist=schhist,inthist=inthist,basmat=basmat))

}

