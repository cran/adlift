"denoisehetero" <-
function(x,f,pred,neigh,int,clo,keep,rule="median"){


nonorcoeff<-list()
nordetlist<-list()
nordetlist1<-list()
nordetlist2<-list()
nortclist<-list()
nortclist1<-list()
nortclist2<-list()
al<-list()
temp<-list()
sdvec<-NULL
sdvec1<-NULL
sdvec2<-NULL
newcoeff<-NULL
newcoeff1<-NULL
newcoeff2<-NULL
fhat<-list()
fhat1<-list()
fhat2<-list()

temp<-transmatdual(x,f,pred,neigh,int,clo,keep)
w<-temp$Wnew
out<-temp$out
x<-temp$x

print(length(x))
print(length(out$pointsin))
nonorcoeff<-out$coeff
lr<-out$lengthsremove    #vector deciding how to divide up coefficients into artificial levels
rem<-out$removelist      #used to convert output to original lr,rem)

al<-artlev1(lr,rem)      #the list of indices of removelist separated into levels
levno<-length(al)


y<-matrix(0,1,(length(x)-keep))
detail<-matrix(0,1,(length(x)-keep))

y<-x[setdiff((1:length(x)),out$pointsin)]
detail<-nonorcoeff[setdiff((1:length(x)),out$pointsin)]
h<-heterovar(y,detail,al)


sdvec<-h$varvec
sdvec1<-h$varvec1
sdvec2<-h$varvec2

sd<-NULL
for (i in 1:length(x)){
if (i<min(out$pointsin)){
	sd[i]<-sdvec[i]}
if (i==min(out$pointsin)){
	sd[i]<-NA}
if ((i>min(out$pointsin)) & (i<max(out$pointsin))){
	sd[i]<-sdvec[i-1]}
if(i==max(out$pointsin)){
	sd[i]<-NA}
if(i>max(out$pointsin)){
	sd[i]<-sdvec[i-2]}
}

sd1<-NULL
for (i in 1:length(x)){
if (i<min(out$pointsin)){
	sd1[i]<-sdvec1[i]}
if (i==min(out$pointsin)){
	sd1[i]<-NA}
if ((i>min(out$pointsin)) & (i<max(out$pointsin))){
	sd1[i]<-sdvec1[i-1]}
if(i==max(out$pointsin)){
	sd1[i]<-NA}
if(i>max(out$pointsin)){
	sd1[i]<-sdvec1[i-2]}
}

sd2<-NULL
for (i in 1:length(x)){
if (i<min(out$pointsin)){
	sd2[i]<-sdvec2[i]}
if (i==min(out$pointsin)){
	sd2[i]<-NA}
if ((i>min(out$pointsin)) & (i<max(out$pointsin))){
	sd2[i]<-sdvec2[i-1]}
if(i==max(out$pointsin)){
	sd2[i]<-NA}
if(i>max(out$pointsin)){
	sd2[i]<-sdvec2[i-2]}
}

for (i in 1:levno){
nordetlist[[i]]<-nonorcoeff[al[[i]]]/(sd[al[[i]]])
nordetlist1[[i]]<-nonorcoeff[al[[i]]]/(sd1[al[[i]]])
nordetlist2[[i]]<-nonorcoeff[al[[i]]]/(sd2[al[[i]]])

	}

for (i in 1:levno){
nortclist[[i]]<-ebayesthresh(nordetlist[[i]],prior="cauchy",a=NA,sdev=1,threshrule=rule)
nortclist1[[i]]<-ebayesthresh(nordetlist1[[i]],prior="cauchy",a=NA,sdev=1,threshrule=rule)
nortclist2[[i]]<-ebayesthresh(nordetlist2[[i]],prior="cauchy",a=NA,sdev=1,threshrule=rule)
	}


for (i in 1:levno){
newcoeff[al[[i]]]<-nortclist[[i]]*(sd[al[[i]]])
newcoeff1[al[[i]]]<-nortclist1[[i]]*(sd1[al[[i]]])
newcoeff2[al[[i]]]<-nortclist2[[i]]*(sd2[al[[i]]])
		}
newcoeff[out$pointsin]<-out$coeff[out$pointsin]
newcoeff1[out$pointsin]<-out$coeff[out$pointsin]
newcoeff2[out$pointsin]<-out$coeff[out$pointsin]

fhat<-invtnp(x,newcoeff,out$lengths,lr,out$pointsin,rem,out$neighbrs,out$schemehist,out$interhist,length(x)-keep,int,neigh,clo,pred)
fhat1<-invtnp(x,newcoeff1,out$lengths,lr,out$pointsin,rem,out$neighbrs,out$schemehist,out$interhist,length(x)-keep,int,neigh,clo,pred)
fhat2<-invtnp(x,newcoeff2,out$lengths,lr,out$pointsin,rem,out$neighbrs,out$schemehist,out$interhist,length(x)-keep,int,neigh,clo,pred)

return(list(out=out,w=w,al=al,sd=sd,fhat=fhat,fhat1=fhat1,fhat2=fhat2))
}
