"denoisehetero2"<-function(x,f,pred=1,neigh=1,int=TRUE,clo=FALSE,keep=2,rule="median",verbose=TRUE,index=1:length(x)){

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

n<-length(x)
temp<-fwtnp2(x,f,pred,neigh,int,clo,keep)
#w<-temp$Wnew
out<-temp$out
x<-temp$x
lca<-out$lca

po<-out$pointsin
nonorcoeff<-out$coeff
lr<-out$lengthsremove    #vector deciding how to divide up coefficients into artificial levels
rem<-lca[,1] 		     #used to convert output to original lr,rem)

al<-artlev(lr,rem)      #the list of indices of removelist separated into levels
levno<-length(al)


detail<-y<-matrix(0,1,n-keep)

y<-x[setdiff(1:n,po)]
detail<-nonorcoeff[setdiff(1:n,po)]
h<-heterovar(y,detail,al)

sdvec<-h$varvec
sdvec1<-h$varvec1
sdvec2<-h$varvec2

sd<-sd1<-sd2<-NULL
m<-min(po)
M<-max(po)

#i<m:
sd[1:(m-1)]<-sdvec[1:(m-1)]
sd1[1:(m-1)]<-sdvec1[1:(m-1)]
sd2[1:(m-1)]<-sdvec2[1:(m-1)]

#i=m, i=M:
sd[c(m,M)]<-sd1[c(m,M)]<-sd2[c(m,M)]<-NA

#i>m & i<M:
sd[(m+1):(M-1)]<-sdvec[m:(M-2)]
sd1[(m+1):(M-1)]<-sdvec1[m:(M-2)]
sd2[(m+1):(M-1)]<-sdvec2[m:(M-2)]

#i>M:
sd[(M+1):n]<-sdvec[(M-1):(n-2)]
sd1[(M+1):n]<-sdvec1[(M-1):(n-2)]
sd2[(M+1):n]<-sdvec2[(M-1):(n-2)]

for (i in 1:levno){
	nordetlist[[i]]<-nonorcoeff[al[[i]]]/(sd[al[[i]]])
	nordetlist1[[i]]<-nonorcoeff[al[[i]]]/(sd1[al[[i]]])
	nordetlist2[[i]]<-nonorcoeff[al[[i]]]/(sd2[al[[i]]])

	nortclist[[i]]<-ebayesthresh(nordetlist[[i]],prior="cauchy",a=NA,sdev=1,threshrule=rule)
	nortclist1[[i]]<-ebayesthresh(nordetlist1[[i]],prior="cauchy",a=NA,sdev=1,threshrule=rule)
	nortclist2[[i]]<-ebayesthresh(nordetlist2[[i]],prior="cauchy",a=NA,sdev=1,threshrule=rule)

	newcoeff[al[[i]]]<-nortclist[[i]]*(sd[al[[i]]])
	newcoeff1[al[[i]]]<-nortclist1[[i]]*(sd1[al[[i]]])
	newcoeff2[al[[i]]]<-nortclist2[[i]]*(sd2[al[[i]]])
}

newcoeff[po]<-out$coeff[po]
newcoeff1[po]<-out$coeff[po]
newcoeff2[po]<-out$coeff[po]

neighbrs<-list()
for(kk in 1:length(rem)){
        neighbrs[[kk]]<-.C("nbrsfromlca",as.double(t(lca)),as.integer(ncol(lca)),as.integer(kk),
        n=as.integer(rep(0,times=lca[kk,2])),PACKAGE="adlift")$n
}
                                                                                
adds <- findadds(rem, neighbrs, po, index)
add <- max(adds)
if (verbose) {
	cat("doing ", add, " steps...\n", sep = "")
}

fhat <- invtnp2(x, newcoeff, out$lengths, lr, po,lca, add)
fhat1 <- invtnp2(x, newcoeff1, out$lengths, lr, po,lca, add)
fhat2 <- invtnp2(x, newcoeff2, out$lengths, lr, po,lca, add)
                                                                                                                    
return(list(out=out,al=al,sd=sd,fhat=fhat,fhat1=fhat1,fhat2=fhat2))
}
