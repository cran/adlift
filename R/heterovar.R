"heterovar" <-
function(y,detail,al){

#works out the individual noise variances

idlist<-list()
newidlist<-list()
dlist<-list()
varvec<-matrix(0,1,length(y))
varvec1<-matrix(0,1,length(y))
varvec2<-matrix(0,1,length(y))

winlength<-(max(y)-min(y))/5

#newdetail<-detail[order(y)]

ch<-matrix(1,1,length(detail))

while(sum(ch)!=0){

ep<-matrix(0,length(y),2)

ep[,1]<-y-winlength/2
ep[,2]<-y+winlength/2

ep1<-ep

q<-which(ep[,1]<min(y))
r<-which(ep[,2]>max(y))

for (i in q){
ep[q,1]<-min(y)
ep[q,2]<-min(y)+winlength
}
for (i in r){
ep[r,2]<-max(y)
ep[r,1]<-max(y)-winlength
}

ep2<-ep

for (i in 1:length(detail)){

idlist[[i]]<-which((y>=ep[i,1])&(y<=ep[i,2]))
newidlist[[i]]<-idlist[[i]][is.element(idlist[[i]],al[[1]])]
dlist[[i]]<-detail[newidlist[[i]]]
varvec[i]<-mad(dlist[[i]])
varvec1[i]<-mad(dlist[[i]],center=0)
varvec2[i]<-median(abs(dlist[[i]]))

}
for (k in 1:length(newidlist)){
ch[k]<-(length(newidlist[[k]])<4)
}

#ch[2]<-(is.element(0,varvec1))
#ch[3]<-(is.element(0,varvec2))

if (sum(ch)!=0){

winlength<-winlength+(max(y)-min(y))/20

}

}

return(list(ep1=ep1,ep2=ep2,idlist=idlist,newidlist=newidlist,dlist=dlist,varvec=varvec,varvec1=varvec1,varvec2=varvec2))

}
