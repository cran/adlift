"artlev1" <-
function(y,rem){

#returns vector of lengthsremove indices (p_i) denoting starting point (detail) of each level 

#y is lengthsremove 

p<-list()
l<-length(y)
q<-.5
p[[1]]<- which(y<=quantile(y,q))
if (length(p[[1]])<10){
p[[1]]<-1:l
	}
else{
all<-p[[1]]
n<-2
endpoint<-0

while(!endpoint){
q<-(1+q)/2
r<-setdiff(which(y<=quantile(y,q)),all)

if(length(r)<10){
r<-setdiff(1:length(y),all)
p[[n]]<-sort(r)
endpoint<-1
} #end if (r)

else{
p[[n]]<-r
all<-union(all,r)

#if ((l-length(all))<10){
#p[[n]]<-sort(union(p[[n]],setdiff(1:l,all)))  

#makes sure the last level includes all points not already in a level

#endpoint<-1}
#else{
n<-n+1
} #end else


}  #end while (endpoint)
} #end else

for (i in 1:length(p)){
p[[i]]<-rem[p[[i]]]
}

p
}
