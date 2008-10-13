"Amatdual" <-
function(steps,pointsin,removelist,nbrs,weights,alpha){

#pointsin is pointsin after chosen decomposition
#steps=number of removed points from decomposition

n<-length(pointsin)+length(removelist)

Adual<-matrix(0,n-steps+1,n-steps+1)
Hdual<-matrix(0,n-steps,n-steps+1)
Gdual<-matrix(0,1,n-steps+1)

newpoints<-(c(pointsin,rev(removelist)))[1:(n-steps+1)]

o<-NULL
for (i in 1:length(nbrs)){
o[i]<-which(newpoints[1:(length(newpoints)-1)]==nbrs[i])
	}

lastcol<-matrix(0,length(newpoints)-1,1)
lastcol[o]<-alpha
Hdual<-cbind(diag(length(newpoints)-1),lastcol)
for (i in 1:length(o)){
for (j in 1:length(o)){
Hdual[o[i],o[j]]<- -alpha[i]*weights[j]
		}
		}
for (i in 1:length(o)){
Hdual[o[i],o[i]]<-Hdual[o[i],o[i]]+1
			}

Gdual[o]<- -weights
Gdual[length(newpoints)]<-1

Adual<-rbind(Hdual,as.row(Gdual))

return(list(Adual=Adual,Hdual=Hdual,Gdual=Gdual,o=o,alpha=alpha,weights=weights))

}

