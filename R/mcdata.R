"mcdata" <-
function(){

#reads in motorcycle data from text file

data(motorcycledata)

times<<-as.column(motorcycledata[,1])
accel<<-as.column(motorcycledata[,2])

}

