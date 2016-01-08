library(sp)
library(RPostgreSQL)
library(Spatially)
library(rgeos)
library(raster)

# connect to the DB as pull the demographic information
# when connected from the office replace 'localhost' for 'db.urban4m.com'
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname = "dev", host = "localhost", user="geo", pass="geo") 
rs <- dbSendQuery(con, "select pop_lt_5, pop_5_14, pop_15_17, pop_18_34, pop_35_54, pop_55_64, pop_ge_65, _u ,st_astext(shape) from pavan.acs_for_accumulo")
s <- as.data.frame(fetch(rs, n=-1))
dbDisconnect(con)

#This step converts the demographic data from a traditional dta frame to a spatial data frame
demo <- SptlyWKTToSpatialDataFrame(s, uid='_u', shape='st_astext')

#This is the input (learning) point from which the demographic characeristics will be learned
pt0 <- SpatialPointsDataFrame(cbind(-71.06087923049925,42.360288053001966), 
                              data=data.frame(id='a'), 
                              proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
#Now we transform the learning point to our coordinate system
pt0 <- spTransform(pt0,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))

#create a 500m buffer around the learning point
buf0 <- SpatialPolygonsDataFrame(gBuffer(pt0,width=500,byid=T), data=data.frame(id='1'))

#A list describing the variables of the desired demographic characteristic.
#Must match the names in the demography spatial data frame
groups <- list(c('pop_lt_5', 'pop_5_14', 'pop_15_17', 'pop_18_34', 'pop_35_54', 'pop_55_64', 'pop_ge_65'))

#Adds columns with the demographic variables presented as percentage of the total
#(percentage of the sum of values in the demographic variables defined in "groups")
for(i in groups){
  nn <- paste0('perc_',unlist(i))
  j <- unlist(lapply(i,function(x){which(x==colnames(demo@data))}))
  x <- as.data.frame(demo@data[,j]/rowSums(as.data.frame(demo@data[,j])))
  colnames(x) <- nn
  demo@data <- cbind(demo@data,x)
}

#Empty entity to dump the results 
XX <- NULL
#Intersection of the demographic polygons with the buffer of the learning point
o1 <- gIntersection(demo, buf0, byid=TRUE, drop_lower_td=TRUE)
#Area of each intersected element
a <- gArea(o1, byid=TRUE)
#Id of the demographic polygons that intersect the buffer
n <- unlist(lapply(names(a),function(x){strsplit(x,' ')[[1]][1]}))
#Some more empty entities do dump results in
xx <- NULL
#Now we do some operations on each demographic polygon that intersect the buffer
for(h in 1:length(n)){
  #What is the h-th demographic polygon that intersects the buffer?
  k <- which(demo@data[,'_u']==n[h])
  #What's the ratio of area intersected by the buffer in that specific demographic polygon?
  aa <- a[h]/gArea(demo[k,])
  #Areal interpolation of the percentages in each demographic group
  x <- demo@data[k,which(unlist(gregexpr('perc',colnames(demo@data)))==1)]*aa
  nn <- paste0('perc_norm_',colnames(x))
  colnames(x) <- nn
  rbind(xx,x) -> xx
}
#summarize the normalized percentages (colmeans)
#This is THE spatial signature used as learning
x0 <- colMeans(xx,na.rm=T)

### NOW THAT WE HAVE LEARNED SOMETHING ABOUT ONE LOCATION, WE WILL SCOUT A GRID OF POINTS
##SEARCHING FOR SIMILARITIES TO WHAT WE'VE LEARNED

#read a file that contains the coordinates of the search grid
r <- read.table('/Users/sebastian/Documents/pointsMatrix',sep='\t')
r$id <- seq(1:nrow(r))
colnames(r)[1] <- 'shape'
grid <- SptlyWKTToSpatialDataFrame(r, uid='id', shape='shape',
                                crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
grid <- spTransform(grid,CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
grid <- gBuffer(grid,width=500,byid=T)

XX <- NULL
for(i in 1:length(grid)){
  o1 <- gIntersection(demo, grid[i,], byid=TRUE, drop_lower_td=TRUE)
  a <- gArea(o1, byid=TRUE)
  n <- unlist(lapply(names(a),function(x){strsplit(x,' ')[[1]][1]}))
  xx <- NULL
  #Now we do some operations on each demographic polygon that intersect the buffer
  for(h in 1:length(n)){
    #What is the h-th demographic polygon that intersects the buffer?
    k <- which(demo@data[,'_u']==n[h])
    #What's the ratio of area intersected by the buffer in that specific demographic polygon?
    aa <- a[h]/gArea(demo[k,])
    #Areal interpolation of the percentages in each demographic group
    x <- demo@data[k,which(unlist(gregexpr('perc',colnames(demo@data)))==1)]*aa
    nn <- paste0('perc_norm_',colnames(x))
    colnames(x) <- nn
    rbind(xx,x) -> xx
  }
  #This is the estimated "signature" at each point in the grid
  xi <- colMeans(xx,na.rm=T)
  #Now we stack them all together
  rbind(XX,xi) -> XX
}
XX <- as.data.frame(XX)

#x0 and xi (difference between observed and estimated signature)
dX <- t(apply(XX,1,function(x){x-x0}))

dX <- rowSums(abs(dX))
dX <- (dX-min(dX))/(max(dX)-min(dX))
dX <- max(dX)-dX

grid@data$similarityScore <- dX
grid@data$coordsWKT <- r$shape
