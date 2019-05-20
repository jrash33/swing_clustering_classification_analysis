rm(list=ls())

#functions
LoadLibraries=function(){
  library(caret)
  library(AppliedPredictiveModeling)
  library(corrplot)
  library(rpart) 
  library(C50)
  library(sfsmisc)
  library(reshape2)
  library(stats)
  library(cluster)
  library(factoextra)
  library(timeDate)
  library(dendextend)
  print("The Libraries have been loaded!")
}

LoadLibraries()


#--------------------- Load Data and configure for modeling -----------------------------------
#open csv file
data=read.csv("data_FR.csv")

#some small transformations

#remove any missing data (NaN)
data = na.omit(data)

#delete players/club columns
data=data[,-1:-2]
col=ncol(data)


#------------------------- Pre-Processing: Normalize/Center --------------------------------
#skewness
skew=rep(NA,col)

for (i in 1:col){
  skew[i]=skewness(data[,i])
  i=i+1
}

#skewness isn't actually that bad
hist(skew)
plot(skew)

skewValues=apply(data,2,skewness)
head(skewValues)

#center and scale data
data_scaled = scale(data,center=TRUE,scale=TRUE)


#------------------------ Clustering Methods ------------------------------

#******HIERARCHIAL CLUSTERING******

#agglomerative hierarchical clustering (HC): computes all pairwise dissimilarities between the elements
# in cluster 1 and the elements in cluster 2, and considers the largest value of these dissimilarities
#as the distance between the 2 clusters


#dissimilarity matrix
d=dist(data_scaled,method="euclidean")

#hierarchical clustering using complete linkage 
HC4=hclust(d,method="ward.D")
plot(HC4)
HC5=hclust(d,method="ward.D2")
plot(HC5)


#compute with agnes method to get the agglomerative function 
#which measures the amount of clustering structure found (values close to 1 are best)
#best method based on agglomerative coefficient
HC14=agnes(data_scaled,method="ward")
HC14$ac
pltree(HC14,cex=.6,hang=-1,main="Dendrogram of agnes (ward)")


#identify clusters/sub-groups with cutree
#cut tree into 3 groups
sub_grp4=cutree(HC4,k=3)
sub_grp5=cutree(HC5,k=3)

#number of members in each cluster
table(sub_grp4)
table(sub_grp5)

#******Moving forward  we will analyze Sub Groups of the ward method******
#way to see all sub-groups!
counts=sapply(2:8,function(ncl)table(cutree(HC4,ncl)))
names(counts)=2:8
counts

counts=sapply(2:8,function(ncl)table(cutree(HC5,ncl)))
names(counts)=2:8
counts

#conclusion: there is no obvious number of clusters--> trial and error from here


#we can add a cluster column to our original data
#AVGdata %>% mutate(cluster=sub_grp1) %>% head

#observe clusters
plot(HC4)
rect.hclust(HC4,k=4,border=2:5)
fviz_cluster(list(data=data_scaled,cluster=sub_grp4))
clusplot(data_scaled,sub_grp4,color=TRUE,shade=TRUE)

plot(HC5)
rect.hclust(HC5,k=4,border=2:5)
fviz_cluster(list(data=data_scaled,cluster=sub_grp5))
clusplot(data_scaled,sub_grp5,color=TRUE,shade=TRUE)


#we can also compare 2 dendrograms: compare hierarchical clustering with complete linkage vs ward's method
#tanglegram plots 2 dendrograms side by side with their labels connected by lines

dend1=as.dendrogram(HC4)
dend2=as.dendrogram(HC5)

#tanglegram plot as visual
tanglegram(dend1,dend2)

#we can measure the quality of the alignment of these 2 trees (1 is no match/0 is perfect match)
dend_list=dendlist(dend1,dend2)
tanglegram(dend1,dend2,
           highlight_distinct_edges=FALSE, #turn off dashed lines
           common_subtrees_color_lines = FALSE, #turn off line colors
           common_subtrees_color_branches = TRUE, #color common branches
           main=paste("entanglement=",round(entanglement(dend_list),2)))


#elbow method
fviz_nbclust(data_scaled,FUN=hcut,method="wss")

#avg silhouette method
fviz_nbclust(data_scaled,FUN=hcut,method="silhouette")

#gap statistic method
gap_stat = clusGap(data_scaled, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)



#look at each cluster group
data_scaled[sub_grp4==1,]




#cluster characterization

#look at the medians of each parameter for each cluster
CI_upper = function(matrix){
  
  #get parameters
  avg = mean(matrix)
  st_dev = sd(matrix)
  n = length(matrix)
  
  #calculate a confidence interval from a t distribution
  interval = qt(.975, df=n-1)*st_dev/sqrt(n)
  
  return(avg+st_dev)
}

#look at the medians of each parameter for each cluster
CI_lower = function(matrix){
  
  #get parameters
  avg = mean(matrix)
  st_dev = sd(matrix)
  n = length(matrix)
  
  #calculate a confidence interval from a t distribution
  interval = qt(.975, df=n-1)*st_dev/sqrt(n)
  
  return(avg-st_dev)
}

Cluster1 = aggregate(data,list(sub_grp4),mean)
Cluster1_upper = aggregate(data,list(sub_grp4),CI_upper)
Cluster1_lower = aggregate(data,list(sub_grp4),CI_lower)
Cluster2=aggregate(data,list(sub_grp5),mean)
Cluster1
Cluster2

#plot these clusters
Cluster1=as.matrix(Cluster1)
Cluster1_upper = as.matrix(Cluster1_upper)
Cluster1_lower = as.matrix(Cluster1_lower)

Cluster1=Cluster1[,-1]
Cluster1_upper = Cluster1_upper[,-1]
Cluster1_lower = Cluster1_lower[,-1]


columns = c(1,4)
test_center = Cluster1[,columns]
test_lower = Cluster1_lower[,columns]
test_upper = Cluster1_upper[,columns]



joey=barplot(test_center, col=colors()[c(125,258,365)] , border="white", font.axis=2, 
        beside=T, legend=rownames(Cluster1), xlab="group", font.lab=2)

arrows(joey, test_lower, joey,
       test_upper, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

  
  






#add in cluster and freq column to above matrix
#a3 = aggregate(AVGdata,list(sub_grp5),median)
#data.frame(Cluster=a3[,1],Freq=as.vector(table(sub_grp5)),a3[,-1])


#let's see how our decisive tree compares to K-means clustering/mediods partioning methods
library(cluster)
PAMa=pam(data_scaled,3)
names(PAMa)

#compare the results to hclust 
table(sub_grp4,PAMa$clustering)
table(sub_grp5,PAMa$clustering)

#solutions seem to agree except for a couple 

#feature with PAM known as the silhouette plot

#.71-1 strong structure has been found
#.5-.7 a reasonable structure has been found
#.26-.5 The structure is weak and could be artificial
# No substantial structure has been found
plot(PAMa)

#silhouette plot for heirarchical cluster analysis
plot(silhouette(cutree(HC4,4),d))
plot(silhouette(cutree(HC5,4),d))





































