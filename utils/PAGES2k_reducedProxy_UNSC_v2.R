# Written by Wan He, Dec 2016
# modified by Julien Emile-Geay, Feb 21 2017
# - eliminated "eigenvalue" selection bc results were garbage
# - changed PCA to PCR
# - only fill in composite up from tMin onwards
rm(list=ls())
library(dplyr)
library(fda)
library(ggplot2)
library(reshape2)
library(cowplot)
library(grid)
library(glmnet)
library(parallel)
library(matrixStats)
library(dr)
#library(edrGraphicalTools)
library(stringr)
library(spls)
library(superpc)
library(tidyr)
library(ggmap)
#library(maptools)
library(maps)
library(config)
library(yaml)
library(here)

cfg <- read_yaml(here("config.yml"))

# TO-DOS/ISSUES:
#   1. Breaks down if t1 != 1, make dynamic.
#   2. Fix grabbing desired proxies by type

set.seed(2018)
# specify options to treat the data
t1 <- cfg$partition_years$t1
t2 <- cfg$partition_years$t2
t3 <- cfg$partition_years$t3
tStart = t1 #define start year (remember: the Common Era does not have a year 0).
tEnd   = t3 #define end year for the analysis
tce  = tStart:tEnd
nce = length(tce)

#load the temperature data
temp <- read.csv(here(cfg$folder_paths$instr_temp_path))
colnames(temp) <- c("year","T","l95","u95")
print(temp)

#load the proxy data
proxydata =read.csv(here('data','PAGES2K_proxy_matrix_1-2000.csv'), header=TRUE, sep=",")
#View(proxydata)

#load metadata
metadataproxy =read.csv(here('data','PAGES2K_proxy_metadata_1-2000.csv'), header=TRUE, sep=",")
#View(metadataproxy)

### GRAB THE DESIRED PROXY TYPE FOR RP
# "ALL", "lake", "speleothem", "ice", "borehole", "tree", "documents", "hybrid", "marine", "coral", "bivalve"
if (cfg$ptype=='ALL'){
  NULL # nothing to do
}else{
  desired_proxies <- metadataproxy[sub("\\..*$", "", metadataproxy$ptype) == cfg$ptype, ]
  proxydata <- proxydata[, colnames(proxydata) %in% desired_proxies$pid, drop = FALSE]
}

############################

t=proxydata[,1]
proxy=as.data.frame(proxydata[,-1])[which(t%in%tce),]
ny=dim(proxy)[1]
np=dim(proxy)[2]
#Get the years of the first available data for each proxy
yearMin= matrix(1, 1, np) 
for (i in 1:np) {
  yearMin[i]=min(which(is.finite(proxy[,i])))
}


#####HISTOGRAM FIRST YEAR#############
datahist <- data.frame(Year = c(yearMin))
histogram_first <- ggplot(data = datahist)+
  geom_histogram(mapping = aes(x = Year),bins = 10,fill='white',color='black')+
  ylab('')+
  theme_bw()+
  theme(axis.title.x = element_text(size=18),axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))
histogram_first

# define analysis options
RP_style = "LASSO" #cfg$rp_method  #possible choices: "PCR", "LASSO", "SIR", "SPLS" #### EDIT THIS YML!!!
## DEFINE FUNCTIONS USED LATER
#define function %!in% which is same as ~ismember
'%!in%' <- function(x,y)!('%in%'(y,x))

#define rfind fun ction as similar to find in matlab
rfind <- function(x)seq(along=x)[as.logical(x)] 

# scale the proxy records to unit variance and zero mean
proxy_scaled = scale(proxy)
#View(proxy_scaled)

# NEW
proxy_filled <- apply(proxy_scaled, 2, function(x) {
  seen_first <- cumsum(!is.na(x)) > 0
  x[is.na(x) & seen_first] <- 0
  x
})
rownames(proxy_filled) <- rownames(proxy_scaled)
colnames(proxy_filled) <- colnames(proxy_scaled)
#View(proxy_filled)

# count the number of available datapoints in each record
navl = matrix(1, 1, np) 
for (i in 1:np){  
  navl[i] = sum(is.finite(proxy_filled[,i]))
}

##  MAKE A TIME-DEPENDENT COMPOSITE
chunk = 250 # define  segment length
ns = ny/chunk  # number of segments
calib = c(t2:t3)#Define calibration year
RP = matrix(NA,ny,ns)#Initialise the reduced proxies
nprox= matrix(0,1,ns)#Initialise nprox that is later used to store number of proxies used in each segmentation
temp_lastc <- temp[which(temp[,1]%in%calib),2] #last century temperature

#evalues=matrix(1,1,npcs)#Initialise a matrix to store the eigenvalues of PCs
coefflist <- list()
npcs <- c(11,13,9,11,7,12,10,6) #first PC whose adjusted R2 > 70%
limitup <- c(1,2,2,2,2,1,1,2) ##Marginal Dimension Tests (SIR)
nslicesv <- c(5,7)
proport <- c(1,1,1,1,1,0.35,0.6,0.75)
pZerov <- c(30,30,30,30,30,30,30,30)
#tsholdv <- c(0.53,1.16,0.36,2.62,1.16)
tsholdv <- NULL
npcs2 <- c(3,3,3,3,3,3,3,3)

for (k in 1:ns){   # loop over segments
  show(k)
  tMin = tStart+(k-1)*chunk
  segProxies = rfind(yearMin <=tMin) 
  timeSpan  = rfind(tce>=tMin)
  nt = length(timeSpan)
  one=matrix(1,nt,1)#Initialise a ny*1 one matrix that can be later appended to PC matrix
  proxy_finite= proxy_filled[timeSpan,segProxies]
  cleanNANs <- apply(X = proxy_finite,MARGIN = 2,function(x) sum(is.na(x))/dim(proxy_finite)[1])<0.05
  proxy_finite <- proxy_finite[,cleanNANs]
  time_calib <- rfind(timeSpan %in% calib)
  proxy_finite_calib= proxy_finite[time_calib,]
  cleanNANs_calib <- apply(X = proxy_finite_calib,MARGIN = 2,function(x) sum(is.na(x))/dim(proxy_finite_calib)[1])<0.05
  proxy_finite <- proxy_finite[,cleanNANs_calib]
  proxy_finite[!is.finite(proxy_finite)]=0
  nprox[k] = dim(proxy_finite)[2]
  
  if (RP_style=="PCR") { 
      pca=prcomp(proxy_finite)
      a <- lm(temp_lastc~pca$x[rfind(timeSpan %in% calib),1:npcs[k]])
      coeff <- as.numeric(a$coefficients)
      RP[timeSpan,k]   = cbind(one,pca$x[,1:npcs[k]])%*%coeff
  }
  if (RP_style=="sPCR") {
    datapc <- list(x=t(proxy_finite[rfind(timeSpan %in% calib),]),y=temp_lastc)
    trainobj <- superpc.train(datapc,type='regression')
    cvobj <- superpc.cv(trainobj,datapc)
    superpc.plotcv(cvobj)
    tsholdv[k] <- cvobj$thresholds[which.max(cvobj$scor[1,])]
    datatest <- list(x=t(proxy_finite))
    fitobj <- superpc.predict(trainobj, datapc, datatest, threshold=tsholdv[k], 
                              n.components=npcs2[k], prediction.type="continuous")
    datapred <- data.frame(temp_lastc,fitobj$v.pred[rfind(timeSpan %in% calib),])
    modelopred <- lm(temp_lastc~.,data = datapred)
    RP[timeSpan,k]   <- predict(modelopred,newdata = data.frame(fitobj$v.pred))
  }
  
  if (RP_style == "LASSO") {
    # Performs a lasso regularization technique to avoid the use of PCA in multiple regression
    Xmatrix <- proxy_finite[rfind(timeSpan %in% calib), ]
    cv.lasso <- cv.glmnet(x = Xmatrix, y = temp_lastc, alpha = 1, parallel = TRUE)
    coeff <- coef(cv.lasso, s = "lambda.min")
    coefflist[[k]] <- coeff
    RP[timeSpan, k] <- predict(cv.lasso, newx = proxy_finite, s = "lambda.min")
  }
  if (RP_style=="SIR"){
    ### TYLER: NEED TO FIGURE OUT HOW TO INSTALL edrGraphicalTools!!
    #Performs a Sliced Inverse Regression (SIR)
    Xmatrix <- proxy_finite[rfind(timeSpan %in% calib),]
    show(dim(Xmatrix)[2])
    if(k >= 6){
      selectModelCSS <- edrSelec(Y = temp_lastc,X = Xmatrix,H = 8,K = 1,
                               method = 'CSS',pZero = nprox[5],NZero = 10000,zeta = 0.1)
      show(max(selectModelCSS$vectSqCor))
      #indexselect <- selectModelCSS$scoreVar>quantile(selectModelCSS$scoreVar,proport[k])
      indexselect <- sort(selectModelCSS$scoreVar,decreasing = T,index.return=T)$ix<=nprox[5]
      show(sum(indexselect))
      Xmatrix <- Xmatrix[,indexselect]
    }
    SIRdr <- dr(temp_lastc ~ Xmatrix,nslices=8)
    namesafterSIR <- names(SIRdr$evectors[,1])%>%str_replace('Xmatrix','')
    Xmatrix <- Xmatrix[,namesafterSIR]
    xtemp <- as.matrix(Xmatrix %*% SIRdr$evectors[,1:(limitup[k])]) 
    modeltemp <- lm(temp_lastc~xtemp)
    xtemptot <- as.matrix(proxy_finite[,namesafterSIR] %*% SIRdr$evectors[,1:(limitup[k])])
    colnames(xtemptot) <- names(modeltemp$coefficients)[(1:(limitup[k]))+1]
    RP[timeSpan,k] <- cbind(one,xtemptot)%*%modeltemp$coefficients
  }
    
  if (RP_style=="SPLS"){
    Xmatrix <- proxy_finite[rfind(timeSpan %in% calib),]
    show(dim(Xmatrix)[2])
    Kint <- seq(1,min(dim(Xmatrix)[2],floor(0.9*dim(Xmatrix)[1]-1)))
    #Kint <- seq(1,70)
    cvspls <- cv.spls(x = Xmatrix,y = temp_lastc,K = Kint,eta = seq(0.1,0.9,0.1))
    SPLSmodel <- spls(x = Xmatrix,y = temp_lastc,K = cvspls$K.opt,eta = cvspls$eta.opt)
    RP[timeSpan,k] <- predict.spls(SPLSmodel,newx = proxy_finite)
  }
}


colnames(RP) <- paste0('RP',seq(1,8))
RP <- data.frame(Year=t1:t3, RP)

RPn <- RP %>% gather(key = 'RPNumber',value = 'Value',RP1:RP8)

plotcombined <- ggplot(data = RPn,mapping = aes(x = Year,y = Value))+
  geom_line(mapping = aes(color=RPNumber))+
  scale_color_brewer(palette = 4,type = 'div')+
  ggtitle(RP_style)+
  theme_bw()
plotcombined

#save(RPn,file=paste0('Plot_RP',RP_style,'.RData'))
#dataset <- 'All'
#save(RP,file = paste0('./results/RPs/RP_new_',dataset,'_',RP_style,'.RData'))


###SINGLE RP SERIES 
RPind <- data.frame(
  Year = RP$Year,
  RP1 = apply(RP[, -1], 1, function(x) mean(x, na.rm = TRUE))
)
write.csv(RPind, here("data", "RPind.csv"), row.names = FALSE)
plot(RPind$Year, RPind$RP1, type='l')
head(RPind)



