summSeg.cgm <- function(file,max.nstate=15) {
  require(HiddenMarkov)
  require(plotly)
  output <- list()
  # read data in
  data <- read.csv(file, skip=11, header=T)
  # get recorded meal time
  meal <- data[data$Meal=='Meal',]
  # get data with recorded glucose sensor
  data.glu <- data[!is.na(data$Sensor.Glucose..mmol.L.),]
  
  # get date
  date <- as.Date(data.glu$Date, "%d/%m/%Y")
  date
  # get hour, minute and second
  timesplit <- unlist(strsplit(as.character(data.glu$Time), split=":"))
  hour <- as.numeric(timesplit[seq(1,length(timesplit),3)])
  mins <- as.numeric(timesplit[seq(2,length(timesplit),3)])
  secs <- as.numeric(timesplit[seq(3,length(timesplit),3)])
  
  uniq.date <- unique(date)
  cont.24hr <- hour + mins/60 + secs/3600
  
  # fit HMM to differentiate states: decreasing, neutral and rising glucose
  # using log BG
  o.x <- log(data.glu$Sensor.Glucose..mmol.L.)
  # smoothed the log BG data (to rid of minor peaks)
  x <- smooth.spline(o.x,cv=TRUE)$y
  # calculate difference 
  d <- diff(x)
  
  # start fitting HMM to find the best model
  BIC <- NULL
  ncomp<- NULL
  conv <- FALSE
  require(doParallel)
  cl <- parallel::makeCluster(detectCores()-1,outfile='out.txt')
  doParallel::registerDoParallel(cl)
  
  nstate <- seq(6,20,2)
  
  temp <- foreach(i = 1:length(nstate), .combine = 'rbind', .packages=c('HiddenMarkov')) %dopar% {  
    Pi <- matrix(1/nstate[i], byrow=TRUE, nrow=nstate[i], ncol=nstate[i])
    delta <- rep(0,nstate[i]); mean.init=runif(nstate[i],min(d),max(d))
    delta[which.min(abs(d[1]-mean.init))] <- 1
    hmm.model <- dthmm(x=na.omit(d), Pi, delta, "norm", list(mean=mean.init, sd=abs(rnorm(nstate[i],sd=0.5))))
    # use above parameter values as initial values 
    y <- tryCatch({BaumWelch(hmm.model,control=bwcontrol(prt=FALSE))}, error=function(e) { return(NA)})
    # number of est par
    k <- nstate[i]^2 + 2*nstate[i] - nstate[i]
    if(!is.na(y)) {
      BIC  <- -2*y$LL + log(summary(y)$n)*k
      ncomp<- nstate[i]
    }
    return(c(BIC,ncomp))
  }
  BIC <- temp[,1]
  ncomp <- temp[,2]
  stopCluster(cl)
  
  # get the segment for HMM model with optimal nstate (with min BIC)
  nstate <- ncomp[which.min(BIC)]
  Pi <- matrix(1/nstate, byrow=TRUE, nrow=nstate, ncol=nstate)
  delta <- rep(0,nstate); mean.init=runif(nstate,min(d),max(d))
  delta[which.min(abs(d[1]-mean.init))] <- 1
  HMM.model <- dthmm(x=na.omit(d), Pi, delta, "norm", list(mean=mean.init, sd=abs(rnorm(nstate,sd=0.5))))
  # use above parameter values as initial values 
  y <- BaumWelch(HMM.model,control=bwcontrol(prt=FALSE))
  state <- Viterbi(y)
  
  # summarize each segments based on their segment characteristics
  desc.stat <- NULL
  dlog.BG   <- c(NA,d)
  d.BG      <- c(NA,diff(exp(o.x)))
  state     <- c(NA,state)
  dlog.dstate<- c(NA,diff(state))
  pos.dlog.dstate <- which(abs(dlog.dstate)>0)
  st <- pos.dlog.dstate[-length(pos.dlog.dstate)]
  en <- pos.dlog.dstate[-1] - 1
  
  
  
  # we can add/remove summary stats here. Summary stats can be used to filter segments later.
  registerDoSEQ()
  for(seg in 1:length(st)){
    seg.data <- exp(o.x[st[seg]:en[seg]])
    seg.dlog <- d[st[seg]:en[seg]]
    seg.d    <- d.BG[st[seg]:en[seg]]
    
    width   <- en[seg]-st[seg] + 1
    mean.BG  <- mean(seg.data,na.rm=T)
    max.BG   <- max(seg.data,na.rm=T)
    start.BG <- max(seg.data[1],na.rm=T)
    sd.BG   <-  sd(seg.data,na.rm=T)
    seg.st   <- st[seg]
    seg.en   <- en[seg]
    
    mean.d   <- mean(seg.d, na.rm=T)
    mean.dlog<- mean(seg.dlog,na.rm=T)
    seg.state <- state[st[seg]]
    desc.stat <- rbind(desc.stat,c(seg.st,seg.en,width,mean.BG,mean.d,mean.dlog,max.BG,start.BG,sd.BG,seg.state))
  }
  
  desc.stat <- data.frame(desc.stat)
  names(desc.stat) <- c('start','end','width','MBG','MdBG','MdlogBG','maxBG','stBG','sdBG','state')
  
  #compile output 
  output <- list()
  output$data <- exp(o.x)
  output$sdata <- exp(x)
  output$date <- date
  output$time <- cont.24hr
  output$hmm <- HMM.model
  output$BIC <- BIC
  output$state <- state
  output$d <- d
  output$stat  <- desc.stat
  # return output 
  saveRDS(output,file='output.rds')
  output
}

summSeg.cgm2 <- function(file,max.nstate=20) {
  require(HiddenMarkov)
  require(lubridate)
  output <- list()
  # read data in
  data <- read_csv(file)
  
  # get data with recorded glucose sensor
  # get data with recorded glucose sensor
  data.glu1 <- data[!is.na(data$`Scan Glucose(mmol/L)` ),]
  
  # get date
  
  data.glu1$time <- format(as.POSIXct (strptime(data.glu1$`Meter Timestamp`,"%d/%m/%Y %H:%M",tz="")) ,format = "%H:%M")
  data.glu1$Dates <- format(as.POSIXct (strptime(data.glu1$`Meter Timestamp`,"%d/%m/%Y %H:%M",tz="")) ,format = "%d/%m/%Y")
  
  date <- as.Date(data.glu1$Dates, "%d/%m/%Y")
  
  timesplit <- unlist(strsplit(as.character(data.glu1$time), split=":"))
  hour <- as.numeric(timesplit[seq(1,length(timesplit),2)])
  mins <- as.numeric(timesplit[seq(2,length(timesplit),2)])
  
  
  uniq.date <- unique(date)
  cont.24hr <- hour + mins/60 
  
  # fit HMM to differentiate states: decreasing, neutral and rising glucose
  # using log BG
  o.x <- log(data.glu1$`Scan Glucose(mmol/L)` )
  # smoothed the log BG data (to rid of minor peaks)
  x <- smooth.spline(o.x,cv=TRUE)$y
  # calculate difference 
  d <- diff(x)
  
  
  # start fitting HMM to find the best model
  BIC <- NULL
  ncomp<- NULL
  conv <- FALSE
  require(doParallel)
  cl <- parallel::makeCluster(detectCores()-1,outfile='outlibre.txt')
  doParallel::registerDoParallel(cl)
  
  nstate <- seq(6,20,2)
  
  temp <- foreach(i = 1:length(nstate), .combine = 'rbind', .packages=c('HiddenMarkov')) %dopar% {  
    Pi <- matrix(1/nstate[i], byrow=TRUE, nrow=nstate[i], ncol=nstate[i])
    delta <- rep(0,nstate[i]); mean.init=runif(nstate[i],min(d),max(d))
    delta[which.min(abs(d[1]-mean.init))] <- 1
    hmm.model <- dthmm(x=na.omit(d), Pi, delta, "norm", list(mean=mean.init, sd=abs(rnorm(nstate[i],sd=0.5))))
    # use above parameter values as initial values 
    y <- tryCatch({BaumWelch(hmm.model,control=bwcontrol(prt=FALSE))}, error=function(e) { return(NA)})
    # number of est par
    k <- nstate[i]^2 + 2*nstate[i] - nstate[i]
    if(!is.na(y)) {
      BIC  <- -2*y$LL + log(summary(y)$n)*k
      ncomp<- nstate[i]
    }
    return(c(BIC,ncomp))
  }
  BIC <- temp[,1]
  ncomp <- temp[,2]
  stopCluster(cl)
  
  
  # get the segment for HMM model with optimal nstate (with min BIC)
  nstate <- ncomp[which.min(BIC)]
  Pi <- matrix(1/nstate, byrow=TRUE, nrow=nstate, ncol=nstate)
  delta <- rep(0,nstate); mean.init=runif(nstate,min(d),max(d))
  delta[which.min(abs(d[1]-mean.init))] <- 1
  HMM.model <- dthmm(x=na.omit(d), Pi, delta, "norm", list(mean=mean.init, sd=abs(rnorm(nstate,sd=0.5))))
  # use above parameter values as initial values 
  y <- BaumWelch(HMM.model,control=bwcontrol(prt=FALSE))
  state <- Viterbi(y)
  
  # summarize each segments based on their segment characteristics
  desc.stat <- NULL
  dlog.BG   <- c(NA,d)
  d.BG      <- c(NA,diff(exp(o.x)))
  state     <- c(NA,state)
  dlog.dstate<- c(NA,diff(state))
  pos.dlog.dstate <- which(abs(dlog.dstate)>0)
  st <- pos.dlog.dstate[-length(pos.dlog.dstate)]
  en <- pos.dlog.dstate[-1] - 1
  
  
  
  # we can add/remove summary stats here. Summary stats can be used to filter segments later.
  registerDoSEQ()
  foreach(seg = 1:length(st)) %dopar% {
    seg.data <- exp(o.x[st[seg]:en[seg]])
    seg.dlog <- d[st[seg]:en[seg]]
    seg.d    <- d.BG[st[seg]:en[seg]]
    
    width   <- en[seg]-st[seg] + 1
    mean.BG  <- mean(seg.data,na.rm=T)
    max.BG   <- max(seg.data,na.rm=T)
    start.BG <- max(seg.data[1],na.rm=T)
    sd.BG   <-  sd(seg.data,na.rm=T)
    seg.st   <- st[seg]
    seg.en   <- en[seg]
    
    mean.d   <- mean(seg.d, na.rm=T)
    mean.dlog<- mean(seg.dlog,na.rm=T)
    seg.state <- state[st[seg]]
    desc.stat <- rbind(desc.stat,c(seg.st,seg.en,width,mean.BG,mean.d,mean.dlog,max.BG,start.BG,sd.BG,seg.state))
  }
  
  desc.stat <- data.frame(desc.stat)
  names(desc.stat) <- c('start','end','width','MBG','MdBG','MdlogBG','maxBG','stBG','sdBG','state')
  
  #compile output 
  output <- list()
  output$data <- exp(o.x)
  output$sdata <- exp(x)
  output$date <- date
  output$time <- cont.24hr
  output$hmm <- HMM.model
  output$BIC <- BIC
  output$state <- state
  output$d <- d
  output$stat  <- desc.stat
  # return output 
  saveRDS(output,file='output.rds')
  output
}




pickSeg.cgm <- function(obj,crit.lower=NULL,crit.upper=NULL,hypo=4,since.hypo=2) {
  # filter first by mean of dBG 
  col.index <- match('MdBG',colnames(obj$stat))
  filtered.stat <- obj$stat
  
  foreach(i = 1:length(col.index))%dopar% 
  {filtered.stat <- filtered.stat[(filtered.stat[,col.index[i]] >= crit.lower[match('MdBG',names(crit.lower))]) & (filtered.stat[,col.index[i]] <= crit.upper[match('MdBG',names(crit.upper))]),]}
  
  # get interval between measurements (in mins)
  int.min <- as.numeric(names(sort(table(diff(obj$time))))[1])*60
  # min.obs that cover 15-min period
  min.obs       <- round(15/int.min) ; min.obs <- ifelse(min.obs<1,1,min.obs)
  # remove segments that start from hypoglycaemic
  from.hypo <- rep(0,nrow(filtered.stat))
  foreach(i = 1:nrow(filtered.stat)) %dopar% {
    # go back to round(60/int.min)*since.hypo before start of segment
    st.hypo  <- filtered.stat$start[i] - round(60/int.min)*since.hypo ; st.hypo <- ifelse(st.hypo<0,1,st.hypo)
    from.hypo[i] <- as.logical(sum(obj$data[st.hypo:filtered.stat$start[i]]< hypo, na.rm=T) >= min.obs)
  } 
  filtered.stat <- filtered.stat[!from.hypo,]
  
  # reduce regions into non-overlapping regions
  library(IRanges)
  filtered.reg <-  IRanges(filtered.stat$start,filtered.stat$end)
  # if two regions are less than 1hour apart, MERGE THEM
  reduced.reg  <- reduce(filtered.reg,min.gapwidth=round(60/int.min))
  st  <- slot(reduced.reg,'start')
  en  <- st + slot(reduced.reg,'width') - 1
  d.BG<- c(NA,diff(obj$data))
  picked <- rep(FALSE,length(obj$data))
  desc.stat <- NULL
  foreach(seg = 1:length(st)) %dopar% {
    seg.data <- obj$data[st[seg]:en[seg]]
    seg.dlog <- obj$d[st[seg]:en[seg]]
    seg.d    <- d.BG[st[seg]:en[seg]]
    
    width   <- en[seg]-st[seg] + 1
    mean.BG  <- mean(seg.data,na.rm=T)
    max.BG   <- max(seg.data,na.rm=T)
    start.BG <- max(seg.data[1],na.rm=T)
    sd.BG   <-  sd(seg.data,na.rm=T)
    seg.st   <- st[seg]
    seg.en   <- en[seg]
    
    mean.d   <- mean(seg.d, na.rm=T)
    mean.dlog<- mean(seg.dlog,na.rm=T)
    # go back to round(60/int.min)*since.hypo before start of segment
    st.hypo  <- seg.st - round(60/int.min)*since.hypo ; st.hypo <- ifelse(st.hypo<0,1,st.hypo)
    # check whether there are at least cumulative 15-min with cgm < hypo threshold
    min.obs       <- round(15/int.min) ; min.obs <- ifelse(min.obs<1,1,min.obs)
    from.hypo     <- sum(obj$data[st.hypo:seg.st]< hypo, na.rm=T) >= min.obs
    desc.stat <- rbind(desc.stat,c(seg.st,seg.en,width,mean.BG,mean.d,mean.dlog,max.BG,start.BG,sd.BG,from.hypo))
    picked[seg.st : seg.en] <- !from.hypo
  }
  desc.stat <- data.frame(desc.stat)
  names(desc.stat) <- c('start','end','width','MBG','MdBG','MdlogBG','maxBG','stBG','sdBG','from.hypo')
  
  # now, apply the remaining filter on the newly merged regions
  crit.lower.oth <- crit.lower[names(crit.lower)!='MdBG']
  crit.upper.oth <- crit.upper[names(crit.upper)!='MdBG']
  col.index <- match(names(crit.lower.oth),colnames(desc.stat))
  foreach(i = 1:length(col.index)) %dopar% {
    desc.stat <- desc.stat[(desc.stat[,col.index[i]] >= crit.lower.oth[i]) & (desc.stat[,col.index[i]] <= crit.upper.oth[i]),]
  }
  # update obj
  obj$stat <- desc.stat
  obj$picked <- picked
  saveRDS(obj,file='output2.rds')
  obj
}
