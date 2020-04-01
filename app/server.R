addResourcePath('static', 'static')

source('static/pickSeg.R')
server <- function(input, output,session) {
  # read data in, pre-proc, summarize using HMM and pick SOI
  summSeg <- reactive({
    if(is.null(input$cgmdata))  return(NULL) 
    
    filename <- input$cgmdata$datapath
    # branching here, depending on device
    
    if (input$device=='ipro2') 
   {
      out <- summSeg.cgm(file=filename)
      int.min <- round(as.numeric(names(sort(table(diff(out$time))))[1])*60)
      out2<- pickSeg.cgm(out,crit.lower=list(MdBG=input$minrate/(60/int.min),maxBG=input$mintopBG,width=60*input$minlen/int.min),crit.upper=list(MdBG=Inf,maxBG=Inf,width=Inf),hypo=input$maxhypo,since.hypo=input$Hprior)
    }
   else if (input$device=='libre') 
      
    {
      out <- summSeg.cgm2(file=filename)
      int.min <- round(as.numeric(names(sort(table(diff(out$time))))[1])*60)
      out2<- pickSeg.cgm(out,crit.lower=list(MdBG=input$minrate/(60/int.min),maxBG=input$mintopBG,width=60*input$minlen/int.min),crit.upper=list(MdBG=Inf,maxBG=Inf,width=Inf),hypo=input$maxhypo,since.hypo=input$Hprior)
    }
    
    #prepare output of this section
    summSeq.out <- list()
    summSeq.out[['out']] <- out ; summSeq.out[['out2']] <- out2
    summSeq.out
  })
  
  
  showEps <- reactive({
    # extract time
    if(is.null(input$cgmdata))  return(NULL) 
    seg <- summSeg() ; out <- seg$out ; out2 <- seg$out2
    # change time to continuous time
    day <- (as.numeric(out$date)-min(as.numeric(out$date))) + 1
    out$time <- (day-1)*24 + out$time
    time.eps <- data.frame(date=out$date[out2$stat$start[out2$stat$from.hypo==0]], start=out$time[out2$stat$start[out2$stat$from.hypo==0]] %% 24, 
                           length= out$time[out2$stat$end[out2$stat$from.hypo==0]]-out$time[out2$stat$start[out2$stat$from.hypo==0]])
    time.eps <- time.eps[!(apply(time.eps, 1, function(y) any(y == 0))),]
    #format time and date nicely
    hour <- trunc(time.eps$start) ; min <-  round( (time.eps$start-trunc(time.eps$start))*60)
    min  <- ifelse(min<10,paste0('0',min),min)
    time.eps$start <- paste0(hour,':',min)
    time.eps$length<- round(time.eps$length*60)
 
    colnames(time.eps) <- c('date','start(24hour_clock)','length(min)')
    time.eps <- format(time.eps,justify='right')   
    time.eps

    
  })
  
  histSeg <- reactive({
    seg <- summSeg() ; out <- seg$out ; out2 <- seg$out2  
    plot(table(out2$time[out2$picked]),ylab='Frequency',xlab='24-hour clock')
  })
  
  # plot data
  plotData <- reactive({
    if(is.null(input$cgmdata))  return(NULL) 
    # extract time
    seg <- summSeg() ; out <- seg$out ; out2 <- seg$out2
    col <- (as.numeric(out$date)-min(as.numeric(out$date))) + 1
    time.append <- (col-1)*24+out$time   
    plot((col-1)*24+out$time,out$data,col=as.numeric(as.factor(out$date)),axes=F,xlim=c(0,max(time.append)),cex=0.3, ylab='BG (mmol/L)', xlab='24 hour time clock')
    axis(1,at = seq(0,max(time.append),6),label = seq(0,max(time.append),6) %% 24)
    axis(2)
    legend('topleft',legend=unique(out$date),col=1:length(unique(out$date)),lty=1,bty='n')
    abline(h=input$maxhypo)
  })
  
  # bubble plot episodes
  bubbleEps <- reactive({
    if(is.null(input$cgmdata))  return(NULL) 
    # extract time
    seg <- summSeg() ; out <- seg$out ; out2 <- seg$out2
    require(ggplot2)
    require(plotly)
    day <- (as.numeric(out$date)-min(as.numeric(out$date))) + 1
    # change time to continuous time
    out$time <- (day-1)*24 + out$time
    out$date <- as.factor(out$date)
    col.pallete <- c('black','red','green','blue','cyan','magenta','grey','orange')
    time.eps <- data.frame(date=out$date[out2$stat$start[out2$stat$from.hypo==0]], start=out$time[out2$stat$start[out2$stat$from.hypo==0]] %% 24, 
                           length= as.numeric(out$time[out2$stat$end[out2$stat$from.hypo==0]]-out$time[out2$stat$start[out2$stat$from.hypo==0]]),
                           day = day[out2$stat$start[out2$stat$from.hypo==0]], maxBG = out2$stat$maxBG[out2$stat$from.hypo==0])
    time.eps$date <- factor(time.eps$date,levels=unique(out$date))
    obj <- ggplot(aes(start,maxBG,size=length,col=date),data=time.eps)+geom_point() + theme_bw() + labs(x='Starting Time of Eps',y='Max BG (mmol/L) during eps',col='Date') +
      lims(x=c(0,24),y=c(0,max(time.eps$maxBG)+0.5)) 
    ggplotly(obj)
  })
  
  #output$EpsFreq<- renderPlot(histSeg())
  output$DataPlot<- renderPlot(plotData())
  output$BubbleEps <- renderPlotly(bubbleEps())
  output$Eps   <- renderDataTable(showEps())
  
  url <- a("User Guide", href="UserGuide.pdf",target="_blank")
  output$UserGuide <- renderUI({ tagList("Help:", url) })
  
}