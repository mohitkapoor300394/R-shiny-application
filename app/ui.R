ui <- fluidPage(
  pageWithSidebar(
    headerPanel('CGM Data Visualization'),
    sidebarPanel(
      fileInput('cgmdata', 'Upload File Containing CGM Data',multiple=FALSE),
      selectInput('device','Device Types',c('ipro2','libre','medtronic','any other devices')),
      numericInput('maxhypo', 'Hypoglycaemic (mmol/L)<=',4.0),
      numericInput('minrate', 'Abnormal Rates (mmol/L/hr) >=',1.00),
      numericInput('Hprior', 'Number of hours since last hypoglycaemic',2),
      numericInput('mintopBG', 'Minimum Top BG during episode',10),
      numericInput('minlen', 'Minimum Segment Length (hr)',1),
      submitButton(text = "Submit"),
      htmlOutput('UserGuide')
    ),
    mainPanel(
      # plot SOI (segment-of-interest) with date as x-axis
      #plotOutput('PlotSeg'),
      # show the episodes in date format
      dataTableOutput('Eps'),
      # show data with (episodes marked)
      plotOutput('DataPlot'),
      # show bubble plot of episodes where bubble size = length    
      plotlyOutput('BubbleEps')
      #display histogram of time-freq of SOI
      #plotOutput('EpsFreq')
    )
  )
)
