# install R packages

install.packages(c("demography", "rainbow", "dplyr"))

# load R packages

require(demography)
require(rainbow)
require(dplyr)

# set working directory

dir = "/Users/hanlinshang/Dropbox/Todos/Taylor_Law_JHMD/TL_JMD/code/"
setwd(paste(dir, "shiny", sep = ""))

df0 = read.csv("asmr.csv", stringsAsFactors = FALSE)

shinyServer(function(input, output)
{
	# df filtered
	dff <- reactive({
		df1 = df0 %>% filter(Code == input$country)
		if(input$sex == "Female")
		{
      dd1 = demogdata(data = matrix(df1$MortFemale, nrow = 100+1),
						pop = matrix(df1$ExpoFemale, nrow = 100+1),
						ages = 0:100, years = unique(df1$Year),
						type = "mortality",
						label = input$country, name = "Female")

      dd1_fts_rate = log(matrix(df1$SmoothmortFemale, nrow = 100+1))
    }
		if(input$sex == "Male")
		{
      dd1 = demogdata(data = matrix(df1$MortMale, nrow = 100+1),
						pop = matrix(df1$ExpoMale, nrow = 100+1),
						ages = 0:100, years = unique(df1$Year),
						type = "mortality",
						label = input$country, name = "Male")

      dd1_fts_rate = log(matrix(df1$SmoothmortMale, nrow = 100+1))
		}
		if(input$sex == "Total")
		{
      dd1 = demogdata(data = matrix(df1$MortTotal, nrow = 100+1),
						pop = matrix(df1$ExpoTotal, nrow = 100+1),
						ages = 0:100, years = unique(df1$Year),
						type = "mortality",
						label = input$country, name = "Total")

      dd1_fts_rate = log(matrix(df1$SmoothmortTotal, nrow = 100+1))
		}
		colnames(dd1_fts_rate) = 1975:2018
		rownames(dd1_fts_rate) = 0:100
		dd1_fts = fts(0:100, dd1_fts_rate, xname = "Age")
		list(df1 = df1, dd1 = dd1, dd1_fts = dd1_fts)
	})

	output$plot1 <- renderPlot({
		dff1 <- dff()
		par(mar=c(2,2,1,0.1))
    if(input$smooth2 == FALSE)
		  plot(e0(dff1$dd1), ylab = "", main = "")
    if(input$smooth2 == TRUE)
      plot(e0(smooth.demogdata(dff1$dd1)), ylab="", main="")
	})

	output$plot2 <- renderPlot({
		dff1 <- dff()
		par(mar=c(2,2,1,0.1))
		if(input$smooth2 == FALSE)
			 plot(dff1$dd1, ylab = "", main = "", plotlegend = TRUE, legendpos = "topleft")
		if(input$smooth2 == TRUE)
			 plot(smooth.demogdata(dff1$dd1), ylab = "", main = "", plotlegend = TRUE, legendpos = "topleft")
	})

	output$plot3 <- renderPlot({
	  dff1 <- dff()
	  par(mar=c(2,2,1,0.1))
	  fboxplot(dff1$dd1_fts, type = "bag", factor = as.numeric(input$factor), legendpos = "topleft", xlab="Age", ylab = "Log mortality rate")
	})

  output$plot4 <- renderPlot({
    dff1 <- dff()
    par(mar=c(2,2,1,0.1))
    fboxplot(dff1$dd1_fts, type = "hdr", alpha = c(as.numeric(input$alpha[1]), as.numeric(input$alpha[2])), legendpos = "topleft", xlab="Age", ylab = "Log mortality rate")
  })

  # output$plot5 <- renderPlot({
  #  dff1 <- dff()
  #  par(mar=c(2,2,1,0.1))
  #  fbplot(dff1$dd1_fts$y)
  # })

	output$plot_t1 <- renderText({
		dff1 <- dff()
		paste0(input$sex, " Life Expectancy ", range(dff1$dd1$year)[1], "-", range(dff1$dd1$year)[2])
	})

	output$plot_t2 <- renderText({
		dff1 = dff()
		paste0(input$sex, " ASMR ", range(dff1$dd1$year)[1], "-", range(dff1$dd1$year)[2])
	})

	output$plot_t3 <- renderText({
	  dff1 = dff()
	  paste0("Functional Bagplot of Hyndman and Shang (2010)")
	})

	output$plot_t4 <- renderText({
	  dff1 = dff()
	  paste0("Functional HDR Boxplot of Hyndman and Shang (2010)")
	})

  # output$plot_t5 <- renderText({
  #  dff1 = dff()
  #  paste0("Functional boxplot of Sun and Genton (2011)")
  # })
})
