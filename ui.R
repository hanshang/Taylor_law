# install R packages

install.packages(c("shiny", "demography", "rainbow", "fda"))

# load R packages

library(shiny)
library(demography)
library(rainbow)
library(fda)

# set working directory

dir = "/Users/hanlinshang/Dropbox/Todos/Taylor_Law_JHMD/TL_JMD/code/"
setwd(paste(dir, "shiny", sep = ""))

load("label.RData")

shinyUI(fluidPage(
  # Application title
  titlePanel("Japanese Mortility Database"),
  sidebarLayout(position = "left",
    sidebarPanel(
    	h4("Functional time series plot"),
  		br(),
		fluidRow(
		    column(12, selectInput(inputId = "country", label = "Select Region", choices = cn0, width="100%" )),
		    column(12, selectInput(inputId = "sex", label = "Select Sex", choices = c("Female" = "Female", "Male" = "Male", "Total" = "Total"), width="100%" )),
		    column(12, checkboxInput(inputId = "smooth2", label = "Smooth schedule of Age-specific mortality rate"))
		),
		br(),
		h4("Functional outlier detection"),
		br(),
		fluidRow(
		    column(12, selectInput(inputId = "factor", label = "Functional bagplot outlier percent",
		                        choices = list("20% of observations" = 1.19,
		                                       "10% of observations" = 1.81,
		                                       "5% of observations" = 2.33,
		                                       "1% of observations" = 3.29,
		                                       "0.5% of observations" = 3.64), selected = 3)),
		    column(12,
		           sliderInput(inputId = "alpha", label = "Functional highest density regions",
		                       min = 0, max = 0.5, value = c(0.05, 0.5))
		    )
		)
	  ),
    mainPanel(
    	fluidRow(
		    column(6, strong(textOutput("plot_t1")),  plotOutput("plot1")),
		    column(6, strong(textOutput("plot_t2")),  plotOutput("plot2"))
	  	),
	  	fluidRow(
	  	  column(6, strong(textOutput("plot_t3")),  plotOutput("plot3")),
	  	  column(6, strong(textOutput("plot_t4")),  plotOutput("plot4"))
	  	  # column(4, strong(textOutput("plot_t5")),  plotOutput("plot5"))
	  	)
    )
  )
))
