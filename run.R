# set up working directory (one level above shiny subfolder)

dir = "/Users/hanlinshang/Dropbox/Todos/Taylor_Law_JHMD/TL_JMD/code/"
setwd(dir)

# install R package

install.packages("shiny")

# require R package

require(shiny)

# run shiny app

runApp("shiny")
