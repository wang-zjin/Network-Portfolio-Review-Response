# Install and load packages
#install.packages("IntroCompFinR", repos="http://R-Forge.R-project.org")
libraries = c("matrixcalc", "MASS", "PortfolioAnalytics", 
              "PerformanceAnalytics", "zoo","xts","timeSeries",
              "plotly", "RiskPortfolios", "devtools", "PMwR","Jmisc",
              "igraph","readxl","quadprog",
              "viridis","hrbrthemes",
               "FRAPO",  "R.utils", "ade4", "grDevices", "foreach", 
               "doParallel","StepwiseTest", "stringr",  "DescTools","DT",
              "ggpubr","knitr", "tidyr", "plotly", "rmarkdown", "gridExtra", "reticulate",
              "fPortfolio", 
              "xtable", "DEoptim", "ROI", "ROI.plugin.quadprog", "ROI.plugin.glpk",
              "nlshrink", "psych", 
              "quantmod", 
              "tcltk2","limSolve",
              "MTS","Matrix","pracma","glasso","MST", "IntroCompFinR",
              "RColorBrewer", "broom", "egg", "nloptr", "fAssets",
              "magick") 
              
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})

lapply(libraries, library, quietly = TRUE, character.only = TRUE)

