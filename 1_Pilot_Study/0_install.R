install.packages("devtools")
install.packages("mvtnorm")
install.packages("dequer")
install.packages("knitr")
install.packages("rmarkdown")

library(devtools)

install_github("katzfuss-group/GPvecchia")

# sessionInfo()
# Sys.getenv()
# 
# Sys.setenv(PATH = "C:/Rtools/bin")
# Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# 
# options(buildtools.check = function(action) TRUE )
# Sys.setenv("USE_CXX11" = "yes")