.First.lib <-function(lib,pkg)
{
ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
     ver <- as.character(ver)	

library.dynam("adlift",pkg,lib)

cat("adlift", ver, "loaded\n")
}