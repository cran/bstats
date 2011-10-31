.onLoad <- function(lib, pkg){
  packageStartupMessage("bstats 1.0-1 loaded\nCopyright B. Wang 2011")
  assign('.bstatsConnect',NULL,pos=.GlobalEnv) 
}

.onUnload <- function(libpath)
    library.dynam.unload("bstats",  libpath)

.bstatsConnect <- NULL
