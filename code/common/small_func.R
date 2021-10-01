library(lubridate)
da = function(x, year=2020){
  if( year< 1000) year=year+2000
        m = as.integer( floor(x) )
        d = as.integer(round( (x-m)*100 ))
        ymd( sprintf("%d-%02d-%02d",year,m,d))
}

daS = function( from, to=c(), by=c(),len=c(), year=2020,...) {
    x = c(from, to)
    x=sort(x)
    if( is.numeric(x) ) {
        if( x[1]<0) {
            dx = c( da(x[2], year)+x[1], da(x[2],year) )
        } else if( length(x)>1  & x[2] == floor(x[2]) ) {
            dx = c( da(x[1], year), da(x[1])+x[2], year )
        } else {
            dx = da(x, year)
        }
    } else {
        dx = x
    }
    if( !is.null(by) & !is.null(len))
      res = seq( dx[1], by=by, len=len)
    else if( !is.null(by))
      res= seq( dx[1], dx[2], by=by )
    else if( !is.null( len)) {
      if( length(x)> 1) 
         res = seq( dx[1], dx[2], len=len )
      else {
         res = seq( dx[1], dx[1]+len-1, by=1)
      }
    } else  
     res = seq( from=dx[1], to=dx[2], by=1 )
    sort(res)
}

save.RDS.dir = function( list, dir, overwrite=F) {
  dir.create(dir, showWarnings=F)
  for(f in list) {
    file = paste0(dir,"/",f,".Rds")
    if( overwrite || !file.exists(file))
      saveRDS(get(f),file=file)
  }
}




my.restore.session = function (file = ".RSession", ...) 
{
    cat("Loading all data...\n")
    load(file, envir = .GlobalEnv)
    cat("Loading packages...\n")
    sapply(rev(.save.session.packages), library, character.only = TRUE)
    cat("Restoring search path...\n")
    pad <- function(x, n) c(rep(NA, n - length(x)), x)
    current.search <- search()[-1]
    saved.search <- .save.session.search[-1]
    identical <- saved.search %in% current.search
    for (i in saved.search[!identical]) {
        if (charmatch("file:", i, nomatch = FALSE)) 
            attach(sub("file:", "", i))
        else if (charmatch("package:", i, nomatch = FALSE)) 
            stop(paste("Somehow we missed loading package", i))
        else {
            do.call("attach", list(as.name(i)))
        }
    }
    rm(list = c(".save.session.packages", ".save.session.search"), 
        envir = .GlobalEnv)
    cat("Done.\n")
}
