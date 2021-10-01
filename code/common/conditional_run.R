check.files.for.rerun.needed=function(files, data.file) {
  if( ! file.exists(data.file))
    return( T )
  if( ! all(file.exists(files) ) )
    return( F )
  newest.t = file.info(files)[,"mtime"] %>%  max
  data.file.t = file.info( data.file)[1,"mtime"]
  return ( data.file.t < newest.t )
}