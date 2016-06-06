


if('.mpga' %in% search()) detach('.mpga')
if(exists('.mpga')) rm(.mpga)
if(exists('.tempVars')) rm(.tempVars)
.mpga = new.env()
.tempVars = new.env()
####
assign('pathIntrospection', function(){
  frameFiles = Filter(Negate(is.null),
                      lapply(sys.frames(),function(x)x$ofile))
  return(dirname(frameFiles[[length(frameFiles)]]))},
  envir = .tempVars
)
.tempVars$loadpath = .tempVars$pathIntrospection()
####
.tempVars$templs = list.files(.tempVars$loadpath)
.tempVars$rfiles = grep("(?<!init)[.][rR]",.tempVars$templs,perl = T,value = T)
.tempVars$path2rfiles = sapply(.tempVars$rfiles,function(x){
  ifelse(.tempVars$loadpath!='',file.path(.tempVars$loadpath,x),x)})

sapply(.tempVars$path2rfiles,function(x)sys.source(x,.mpga))
# compiled version of functions
# .tempVars$funname = sub("[.][Rr]$","",.tempVars$rfiles)
# m_ply(.tempVars$funname,function(x)assign(x,value = cmpfun(get(x,envir = .garf)),envir= .garf))
attach(.mpga,pos=2,warn.conflicts=FALSE)
rm(.tempVars)

