#' @title parsing KEGG xmls ready for pathway enrichment analysis
#' @usage makeSPIAdata(kgml.path=NULL,organism="mmu",out.path=".")
#' @description Reference: Tarca A. L. et al. (2009) A novel signaling pathway impact analysis. Bioinformatics, 25, 75–82.
#' @param kgml.path path to KEGG xml files
#' @param organism three letter organism code from KEGG database
#' @param out.path path to which generated RData file will be placed 
#' @return TRUE on success, FALSE otherwise
#' @author  Original implementation by Adi Laurentiu Tarca <atarca AT med.wayne.edu>, Purvesh Kathri
#'    <purvesh AT cs.wayne.edu> and Sorin Draghici <sorin AT wayne.edu>, source code shown at the bottom.
#'    Modification by Hannan Ai
#' @export

makeSPIAdata<-function(kgml.path=NULL,organism="mmu",out.path="."){
if(is.null(kgml.path)){
  mdir=system.file("extdata/mmuKEGGxml",package="pathwayko")
}else{
  mdir=kgml.path
}


alr<-NULL

L<-list()

info<-NULL
path.nodes<-list()

paths<-list.files(mdir,pattern="*.xml")

if (length(paths)>0){


#new rel
rel<-c("activation",
"compound",
"binding/association",
"expression",
"inhibition",
"activation_phosphorylation",
"phosphorylation",
"inhibition_phosphorylation",
"inhibition_dephosphorylation",
"dissociation",
"dephosphorylation",
"activation_dephosphorylation",
"state change",
"activation_indirect effect",
"inhibition_ubiquination",
"ubiquination",
"expression_indirect effect",
"inhibition_indirect effect",
"repression",
"dissociation_phosphorylation",
"indirect effect_phosphorylation",
"activation_binding/association",
"indirect effect",
"activation_compound",
"activation_ubiquination"
)
betas=c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
names(betas)=rel

for(p in 1:length(paths)){
L[[p]]<-NA
info<-rbind(info,c(NA,NA))
mapkpathway<-try(KEGGgraph::parseKGML(paste(mdir,paths[p],sep="/")),TRUE)

if (class(mapkpathway)!="try-error"){

info[p,]<-c(mapkpathway@pathwayInfo@number,mapkpathway@pathwayInfo@title)

if(length(mapkpathway@edges)>=1){


nodlst<-nodes(mapkpathway)
G<-list()
edg<-mapkpathway@edges
refty<-unlist(lapply(nodlst,function(x){x@type}))



for(kk in 1:length(nodlst)){
if((nodlst[[kk]])@type%in%c("gene","group")){
 if((nodlst[[kk]])@type=="gene"){
  tmpp<-(nodlst[[kk]])@name 
  G[[kk]]<-unlist(lapply(strsplit(tmpp,split=":"),function(x){x[2]}))
 } 
 if((nodlst[[kk]])@type=="group"){
  gc=na.omit((nodlst[[kk]])@component[refty[(nodlst[[kk]])@component]=="gene"])
  if(length(gc)>=1){
  tmpp<- unlist(lapply((nodlst[gc]),function(x){x@name}))
  G[[kk]]<-unlist(lapply(strsplit(tmpp,split=":"),function(x){x[2]}))
  }else{G[[kk]]<-"other"}
 } 

 }else{G[[kk]]<-"other"}
} 

nods<-names(nodes(mapkpathway))
names(G)<-nods



edg<-edges(mapkpathway)
N<-length(nods)



LT<-list()
for ( j in 1:length(rel)){
LT[[j]]<-matrix(0,N,N)
rownames(LT[[j]])<-colnames(LT[[j]])<-nods
}
names(LT)<-rel



for(i in 1:length(edg)){
 tmp<-edg[i]$rel
 tmp2<-edg[i]$rel@subtype
if(length(tmp2)>0){
 zz<-NULL

for(k in 1:length(tmp2)){
zz<-c(zz,tmp2[k]$subtype@name)
}
zz<-sort(zz)
zz<-zz[!zz%in%c(""," ")]
zz=paste(zz,collapse="_")

if(zz=="dephosphorylation_inhibition"){zz="inhibition_dephosphorylation"}
if(zz=="indirect effect_inhibition"){zz="inhibition_indirect effect"}
if(zz=="activation_indirect effect_phosphorylation"){zz=c("activation_indirect effect","activation_phosphorylation")}

for(jj in 1:length(rel)){
LT[[jj]][tmp@entry2ID,tmp@entry1ID]<-ifelse(rel[jj]%in%zz,1,0)
}
alr <-c(alr,paste(zz,collapse="_"))
}#end for edges


LT[["nodes"]]<-G
LT[["title"]]<-mapkpathway@pathwayInfo@title
LT[["NumberOfReactions"]]<-length(mapkpathway@reactions)

}
L[[p]]<-LT


} #end if no relations

} # end if error in xml parsing

}

names(L)<-paths

nok<-is.na(L)
L[nok]<-NULL
info<-info[!nok,]
info=rbind(info)
names(L)<-info[,1]


if (length(L)>0){

sumrel<-unlist(lapply(L,function(x){sum(unlist(lapply(x[c(1,3:23)],sum)))>0}))
okk=unlist(lapply(L,function(x){x$NumberOfReactions}))<1 & sumrel
info<-info[okk,]
L<-L[okk]
}else{return(FALSE)}


if (length(L)>0){

#prune pathways
nL<-L

for (ll in 1:length(L)){

cornodes<-L[[ll]]$nodes[!is.na(L[[ll]]$nodes) & L[[ll]]$nodes!="other"]
allgns<- unique(unlist(cornodes))

for(re in rel){
L[[ll]][[re]]<-L[[ll]][[re]][names(cornodes),names(cornodes)]
nL[[ll]][[re]]<-matrix(0,length(allgns),length(allgns))
rownames(nL[[ll]][[re]])<-colnames(nL[[ll]][[re]])<-allgns
for(n in 1:length(cornodes)){
if(cornodes[n]!="other"){
 nd<- names(cornodes)[n]
 s=names(L[[ll]][[re]][,nd][L[[ll]][[re]][,nd]!=0])
 length(s)

 if (length(s)>=1){
  for (jb in 1:length(s)){
   frm=unlist(cornodes[n]) 
   to=unlist(cornodes[s[jb]])
  if(!is.null(frm)& !is.null(to)){
   IND<-cbind(rep(frm,each=length(to)),rep(to,length(frm)))
   for(kb in 1:dim(IND)[1]){
   nL[[ll]][[re]][IND[kb,2],IND[kb,1]]<-1
   }
  }#end if
  }
 }#if downstream genes
}#if valid node
}#for nodes
}#for reltion
nL[[ll]]$nodes<-allgns

}#for path
#save(nL,file=paste(organism,"nL.RData",sep=""))

sumrel<-unlist(lapply(nL,function(x){sum(unlist(lapply(x[names(betas[betas!=0])],sum)))>0}))

path.info<-nL[unlist(lapply(nL,function(x){x$NumberOfReactions}))<1 & sumrel]


if (length(path.info)>=1 ){
if(is.null(out.path)){
save(path.info,file=paste("./",organism,"SPIA.RData",sep=""), compress="xz")
}else{
save(path.info,file=paste(out.path,"/",organism,"SPIA.RData",sep=""), compress="xz")
}
#save(alr,file=paste(mydir,"/",organism,"SPIA.RData",sep=""))

return(TRUE)
}else{cat (paste("there are no valid pathways in",organism) );
 return(FALSE)}

}else{return(FALSE)}


}else{cat (paste("there are no xml files in",organism) )
 return(FALSE)
} #end what if no xml files

}#end function


###########################################################################


# makeSPIAdata<-function(kgml.path="./hsa",organism="hsa",out.path="."){

# mdir=kgml.path

# alr<-NULL

# L<-list()

# info<-NULL
# path.nodes<-list()

# paths<-dir(mdir,pattern=".xml")

# if (length(paths)>0){


# #new rel
# rel<-c("activation",
# "compound",
# "binding/association",
# "expression",
# "inhibition",
# "activation_phosphorylation",
# "phosphorylation",
# "inhibition_phosphorylation",
# "inhibition_dephosphorylation",
# "dissociation",
# "dephosphorylation",
# "activation_dephosphorylation",
# "state change",
# "activation_indirect effect",
# "inhibition_ubiquination",
# "ubiquination",
# "expression_indirect effect",
# "inhibition_indirect effect",
# "repression",
# "dissociation_phosphorylation",
# "indirect effect_phosphorylation",
# "activation_binding/association",
# "indirect effect",
# "activation_compound",
# "activation_ubiquination"
# )
# betas=c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
# names(betas)=rel

# for(p in 1:length(paths)){
# L[[p]]<-NA
# info<-rbind(info,c(NA,NA))
# mapkpathway<-try(parseKGML(paste(mdir,paths[p],sep="/")),TRUE)

# if (class(mapkpathway)!="try-error"){

# info[p,]<-c(mapkpathway@pathwayInfo@number,mapkpathway@pathwayInfo@title)

# if(length(mapkpathway@edges)>=1){


# nodlst<-nodes(mapkpathway)
# G<-list()
# edg<-mapkpathway@edges
# refty<-unlist(lapply(nodlst,function(x){x@type}))



# for(kk in 1:length(nodlst)){
# if((nodlst[[kk]])@type%in%c("gene","group")){
#  if((nodlst[[kk]])@type=="gene"){
#   tmpp<-(nodlst[[kk]])@name 
#   G[[kk]]<-unlist(lapply(strsplit(tmpp,split=":"),function(x){x[2]}))
#  } 
#  if((nodlst[[kk]])@type=="group"){
#   gc=na.omit((nodlst[[kk]])@component[refty[(nodlst[[kk]])@component]=="gene"])
#   if(length(gc)>=1){
#   tmpp<- unlist(lapply((nodlst[gc]),function(x){x@name}))
#   G[[kk]]<-unlist(lapply(strsplit(tmpp,split=":"),function(x){x[2]}))
#   }else{G[[kk]]<-"other"}
#  } 

#  }else{G[[kk]]<-"other"}
# } 

# nods<-names(nodes(mapkpathway))
# names(G)<-nods



# edg<-edges(mapkpathway)
# N<-length(nods)



# LT<-list()
# for ( j in 1:length(rel)){
# LT[[j]]<-matrix(0,N,N)
# rownames(LT[[j]])<-colnames(LT[[j]])<-nods
# }
# names(LT)<-rel



# for(i in 1:length(edg)){
#  tmp<-edg[i]$rel
#  tmp2<-edg[i]$rel@subtype
# if(length(tmp2)>0){
#  zz<-NULL

# for(k in 1:length(tmp2)){
# zz<-c(zz,tmp2[k]$subtype@name)
# }
# zz<-sort(zz)
# zz<-zz[!zz%in%c(""," ")]
# zz=paste(zz,collapse="_")

# if(zz=="dephosphorylation_inhibition"){zz="inhibition_dephosphorylation"}
# if(zz=="indirect effect_inhibition"){zz="inhibition_indirect effect"}
# if(zz=="activation_indirect effect_phosphorylation"){zz=c("activation_indirect effect","activation_phosphorylation")}

# for(jj in 1:length(rel)){
# LT[[jj]][tmp@entry2ID,tmp@entry1ID]<-ifelse(rel[jj]%in%zz,1,0)
# }
# alr <-c(alr,paste(zz,collapse="_"))
# }#end for edges


# LT[["nodes"]]<-G
# LT[["title"]]<-mapkpathway@pathwayInfo@title
# LT[["NumberOfReactions"]]<-length(mapkpathway@reactions)

# }
# L[[p]]<-LT


# } #end if no relations

# } # end if error in xml parsing

# }

# names(L)<-paths

# nok<-is.na(L)
# L[nok]<-NULL
# info<-info[!nok,]
# info=rbind(info)
# names(L)<-info[,1]


# if (length(L)>0){

# sumrel<-unlist(lapply(L,function(x){sum(unlist(lapply(x[c(1,3:23)],sum)))>0}))
# okk=unlist(lapply(L,function(x){x$NumberOfReactions}))<1 & sumrel
# info<-info[okk,]
# L<-L[okk]
# }else{return(FALSE)}


# if (length(L)>0){

# #prune pathways
# nL<-L

# for (ll in 1:length(L)){

# cornodes<-L[[ll]]$nodes[!is.na(L[[ll]]$nodes) & L[[ll]]$nodes!="other"]
# allgns<- unique(unlist(cornodes))

# for(re in rel){
# L[[ll]][[re]]<-L[[ll]][[re]][names(cornodes),names(cornodes)]
# nL[[ll]][[re]]<-matrix(0,length(allgns),length(allgns))
# rownames(nL[[ll]][[re]])<-colnames(nL[[ll]][[re]])<-allgns
# for(n in 1:length(cornodes)){
# if(cornodes[n]!="other"){
#  nd<- names(cornodes)[n]
#  s=names(L[[ll]][[re]][,nd][L[[ll]][[re]][,nd]!=0])
#  length(s)

#  if (length(s)>=1){
#   for (jb in 1:length(s)){
#    frm=unlist(cornodes[n]) 
#    to=unlist(cornodes[s[jb]])
#   if(!is.null(frm)& !is.null(to)){
#    IND<-cbind(rep(frm,each=length(to)),rep(to,length(frm)))
#    for(kb in 1:dim(IND)[1]){
#    nL[[ll]][[re]][IND[kb,2],IND[kb,1]]<-1
#    }
#   }#end if
#   }
#  }#if downstream genes
# }#if valid node
# }#for nodes
# }#for reltion
# nL[[ll]]$nodes<-allgns

# }#for path
# #save(nL,file=paste(organism,"nL.RData",sep=""))

# sumrel<-unlist(lapply(nL,function(x){sum(unlist(lapply(x[names(betas[betas!=0])],sum)))>0}))

# path.info<-nL[unlist(lapply(nL,function(x){x$NumberOfReactions}))<1 & sumrel]


# if (length(path.info)>=1 ){
# if(is.null(out.path)){
# save(path.info,file=paste(system.file("extdata",package="SPIA"),"/",organism,"SPIA.RData",sep=""))
# }else{
# save(path.info,file=paste(out.path,"/",organism,"SPIA.RData",sep=""))
# }
# #save(alr,file=paste(mydir,"/",organism,"SPIA.RData",sep=""))

# return(TRUE)
# }else{cat (paste("there are no valid pathways in",organism) );
#  return(FALSE)}

# }else{return(FALSE)}


# }else{cat (paste("there are no xml files in",organism) )
#  return(FALSE)
# } #end what if no xml files

# }#end function


#########

















