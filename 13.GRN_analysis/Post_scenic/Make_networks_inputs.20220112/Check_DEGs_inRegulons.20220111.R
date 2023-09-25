library(AUCell)
library(SCENIC)
library(SCopeLoomR)
library(GSEABase)

disease='Pancan'

#Set the path to the wdir:
wdir=''
loom=open_loom(paste(wdir,'/',disease,'_obj/',disease,'_pyscenic_output.v2.loom',sep=''),mode='r')

#Read information from loom file:
regulons_incidMat <- get_regulons(loom,column.attr.name = "Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)

sel_tfs=c('ZIC1','TBX5','NEUROD2','PAX6','NEUROD1','POU3F2','ZEB1','SOX8','SOX2')
sel_tfs=c('ZIC1','TBX5','NEUROD2','PAX6','NEUROD1')
#sel_tfs=c('SOX2','NFE2L3','TFCP2','NFIX','FOXP2')


all_links=NULL
all_nodes=NULL
for (i in 1:length(regulons)){
    tf=names(regulons[i])
    tf=gsub('(.*)\\(.*','\\1',tf)
    if (tf %in% sel_tfs){
    node_1=cbind(tf,tf,"tf")
    all_nodes=rbind(all_nodes,node_1)

    for (g in 1:length(regulons[i][[1]])){
    	gene=regulons[i][[1]][g]
	link=cbind(tf,gene,1,"affected")
	node_2=cbind(gene,gene,"gene")
	all_links=rbind(all_links,link)
	all_nodes=rbind(all_nodes,node_2)
}
print(i)
}
}

all_links=as.data.frame(all_links)
all_nodes=as.data.frame(all_nodes)
colnames(all_links)=c('from', 'to', 'weight', 'type')
colnames(all_nodes)=c('id','gene','type')

write.table(all_links,'out/links.tsv',sep='\t',quote=F,row.names=F)
write.table(all_nodes,'out/nodes.tsv',sep='\t',quote=F,row.names=F)