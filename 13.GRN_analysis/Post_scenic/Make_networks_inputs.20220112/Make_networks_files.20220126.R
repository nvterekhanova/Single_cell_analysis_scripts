export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(AUCell)
library(SCENIC)
library(SCopeLoomR)
library(GSEABase)

disease='Pancan'
wdir='/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/9.GRN_analysis/Run.v.20220123'
loom=open_loom(paste(wdir,'/',disease,'_obj/',disease,'_pyscenic_output.loom',sep=''),mode='r')

###Regulons_ranks:
dis_s='BRCA_Basal'


regs_ranks=read.table(paste("/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/9.GRN_analysis/",
"Run.v.20220123/Post_scenic/Pancan_obj/out/Regulons_BinAct_difference.20220124.tsv",sep=''),header=T,sep='\t')

###Apply a bit different filtering (latest):
regs_ranks_s=regs_ranks[regs_ranks$cell_t1==dis_s & regs_ranks$diff>0.1 & regs_ranks$mean_score1>0.1 
& regs_ranks$FDR<0.05,]
regs_ranks_s=regs_ranks_s[order(-regs_ranks_s$diff),]

###Also, try filtering Regulons with few targets:
annot_r=read.table("../out/Regulon_annot.20220124.tsv",sep='\t',header=T)
rownames(annot_r)=annot_r$Regulon
annot_r=annot_r[regs_ranks_s$Regulon,]
annot_r$diff=regs_ranks_s$diff
annot_r$mean_score1=regs_ranks_s$mean_score1
annot_r$mean_score2=regs_ranks_s$mean_score2
annot_r$Regulon_1=gsub('\\(\\+\\)','',annot_r$Regulon)
annot_r$ID=paste(annot_r$Regulon_1,' (',annot_r$Genes_N,'g',')',sep='')
annot_r=annot_r[annot_r$Genes_N>=20,]
annot_r$Fch=annot_r$mean_score1/annot_r$mean_score2
annot_r=annot_r[order(-annot_r$Fch),]

sel_tfs=annot_r$Regulon[1:15]
sel_tfs=gsub('\\(\\+\\)','',sel_tfs)
sel_tfs=sel_tfs[!is.na(sel_tfs)]

# Read information from loom file:
regulons_incidMat <- get_regulons(loom,column.attr.name = "Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)


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

###using binact-diff approach:
dir.create(paste(dis_s,"/out",sep=''))
write.table(all_links, paste(dis_s,'/out/links_',dis_s,'.v1.tsv',sep=''),sep='\t',quote=F,row.names=F)
write.table(all_nodes, paste(dis_s,'/out/nodes_',dis_s,'.v1.tsv',sep=''),sep='\t',quote=F,row.names=F)


