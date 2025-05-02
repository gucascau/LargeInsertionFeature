###analyze the gene ontology
library(data.table)
library(ALL)
library(GO.db)
library(annotate)
library(genefilter)
library(GOstats)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(Rgraphviz)



                        ### load data
                        dt.final<-fread("/Users/yy/Desktop/ChenLab/Repair-seq/Data/Xin-data/Final_normalizedCount/Final_normlizated_list")
                        dt.final[,Standard_Name:=data.table(mapIds(org.Sc.sgd.db,
                                                                   keys=dt.final$V1, 
                                                                   column="GENENAME", 
                                                                   keytype="ENSEMBL",
                                                                   multiVals="first"))]
                        
                        dt.final.expr.melt<-melt(dt.final,id.vars = c("V1","Standard_Name"))
                        
                        samp.id<-fread("/Users/yy/Desktop/ChenLab/Repair-seq/Data/Xin-data/SampInfo/Sample_IDs.csv")
                        samp.id[Genotype=="fzoq",Genotype:="fzo1"]
                        setnames(samp.id,"Sample ID","variable")
                        
                        setkey(dt.final.expr.melt,variable)
                        setkey(samp.id,variable)
                        
                        ### logFC for each gene
                        dt.xin.expr.melt<-merge(dt.final.expr.melt,samp.id)
                        dt.xinexprs<-dt.xin.expr.melt[, list(ave.expr=mean(value,na.rm=T)),list(V1,Standard_Name,Genotype)]
                        
                        o.xin.expr.commp<-foreach(name=unique(dt.xinexprs$V1),.combine = "rbind")%do%{
                                  dt.tmp<-dt.xinexprs[V1==name]
                                  o.dt = data.table(V1=name,
                                                  Standard_Name=unique(dt.tmp$Standard_Name),
                                                  mgm1.wt=log2(dt.tmp[Genotype=="mgm1"]$ave.expr/dt.tmp[Genotype=="WT"]$ave.expr),
                                                  fzo1.wt=log2(dt.tmp[Genotype=="fzo1"]$ave.expr/dt.tmp[Genotype=="WT"]$ave.expr),
                                                  ugo1.wt=log2(dt.tmp[Genotype=="ugo1"]$ave.expr/dt.tmp[Genotype=="WT"]$ave.expr),
                                                  pif1_m2.wt=log2(dt.tmp[Genotype=="pif1-m2"]$ave.expr/dt.tmp[Genotype=="WT"]$ave.expr),
                                                  pif1_m2_dna2.wt=log2(dt.tmp[Genotype=="pif1-m2 dna2"]$ave.expr/dt.tmp[Genotype=="WT"]$ave.expr)
                                                  )
                        }
                                  
                        
                        
                        ### can skip
                        #mgm1<-fread("/Users/yy/Desktop/ChenLab/Repair-seq/Data/Xin-data/DEGs/Mgm1_WT_final_DEG.txt")
                        #Foz1<-fread("/Users/yy/Desktop/ChenLab/Repair-seq/Data/Xin-data/DEGs/Foz1_WT_final_DEG.txt")
                        #Ugo1<-fread("/Users/yy/Desktop/ChenLab/Repair-seq/Data/Xin-data/DEGs/Ugo1_WT_final_DEG.txt")
                        
                        
                        ### unversal ids  
                        keytypes(org.Sc.sgd.db)
                        
                        #this is the p value cut off
                        hgCutoff <-  0.05
                        
                        ### GO analysis on higher expressed genes in each strain
                        
                        o.go<-foreach(i=c("mgm1","Foz1","Ugo1"),.combine = "rbind")%do%{
                                dt.input<-fread(paste0("/Users/yy/Desktop/ChenLab/Repair-seq/Data/Xin-data/DEGs/",i,"_WT_final_DEG.txt"))
                                ### whether logFC > 0 or use all
                                params <- new("GOHyperGParams",
                                              geneIds=dt.input[logFC>0]$V1,
                                              universeGeneIds=dt.final$V1,
                                              annotation="org.Sc.sgd.db",
                                              ontology="BP",
                                              pvalueCutoff=hgCutoff,
                                              conditional= F,
                                              testDirection="over")
                                
                                hgOver <- hyperGTest(params)
                                df <- as.data.table(summary(hgOver,pvalue=0.05,categorySize=9))
                                
                                df[,strain:=i]
                        }
                        
                        o.go.up.iron<-o.go[,p:=(0-log10(Pvalue))][,logOR:=log2(OddsRatio)][grepl("iron|Iron",Term)]
                        
                        ### visualize
                        ggplot(o.go.up.iron, aes(logOR, factor(Term), color=factor(strain)))+ 
                          geom_point(size =2) + ylab("GO-terms")+xlab("log2 Odds Ratio")+theme_bw()+labs(color="Strains")
                        
                        
                        ### Sort all Go terms of significance and then look at gene expression directions between mutants and wild type
                        dt.mut.wt.xue<-dt.expr2[,BY4741_rep1:=log2(BY4741_rep1+0.1)][,mmm1_rep1:=log2(mmm1_rep1+0.1)][,mdm34_rep1:=log2(mdm34_rep1+0.1)]
                        dt.mut.wt.xue<-dt.mut.wt.xue[,mm1.wt:=log2(mmm1_rep1/BY4741_rep1)][,mdm34.wt:=log2(mdm34_rep1/BY4741_rep1)][,c("Systematic_Name","Standard_Name","mm1.wt","mdm34.wt")]
                        
                        ### Xue data using dt.expr2
                        dt.xue.mut.melt<-melt(dt.mut.wt.xue,id.vars = c("Systematic_Name","Standard_Name"))
                        
                        o.fcxue<-foreach(g=unique(o.go.up.iron$GOBPID),.combine = "rbind")%do%{
                                  go_id <- g
                                  allegs <- unique(get(go_id, org.Sc.sgdGO2ALLORFS))
                                  esids <- data.table(unlist(mget(allegs, org.Sc.sgdGENENAME)))
                                  
                                  o.strains<-dt.xue.mut.melt[Standard_Name%in%esids$V1]
                                  o.strains[,GO:=go_id]
                          
                        }
                        
                        dt.gons<-o.go[,c("GOBPID","Size","Term","p"),with=F][!duplicated(GOBPID)]
                        setnames(dt.gons,"GOBPID","GO")
                        
                        setkey(dt.gons,GO)
                        setkey(o.fcxue,GO)
                        
                        o.fcxue<-merge(o.fcxue,dt.gons)
                        
                        ### visualize
                        ggplot(o.fcxue, aes(value, factor(Term), color=factor(Term)))+ 
                          geom_boxplot(horizontal=TRUE) + facet_grid(variable~.)+geom_vline(xintercept = 0,linetype="dashed")+
                          ylab("GO-terms")+xlab("log2 Fold Change")+
                          theme_bw()+labs(color="Strains")+theme(legend.position = "non")
                        
                        
                        
                        ### Xin data using dt.final
                        o.fc<-foreach(g=unique(o.go.up.iron$GOBPID),.combine = "rbind")%do%{
                              go_id <- g
                              allegs <- unique(get(go_id, org.Sc.sgdGO2ALLORFS))
                              esids <- data.table(unlist(mget(allegs, org.Sc.sgdGENENAME)))
                              
                              o.strains<-o.xin.expr.commp[Standard_Name%in%esids$V1]
                              o.strains[,GO:=go_id]
                              
                        }
                        
                        o.fc$V1<-NULL
                        o.fcmelt<-melt(o.fc[!is.na(Standard_Name)],id.vars = c("GO","Standard_Name"))
                        
                        dt.gons<-o.go[,c("GOBPID","Size","Term","p"),with=F][!duplicated(GOBPID)]
                        setnames(dt.gons,"GOBPID","GO")
                        
                        setkey(dt.gons,GO)
                        setkey(o.fcmelt,GO)
                        
                        o.fcmelt<-merge(o.fcmelt,dt.gons)
                        
                        ### visualize
                        ggplot(o.fcmelt, aes(value, factor(Term), color=factor(Term)))+ 
                          geom_boxplot(horizontal=TRUE) + facet_grid(variable~.)+geom_vline(xintercept = 0,linetype="dashed")+
                          ylab("GO-terms")+xlab("log2 Fold Change")+
                          theme_bw()+labs(color="Strains")+theme(legend.position = "non")
                        
                        ### Xue data same GO terms grouped genes
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        ### make sample.id table for Xin data self remapping
                        xin.sample<-fread("/Users/yy/Desktop/ChenLab/Repair-seq/Data/Xin-data/SampInfo/Sample_IDs.csv")
                        xin.sample<-xin.sample[,c(1)]
                        
                        write.table(xin.sample,file="/Users/yy/Desktop/ChenLab/Repair-seq/Data/Xin-data/SampInfo/XinSampleIDMapping120422.txt",
                                    quote=F,row.names=F,col.names=F,sep="\t")
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        