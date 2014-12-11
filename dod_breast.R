library(minfi)
library(doMC)
library(reshape2)
library(ggplot2)
registerDoMC()
options(cores=6)


datapath="/mithril/Data/Infinium/CancProg_Feinberg"
plotdir="~/Dropbox/Data/Genetics/Infinium/120813_dodgrant"


if (!file.exists(file.path(plotdir, "breastdcis.rda"))) {

    plates=paste0("IL0", c("02", "05", "07", "08", "09", "10", "12", "45", "46"))
    
    baseDir=file.path(datapath, plates)
    
    ##Deleted all but v2
    ##Fix IL012 extra line at end
    
    targets=read.450k.sheet(baseDir)
    
    breast.targets=targets[targets$Tissue=="breast",]
    
    breast.raw=read.450k.exp(base=baseDir, targets=breast.targets)
    breast.pd=pData(breast.raw)
    
    
    qcReport(breast.raw, sampNames=breast.pd$Sample.ID, sampGroups=breast.pd$Phenotype, pdf=file.path(plotdir, "qcReport.pdf"))
    
    ##Ok -Breast_DCIS4 looks bad - exclude
    
    breast.raw=breast.raw[,breast.raw$Sample.ID!="Breast_DCIS_4"]
    breast.pd=breast.pd[breast.raw$Sample.ID!="Breast_DCIS_4",]
    
    ##Preprocess with SWAN(why not)?  Let's do Illumina for now

    ##Let's exclude hyperplastic, it's basically just normal with the samples we got from Argani
    breast.raw=breast.raw[,breast.raw$Phenotype!="hyperplastic"]
    breast.pd=breast.pd[breast.raw$Phenotype!="hyperplastic",]

    
    breast.dat=preprocessIllumina(breast.raw)
    
    breast.dat=mapToGenome(breast.dat)
    breast.ratio=ratioConvert(breast.dat)    
    breast.collapse=cpgCollapse(breast.ratio, what="Beta", returnBlockInfo=F)
    
    minfi.cluster=clusterMaker(chr=seqnames(breast.ratio), pos=start(breast.ratio))
    
    ##Find blocks 450k
    
    ##derived from amy email
    ##Just normal and DCIS for now
    breast.ratio=ratioConvert(breast.dat)[,breast.dat$Phenotype %in% c("normal", "adenoma")]
    breast.collapse=cpgCollapse(breast.ratio, what="Beta", returnBlockInfo=F)[,breast.dat$Phenotype %in% c("normal", "adenoma")]
    design=model.matrix(~breast.ratio$Phenotype)
    
    
    dcis.dmrs=bumphunter(breast.ratio, design=design, B=500, smooth=F, type="Beta", pickCutoff=T)        
    dcis.blocks=blockFinder(breast.collapse, design=design, B=500, what="Beta", pickCutoff=T)
    
    ##Now normal and cancer
    breast.ratio=ratioConvert(breast.dat)[,breast.dat$Phenotype %in% c("normal", "cancer")]
    breast.collapse=cpgCollapse(breast.ratio, what="Beta", returnBlockInfo=F)
    design=model.matrix(~breast.ratio$Phenotype)
    
    
    cancer.dmrs=bumphunter(breast.ratio, design=design, B=500, smooth=F, type="Beta", pickCutoff=T)        
    cancer.blocks=blockFinder(breast.collapse, design=design, B=500, what="Beta", pickCutoff=T)

    ##Now DCIS and cancer
    breast.ratio=ratioConvert(breast.dat)[,breast.dat$Phenotype %in% c("adenoma", "cancer")]
    breast.collapse=cpgCollapse(breast.ratio, what="Beta", returnBlockInfo=F)
    design=model.matrix(~breast.ratio$Phenotype)
    
    
    agg.dmrs=bumphunter(breast.ratio, design=design, B=500, smooth=F, type="Beta", pickCutoff=T)        
    agg.blocks=blockFinder(breast.collapse, design=design, B=500, what="Beta", pickCutoff=T)


    
    ##Load genes
    
    
    ##Ideas - Pie chart of location of DMRs
    ##Plot a DMR
    ##Plot a block
    
    save(list=c("agg.dmrs", "agg.blocks", "cancer.dmrs", "cancer.blocks", "dcis.dmrs", "dcis.blocks", "breast.dat"), compress="gzip", file=file.path(plotdir, "breastdcis.rda"))
} else {
    load(file.path(plotdir, "breastdcis.rda"))
}

pcutoff=.05
dcis.block.gr=GRanges(seqnames=dcis.blocks$table$chr, ranges=IRanges(start=dcis.blocks$table$start, end=dcis.blocks$table$end))
##P val cutoff
dcis.dmr.gr=GRanges(seqnames=dcis.dmrs$table$chr[dcis.dmrs$table$p.value<pcutoff], ranges=IRanges(start=dcis.dmrs$table$start[dcis.dmrs$table$p.value<pcutoff],
                                                                                       end=dcis.dmrs$table$end[dcis.dmrs$table$p.value<pcutoff]))

cancer.block.gr=GRanges(seqnames=cancer.blocks$table$chr, ranges=IRanges(start=cancer.blocks$table$start, end=cancer.blocks$table$end))
cancer.dmr.gr=GRanges(seqnames=cancer.dmrs$table$chr[dcis.dmrs$table$p.value<pcutoff],
    ranges=IRanges(start=cancer.dmrs$table$start[dcis.dmrs$table$p.value<pcutoff], end=cancer.dmrs$table$end[dcis.dmrs$table$p.value<pcutoff]))

breast.beta=getBeta(breast.dat)
probe.gr=granges(breast.dat)

if (FALSE) {
    pdf(file.path(plotdir, "dcis_block1.pdf"), width=11, height=8.5)
    
    for (i in 1:10) {
        rprobes=overlapsAny(probe.gr,resize(dcis.block.gr[i], width=width(dcis.block.gr[i])*3, fix="center"))
        reg=data.frame(x1=start(dcis.block.gr[i]), x2=end(dcis.block.gr[i]), y1=0, y2=1)
        subdata=breast.beta[rprobes,]
        colnames(subdata)=breast.dat$Phenotype
        rownames(subdata)=start(probe.gr[rprobes])
        mdata=melt(subdata)
        
        print(ggplot(mdata)+geom_smooth(aes(x=Var1, y=value, fill=Var2, color=Var2),method="loess")+
              geom_point(aes(x=Var1, y=value, fill=Var2, color=Var2),size=1, alpha=.5)+theme_bw()+scale_y_continuous(limits=c(0,1))+
              geom_rect(data=reg, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="orange", fill="orange", alpha=.2))
        
    }
    dev.off()

}

if (FALSE) {

    ##Histogram of block probes vs. not block probes
    bprobes=overlapsAny(probe.gr,dcis.block.gr)

    subdata=breast.beta[bprobes,]
    colnames(subdata)=breast.dat$Sample.ID
    rownames(subdata)=start(probe.gr[bprobes])
    mdata=melt(subdata)
    mdata$type="Normal"
    mdata$type[grepl("Cancer", mdata$Var2)]="Cancer"
    mdata$type[grepl("DCIS", mdata$Var2)]="DCIS"

    subavg=data.frame(Normal=rowMeans(subdata[,breast.dat$Phenotype=="normal"]),
        DCIS=rowMeans(subdata[,breast.dat$Phenotype=="adenoma"]),
        Cancer=rowMeans(subdata[,breast.dat$Phenotype=="cancer"]))

    msubavg=melt(subavg)

    pdf(file.path(plotdir, "dcis_block_hist1.pdf"))
    print(ggplot(mdata)+geom_density(aes(x=value, group=Var2, color=type),alpha=.5, linetype="dotted")+
          geom_density(data=msubavg, aes(x=value, color=variable, group=variable), size=1.5)+
          theme_bw())

    dev.off()
    
}


if (TRUE) {
    pdf(file.path(plotdir, "cancer_block1.pdf"), width=11, height=8.5)
    
    for (i in 1:50) {
        rprobes=overlapsAny(probe.gr,resize(cancer.block.gr[i], width=width(cancer.block.gr[i])*3, fix="center"))
        reg=data.frame(x1=start(cancer.block.gr[i]), x2=end(cancer.block.gr[i]), y1=0, y2=1)
        subdata=breast.beta[rprobes,]
        colnames(subdata)=breast.dat$Phenotype
        rownames(subdata)=start(probe.gr[rprobes])
        mdata=melt(subdata)
        
        print(ggplot(mdata)+geom_smooth(aes(x=Var1, y=value, fill=Var2, color=Var2),method="loess")+
              geom_point(aes(x=Var1, y=value, fill=Var2, color=Var2),size=1, alpha=.5)+theme_bw()+scale_y_continuous(limits=c(0,1))+
              geom_rect(data=reg, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="orange", fill="orange", alpha=.2))
        
    }
    dev.off()

}

if (TRUE) {
    pdf(file.path(plotdir, "cancer_dmr1.pdf"), width=11, height=8.5)
    
    for (i in 1:50) {
        rprobes=overlapsAny(probe.gr,resize(cancer.dmr.gr[i], width=width(cancer.dmr.gr[i])*3, fix="center"))
        reg=data.frame(x1=start(cancer.dmr.gr[i]), x2=end(cancer.dmr.gr[i]), y1=0, y2=1)
        subdata=breast.beta[rprobes,]
        colnames(subdata)=breast.dat$Phenotype
        rownames(subdata)=start(probe.gr[rprobes])
        mdata=melt(subdata)
        
        print(ggplot(mdata)+geom_smooth(aes(x=Var1, y=value, fill=Var2, color=Var2),method="loess")+
              geom_point(aes(x=Var1, y=value, fill=Var2, color=Var2),size=1, alpha=.5)+theme_bw()+scale_y_continuous(limits=c(0,1))+
              geom_rect(data=reg, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="orange", fill="orange", alpha=.2))
        
    }
    dev.off()

}


if (TRUE) {

    ##Histogram of block probes vs. not block probes
    bprobes=overlapsAny(probe.gr,cancer.block.gr)

    subdata=breast.beta[bprobes,]
    colnames(subdata)=breast.dat$Sample.ID
    rownames(subdata)=start(probe.gr[bprobes])
    mdata=melt(subdata)
    mdata$type="Normal"
    mdata$type[grepl("Cancer", mdata$Var2)]="Cancer"
    mdata$type[grepl("DCIS", mdata$Var2)]="DCIS"

    subavg=data.frame(Normal=rowMeans(subdata[,breast.dat$Phenotype=="normal"]),
        DCIS=rowMeans(subdata[,breast.dat$Phenotype=="adenoma"]),
        Cancer=rowMeans(subdata[,breast.dat$Phenotype=="cancer"]))

    msubavg=melt(subavg)

    pdf(file.path(plotdir, "cancer_block_hist1.pdf"))
    print(ggplot(mdata)+geom_density(aes(x=value, group=Var2, color=type),alpha=.5, linetype="dotted")+
          geom_density(data=msubavg, aes(x=value, color=variable, group=variable), size=1.5)+
          theme_bw())

    dev.off()
   
}

if (TRUE) {
    ##All probes
    
    subdata=breast.beta
    colnames(subdata)=breast.dat$Sample.ID
    rownames(subdata)=start(probe.gr[bprobes])
    mdata=melt(subdata)
    mdata$type="Normal"
    mdata$type[grepl("Cancer", mdata$Var2)]="Cancer"
    mdata$type[grepl("DCIS", mdata$Var2)]="DCIS"
    
    subavg=data.frame(Normal=rowMeans(subdata[,breast.dat$Phenotype=="normal"]),
        DCIS=rowMeans(subdata[,breast.dat$Phenotype=="adenoma"]),
        Cancer=rowMeans(subdata[,breast.dat$Phenotype=="cancer"]))
    
    msubavg=melt(subavg)
    
    pdf(file.path(plotdir, "cancer_hist1.pdf"))
    print(ggplot(mdata)+geom_density(aes(x=value, group=Var2, color=type),alpha=.5, linetype="dotted")+
          geom_density(data=msubavg, aes(x=value, color=variable, group=variable), size=1.5)+
          theme_bw())
    
    dev.off()
    
}

if (TRUE) {
    load("~/Code/timp_illumina/timp_illumina_data/gene_island.rda")

    prom=promoters(refseq.genes)

    dmr.regions=data.frame(prom=c(sum(overlapsAny(dcis.dmr.gr, prom)), sum(overlapsAny(cancer.dmr.gr, prom))),
        genes=c(sum(overlapsAny(dcis.dmr.gr, refseq.genes) & !overlapsAny(dcis.dmr.gr, prom)),
            sum(overlapsAny(cancer.dmr.gr, refseq.genes) & !overlapsAny(cancer.dmr.gr, prom))),
        inter=c(sum(!overlapsAny(dcis.dmr.gr, c(prom, refseq.genes))), sum(!overlapsAny(cancer.dmr.gr, c(prom,refseq.genes)))),
        tots=c(length(dcis.dmr.gr), length(cancer.dmr.gr))
        )
    rownames(dmr.regions)=c("DCIS", "Cancer")

    dmr.regions.freq=dmr.regions/dmr.regions$tots
    mdmrlocsf=melt(as.matrix(dmr.regions.freq[,1:3]))
    mdmrlocs=melt(as.matrix(dmr.regions[,1:3]))

    breast.ratio=ratioConvert(breast.dat)
    breast.collapse=cpgCollapse(breast.ratio, what="Beta", returnBlockInfo=F)   
    probey=granges(breast.collapse)

    dmr.cgi=rbind(table(probey$type[overlapsAny(probey, dcis.dmr.gr)]), table(probey$type[overlapsAny(probey, cancer.dmr.gr)]), table(probey$type))
    rownames(dmr.cgi)=c("DCIS", "Cancer", "Total")

    mdmrcgi=melt(as.matrix(dmr.cgi[1:2,]))
    
    pdf(file.path(plotdir, "dmr_locs1.pdf"))
    print(ggplot(mdmrlocs)+geom_bar(aes(x=Var2, y=value, group=Var1, fill=Var1), stat="identity", position=position_dodge(width = 0.7), width=.6)+
          theme_bw())
    print(ggplot(mdmrlocsf)+geom_bar(aes(x=Var2, y=value, group=Var1, fill=Var1), stat="identity", position=position_dodge(width = 0.7), width=.6)+
          theme_bw())

    print(ggplot(mdmrcgi)+geom_bar(aes(x=Var2, y=value, group=Var1, fill=Var1), stat="identity", position=position_dodge(width = 0.7), width=.6)+
          theme_bw())

    dev.off()



    
}


if (TRUE) {
    ##Venn diagram
    library(VennDiagram)
    ##Number of dmrs passing cuttoff

    pdf(file.path(plotdir, "dmr_venn1.pdf"))
    venn.plot=draw.pairwise.venn(length(dcis.dmr.gr), length(cancer.dmr.gr), sum(overlapsAny(dcis.dmr.gr, cancer.dmr.gr)))
    grid.draw(venn.plot);
    grid.newpage();
    dev.off()

}


##Select samples where 
