library(minfi)
library(doMC)
library(reshape2)
library(ggplot2)
registerDoMC()
options(cores=6)


datapath="/mithril/Data/Infinium/020615_iceland3d"
plotdir="~/Dropbox/Data/Genetics/Infinium/021015_iceland"


if (!file.exists(file.path(plotdir, "iceland.rda"))) {

    plates="IL114"
    
    baseDir=file.path(datapath, plates)
    
    ##Deleted all but v2
    ##Fix IL012 extra line at end
    
    targets=read.450k.sheet(baseDir)
    
    iceland.raw=read.450k.exp(base=baseDir, targets=targets)
    iceland.pd=pData(iceland.raw)
    
    ##QC
    #qcReport(iceland.raw, sampNames=iceland.pd$Sample.ID, sampGroups=iceland.pd$Phenotype, pdf=file.path(plotdir, "qcReport.pdf"))
    ##Says sample too sparse?  ##Report to Kasper/Martin?
    
    ##Preprocess with SWAN(why not)?  Let's do Illumina for now

    iceland.dat=preprocessIllumina(iceland.raw)
    
    iceland.dat=mapToGenome(iceland.dat)
    iceland.ratio=ratioConvert(iceland.dat)    
    iceland.collapse=cpgCollapse(iceland.ratio, what="Beta", returnBlockInfo=F)
    
    minfi.cluster=clusterMaker(chr=seqnames(iceland.ratio), pos=start(iceland.ratio))
    
    ##Find blocks 450k
    
    ##What samples?
    ##2D v 3D
    iceland.dat$Tissue=c(rep("3D", 6), "2D", "2D")
    ##And Epithelial vs. Mesenchymal
    iceland.dat$Phenotype=c("e", "e", "e", "m","m","m", "e","m")
    iceland.ratio=ratioConvert(iceland.dat)[,iceland.dat$Phenotype %in% "e"]
    iceland.collapse=cpgCollapse(iceland.ratio, what="Beta", returnBlockInfo=F)
    design=model.matrix(~iceland.ratio$Tissue)
    
    
    e2d3d.dmrs=bumphunter(iceland.ratio, design=design, B=500, smooth=F, type="Beta", pickCutoff=T)        
    e2d3d.blocks=blockFinder(iceland.collapse, design=design, B=500, what="Beta", pickCutoff=T)

    iceland.ratio=ratioConvert(iceland.dat)[,iceland.dat$Phenotype %in% "m"]
    iceland.collapse=cpgCollapse(iceland.ratio, what="Beta", returnBlockInfo=F)
    design=model.matrix(~iceland.ratio$Tissue)

    m2d3d.dmrs=bumphunter(iceland.ratio, design=design, B=500, smooth=F, type="Beta", pickCutoff=T)        
    m2d3d.blocks=blockFinder(iceland.collapse, design=design, B=500, what="Beta", pickCutoff=T)

    iceland.ratio=ratioConvert(iceland.dat)[,iceland.dat$Tissue %in% "3D"]
    iceland.collapse=cpgCollapse(iceland.ratio, what="Beta", returnBlockInfo=F)
    design=model.matrix(~iceland.ratio$Phenotype)
    
    evm3d.dmrs=bumphunter(iceland.ratio, design=design, B=500, smooth=F, type="Beta", pickCutoff=T)        
    evm3d.blocks=blockFinder(iceland.collapse, design=design, B=500, what="Beta", pickCutoff=T)

    iceland.ratio=ratioConvert(iceland.dat)[,iceland.dat$Tissue %in% "2D"]
    iceland.collapse=cpgCollapse(iceland.ratio, what="Beta", returnBlockInfo=F)
    design=model.matrix(~iceland.ratio$Phenotype)
    
    evm2d.dmrs=bumphunter(iceland.ratio, design=design, B=500, smooth=F, type="Beta", pickCutoff=T)        
    evm2d.blocks=blockFinder(iceland.collapse, design=design, B=500, what="Beta", pickCutoff=T)

    ##Load genes
     
    ##Ideas - Pie chart of location of DMRs
    ##Plot a DMR
    ##Plot a block
    
    save(list=c("e2d3d.blocks", "e2d3d.dmrs", "m2d3d.blocks", "m2d3d.dmrs", "evm3d.blocks", "evm3d.dmrs", "evm2d.blocks", "evm2d.dmrs", "iceland.dat"), compress="gzip", file=file.path(plotdir, "iceland.rda"))
} else {
    load(file.path(plotdir, "iceland.rda"))
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
        #rprobes=overlapsAny(probe.gr,resize(cancer.dmr.gr[i], width=width(cancer.dmr.gr[i])*3, fix="center"))
        rprobes=overlapsAny(probe.gr,resize(cancer.dmr.gr[i], width=1e4, fix="center"))
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
