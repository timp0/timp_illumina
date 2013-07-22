hclust.ward <- function(x) {
  hc = hclust(x, method="ward")
}


bump2grange <- function(dmr) {

  ##rate are the different p values, I think?
  rng=GRanges(seqnames=dmr$chr, ranges=IRanges(start=dmr$start, end=dmr$end))
  values(rng)=dmr[,4:14]
  return(rng)
}

timp.probeanno <- function(dat, codedir="~/Code/timp_illumina") {
  load(file.path(codedir, "timp_illumina_data", "probe_obj_final.rda"))
  probey=rowData(dat)
  values(probey)$dist.island=values(gprobes[match(names(probey), values(gprobes)$name)])$dist.island
  values(probey)$islrelate="OpenSea"
  values(probey)$islrelate[values(probey)$dist.island<4001]="Shelf"
  values(probey)$islrelate[values(probey)$dist.island<2001]="Shore"
  values(probey)$islrelate[values(probey)$dist.island==0]="Island"
  ##Send back
  rowData(dat)=probey
  return(dat)
}

  
dat.melt <- function(dat, logit=T) {
  require(reshape)
  ##This funciton melts data into a ggplot like format

  ##Log2 or beta
  if (logit) {
    melted=melt(getM(dat))
  } else {
    melted=melt(getBeta(dat))
  }

  ##Label columns apporpriately of melted matrix
  names(melted)=c("pid", "sid", "value")
  panno=rowData(dat)
  sanno=as.data.frame(colData(dat))
  ##Add probe annotation - because duplicates in rownames of gprobes->data.frame, row.names=NULL
  melted=cbind(melted, as.data.frame(panno[match(melted$pid, names(panno))], row.names=NULL))
  ##Add sample annotation
  melted=cbind(melted, as.data.frame(sanno[match(melted$sid, rownames(sanno)),]))
  return(melted)
}

range.plot <- function(dat, tab, grp="Status", grp2="Sample.ID", logit=T, num.plot=25) {
  ##Incoming tab is a GRanges

  require(GenomicRanges)
  require(RColorBrewer)
  require(ggplot2)

  ##Plot first 25 blocks
  M=min(length(tab), num.plot)

  ##Plot this far (in %) on either side
  ADD=0.1

  ##Reference for probe coloring
  probe.type=data.frame(lab=c("Island", "Shore", "Shelf", "OpenSea"), color=c("green", "blue", "orange", "red"),stringsAsFactors=F)                           


  ##Reference for sample coloring
  sample.type=data.frame(lab=as.character(unique(colData(dat)[[grp]])), stringsAsFactors=F)
  sample.type$color=brewer.pal(9, "Set1")[seq(along=sample.type$lab)]

  coly=c(probe.type$color, sample.type$color)
  names(coly)=c(probe.type$lab, sample.type$lab)
  
  for (i in 1:M) {
    ##Set range over which we will plot
    extra.width=max(width(tab[i])*(1+ADD), 5e3)
    plot.range=resize(tab[i], width=extra.width, fix="center")

    ##Find probes in that region
    pprobes=which(rowData(dat) %over% plot.range)

    subdat=dat[pprobes,]

    melted=dat.melt(subdat, logit=logit)

    rect=data.frame(xmin=start(tab[i]), xmax=end(tab[i]),
      ymin=-Inf, ymax=Inf)  
    
    to.plot=ggplot()+theme_bw()+theme(panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
      panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())
      labs(title=paste0("Region:", i, " Chromsome:",as.character(seqnames(plot.range))))

    ##region overlay
    to.plot=to.plot+geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), colour="grey", alpha=.05)

    
    ##If too many points, just plot lines(saves ugly block pictures)
    if (length(pprobes)>50) {
      to.plot=to.plot+geom_line(data=melted, stat="smooth", aes_string(x="start", y="value", colour=grp, group=grp2), alpha=.2, method="loess", span=.1)
      to.plot=to.plot+geom_line(data=melted, stat="smooth", aes_string(x="start", y="value", colour=grp), size=1, alpha=1, method="loess", span=.1)     
    } else {
      if (length(pprobes)>2) {
        to.plot=to.plot+stat_smooth(data=melted, aes_string(x="start", y="value", colour=grp, fill=grp), method="loess", alpha=.1, span=100)
        ##, se=F, alpha=.1, method="loess", span=.1)
      }
      to.plot=to.plot+geom_jitter(data=melted, aes_string(x="start", y="value", colour=grp, fill=grp), alpha=0.5)
    }
    
    ##Add rug of probes
    to.plot=to.plot+geom_rug(data=melted, aes(x=start, y=NULL, color=islrelate))

    to.plot=to.plot+scale_fill_manual(values=coly, guide=F)+scale_color_manual(values=coly)

    print(to.plot+scale_y_continuous(breaks=c(-0,.2, .4, .6, .8, 1))+coord_cartesian(ylim=c(0,1)))
    
  }
  
}

anno.region.plot <- function(dat, tab, grp="status", logit=T, num.plot=25, codedir="~/Code/timp_illumina") {
  ##Plot with annotaiton objects
  
  require(Gviz)
  
  ##Plot first 25 blocks
  M=min(length(tab), num.plot)

  ##Plot this far (in %) on either side
  ADD=0.1

  sampy=colData(dat)

  ##Gene Annotation
  load(file.path(codedir, "timp_illumina_data", "gene_island.rda"))
  values(refseq.exons)$exon=paste(values(refseq.exons)$refseq.name, values(refseq.exons)$exon.number, sep=".")

  for (i in 1:M) {

    ##Set range over which we will plot
    extra.width=max(width(tab[i])*(1+ADD), 5e3)
    plot.area=resize(tab[i], width=extra.width, fix="center")
    
    ##Find probes in that region

    subdat=dat[which(rowData(dat) %over% plot.area),]
    
    pprobes=rowData(subdat)
    
    if (logit) {
      yy=getM(subdat)
    } else {
      yy=getBeta(subdat)
    }
    
    chromy=as.character(seqnames(plot.area))
    
    ##TODO: Use size of dots to make small dots, maybe also alpha for dots
    ##Data track        
    if (length(pprobes)>10) {
      mtrack=DataTrack(pprobes, data=t(yy),
        genome="hg19", name="Beta",groups=sampy[[grp]], type="smooth", legend=T)
    } else {
      mtrack=DataTrack(pprobes, data=t(yy),
        genome="hg19", name="Beta",groups=sampy[[grp]], type=c("p", "a"))
    }

    ##Region track
    rtrack=AnnotationTrack(subsetByOverlaps(tab, plot.area),
      collapse=T, showId=F, stacking="dense", name="Regions", genome="hg19")

    ##Probe track
    ptrack=AnnotationTrack(pprobes, genome="hg19", name="Probes",
      feature=values(pprobes)$islrelate, collapse=T, mergeGroups=T, showId=F,
      stacking="dense",Shore="green", Island="blue", OpenSea="red", Shelf="orange")

    ##CpG island track
    isltrack=AnnotationTrack(subsetByOverlaps(ucsc.isl, plot.area),
      collapse=T, showId=F, stacking="dense", name="CpG Islands", genome="hg19")

    genes.overlapped=values(subsetByOverlaps(refseq.genes, plot.area))$refseq.name

    if (length(genes.overlapped)>0) {
      exons.overlapped=refseq.exons[(values(refseq.exons)$refseq.name %in% genes.overlapped)&(seqnames(refseq.exons)==chromy)]
      
      exon.frame=data.frame(start=start(exons.overlapped), end=end(exons.overlapped), symbol=values(exons.overlapped)$gene.name,
        exon=values(exons.overlapped)$exon, strand=as.character(strand(exons.overlapped)),
        feature="cds", transcript=values(exons.overlapped)$refseq.name, gene=values(exons.overlapped)$refseq.name,
        stringsAsFactors=F)
      
      genetrack=GeneRegionTrack(exon.frame, chromosome=chromy, genome="hg19",name="RefSeq Genes", showId=T, stacking="pack")
    } else {
      ##Empty track
      genetrack=AnnotationTrack()
    }
    
    itrack=IdeogramTrack(genome="hg19", chromosome=chromy)
    gtrack=GenomeAxisTrack()

    
    plotTracks(list(itrack, gtrack, genetrack, mtrack, rtrack, isltrack, ptrack) ,background.title="darkblue", from=start(plot.area), to=end(plot.area))
  }

}


st.region.plot <- function(dat, tab, grp="status", logit=T, num.plot=25, codedir="~/Code/timp_illumina") {
  ##Plot using standard package

  require(RColorBrewer)
  
  ##Plot first 25 regions
  M=min(length(tab), num.plot)

  ##Plot this far (in %) on either side
  ADD=0.1

  coly=brewer.pal(9, "Set1")

  ##Gene Annotation
  load(file.path(codedir, "timp_illumina_data", "gene_island.rda"))

  for (i in 1:M) {
    
    ##Set range over which we will plot
    extra.width=max(width(tab[i])*(1+ADD), 5e3)
    plot.area=resize(tab[i], width=extra.width, fix="center")
    
    ##Find probes in that region

    subdat=dat[which(rowData(dat) %over% plot.area),]
    
    pprobes=as.data.frame(rowData(subdat))
    pprobes$col="red"
    pprobes$col[pprobes$islrelate=="Shelf"]="orange"
    pprobes$col[pprobes$islrelate=="Shore"]="green"
    pprobes$col[pprobes$islrelate=="Island"]="blue"
    
    melted=dat.melt(subdat, logit=logit)
    
    if (logit) {
      ylim=range(melted$value)
    } else {
      ylim=c(0,1)
    }
    
    meltsum=ddply(melted, c(grp, "start"), function(x) {data.frame(val=median(x$value), bot=quantile(x$value, .05), top=quantile(x$value, .95))})


    meltsum$col=coly[factor(meltsum$anno)]
    
    
    plot(0,0, xlim=c(start(plot.area), end(plot.area)), ylim=ylim)
    
    d_ply(meltsum, grp, function(x) {lines(x$start, x$val, col=x$col[1], lwd=2);
                                   polygon(c(x$start, rev(x$start)), c(x$top, rev(x$bot)), border=NA, col=paste0(x$col[1], "33"))})

    rug(pprobes$start, col=pprobes$col, lwd=2)


    title(paste0("Chrom:", melted$seqnames[1], " Reg:", i))
    
    
  }
  
}


cg.dendro <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer"), p.thresh=1e-5, r.thresh=2.5) {
  ##This function makes a linkage tree uses unsupervised heirarchical clustering

  require(ggplot2)
  require(limma)
  require(RColorBrewer)
  require(gplots)
  
  ##Do limma t-test on probes
  probes=cg.dmtest(dat, ccomp=ccomp, grps=grps)

  goody=values(probes)$pv < p.thresh & abs(values(probes)$coef)>r.thresh

  sub=dat[goody,]

  y=getM(sub)

  d=dist(t(y), method="euclidean")
  fit=hclust(d, method="ward")

  coly=brewer.pal(9, "Set1")

  sampy=factor(colData(dat)$anno)

  heatmap.2(y, dendrogram="column", trace="none", labCol=colData(sub)$Sample.ID,
                        ColSideColors=coly[sampy],  labRow="", col=brewer.pal(11, "RdYlBu"))

  
  legend(x="bottomleft", legend=levels(sampy), col=coly[factor(levels(sampy))], pch=15)
  
  
  plot(fit, labels=colData(sub)$anno)

  
}

reg.cluster <- function(dat, tab, ccomp="Phenotype", namey="testplot") {
  ##This function does mds and scatter based on regions

  pprobes=which(rowData(dat) %over% tab)

  sub=dat[pprobes,]
  
  fullmd=cmdscale(dist(t(getM(dat[pprobes,]))))
  fullmd.probes=data.frame(x=fullmd[,1], y=fullmd[,2], outcome=colData(dat)[[ccomp]])
  
    
  print(ggplot(fullmd.probes, aes(x=x, y=y, colour=outcome))+geom_point(size=2.5) + 
        theme_bw()+labs(title=namey) + scale_colour_brewer(type="qual", palette="Dark2"))

}

cg.cluster <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer"), p.thresh=1e-5, r.thresh=2.5, volcano=F) {
  ##This function does mds of probes which show a difference, seperated by regional differences

  require(ggplot2)
  require(limma)
  require(plyr)

  
  ##Do limma t-test on probes
  probes=cg.dmtest(dat, ccomp=ccomp, grps=grps)

  ##Make volcano plot
  volcano.probes=as.data.frame(values(probes))
  if (volcano) {
    print(ggplot(volcano.probes, aes(x=coef, y=pv))+geom_point()+facet_wrap(~islrelate, scales="free")+
          theme_bw()+scale_y_log10()+labs(title="Volcano"))
  }
  goody=which(volcano.probes$pv < p.thresh & abs(volcano.probes$coef)>r.thresh)
  
  pass.probes=data.frame(idx=goody, type=volcano.probes$islrelate[goody])
  

  md.probes=ddply(pass.probes, .(type), function(x) {md=cmdscale(dist(t(getM(dat[x$idx,]))));
                                                  y=data.frame(outcome=colData(dat)[[ccomp]], 
                                                    x=md[,1], y=md[,2]); return(y)})
                                                    
  print(ggplot(md.probes, aes(x=x, y=y, colour=outcome))+geom_point() + facet_wrap(~type, scales="free")+
          theme_bw()+labs(title=paste(grps[1], grps[2], sep="-"))+ scale_colour_brewer(type="qual", palette="Set1"))
  
  fullmd=cmdscale(dist(t(getM(dat[pass.probes$idx,]))))
  fullmd.probes=data.frame(x=fullmd[,1], y=fullmd[,2], outcome=colData(dat)[[ccomp]])

  print(ggplot(fullmd.probes, aes(x=x, y=y, colour=outcome))+geom_point() +
          theme_bw()+labs(title=paste(grps[1], grps[2], sep="-")) + scale_colour_brewer(type="qual", palette="Set1"))
  
}


