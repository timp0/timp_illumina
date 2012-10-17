rafablocks="/home/bst/faculty/ririzarr/projects/cegs/winston/storage/rdas"
mutdir="~/Dropbox/Data/Genetics/Mutations/092512_dl"

plotdir="~/Dropbox/Data/Genetics/Mutations/101112_mutationblock"

#Mutation tables
load(file.path(mutdir, "muts.rda"))

bmut=tcga.mutation$breast.full
muts=list(breast.mut=tcga.mutation$breast.full, colon.mut=tcga.mutation$colon.full, lung.mut=tcga.mutation$lung.full)

blocks=list()
#Breast try first
load(file.path(rafablocks, "blocks_breast_cancer_normal.rda"))
blocks$breast.block=b$tab[order(values(b$tab)$pvals)]

load(file.path(rafablocks, "blocks_colon_cancer_normal.rda"))
blocks$colon.block=b$tab[order(values(b$tab)$pvals)]

load(file.path(rafablocks, "blocks_lung_cancer_normal.rda"))
blocks$lung.block=b$tab[order(values(b$tab)$pvals)]

z=matrix(nrow=length(muts), ncol=length(blocks))
rownames(z)=names(muts)
colnames(z)=names(blocks)

for (i in 1:dim(z)[1]) {
  for (j in 1:dim(z)[2]) {
    this.block=blocks[[i]]
    this.block=this.block[values(this.block)$pvals<.05]
    this.mut=muts[[j]]
    z[i,j]=sum(countOverlaps(this.mut, this.block))/length(this.mut)
  }
}





