function [] = classcompare_sub1(data,classes,subset,opty,modif)
%Plot the tumor normal pairs for a subset of genes
%Also plot all samples
%opty.huge decides whether the plot should be large or not
%opty.cluster decides how to cluster(rows, columns, or both)

color_set='brgmkcybrgmkcybrgmkcy';


%Make the sample_class index
sampy=ismember(data.sample_class,classes);
    
sample_classes=data.sample_class(sampy);
sample_names=data.samples(sampy);
sample_dat=data.avg(subset,sampy);


%Split the first sample name to extract the tissue type we are playing
%with. - and to extract whether tumor or normal

sam_colors=cell2struct(sample_names,'Labels', 2);


for j=1:length(sample_names)
    sam_colors(j).Colors=color_set(find(classes==sample_classes(j),1));
    
end

%Make the clustergram
cg=clustergram(sample_dat, 'ColumnLabels', sample_names, 'RowLabels', ...
    data.genes.label(subset),'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2,'cluster',opty.cluster);
set(cg, 'ColumnLabelsColor', sam_colors);
clustergram_plot2(['tn_' modif], opty.huge);

close all hidden;
end
