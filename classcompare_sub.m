function [] = classcompare_sub2(data,opty)
%Plot the tumor normal pairs for a subset of genes
%Also plot all samples
%opty.huge decides whether the plot should be large or not
%opty.cluster decides how to cluster(rows, columns, or both)

if nargin < 2
    opty.not=1;
end

if ~isfield(opty, 'cluster')
    opty.cluster=3;
end

if ~isfield(opty, 'modif')
    opty.modif='a';
end

if ~isfield(opty, 'huge')
    opty.huge=0;
end

if ~isfield(opty,'classes')
    opty.classes=[2 3 4 5 6 7 8 9 10 11];
end

if ~isfield(opty,'subset')
    opty.subset=true(length(data.genes.label),1);
end

if ~isfield(opty,'id_set')
    opty.id_set=unique(data.sample_id);
end



color_set='brgmkcybrgmkcybrgmkcy';


%Make the sample_class index
sampy=ismember(data.sample_class,opty.classes)&ismember(data.sample_id,opty.id_set);
    
sample_classes=data.sample_class(sampy);
sample_names=data.samples(sampy);
sample_dat=data.avg(opty.subset,sampy);


%Split the first sample name to extract the tissue type we are playing
%with. - and to extract whether tumor or normal

sam_colors=cell2struct(sample_names,'Labels', 2);


for j=1:length(sample_names)
    sam_colors(j).Colors=color_set(find(opty.classes==sample_classes(j),1));
end

%Make the clustergram
cg=clustergram(sample_dat, 'ColumnLabels', sample_names, 'RowLabels', ...
    data.genes.label(opty.subset),'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2,'cluster',opty.cluster);
set(cg, 'ColumnLabelsColor', sam_colors,'LabelsWithMarkers', 1);
clustergram_plot2(['class_' opty.modif], opty.huge);

close all hidden;
end
