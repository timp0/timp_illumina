function [] = plot_tissue_diff1(data,tn)
%Plot the tumor normal pairs
%For Charm - i=2


close all hidden;
%Extract the tissue types using sample class info

sampy=ismember(data.sample_class,[3-tn,5-tn,7-tn,9-tn,11-tn]);
sample_names=data.samples(sampy);
sample_tiss=data.avg(:,sampy);
%Split the first sample name to extract the tissue type we are playing
%with. - and to extract whether tumor or normal
for j=1:length(sample_names)
    broken=regexp(sample_names{j}, '_', 'split');
    tiss_colors(j).Labels=sample_names{j};
    switch broken{1}
        case 'Colon'
            tiss_colors(j).Colors='k';
        case 'Breast'
            tiss_colors(j).Colors='c';
        case 'Lung'
            tiss_colors(j).Colors='g';
        case 'Ovary'
            tiss_colors(j).Colors='r';
        case 'Wilms'
            tiss_colors(j).Colors='m';
    end
      
end
%Make the clustergram
cg=clustergram(sample_tiss, 'ColumnLabels', sample_names, 'RowLabels', ...
    data.genes,'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2);
set(cg, 'ColumnLabelsColor', tiss_colors);
if tn
    clustergram_plot2('tissue_tumor',1);
else
    clustergram_plot2('tissue_normal', 1);
end


end


    
    