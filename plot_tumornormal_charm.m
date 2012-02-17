function [] = plot_tumornormal_charm(data,modif)
%Plot the tumor normal pairs
%For Charm - i=2

i=2;

%Extract the tumor normal pairs using sample_class index
sampy=(data.sample_class==(i*2)|data.sample_class==(i*2+1))&(data.sample_id<6);
sample_names=data.samples(sampy);
sample_tn=data.avg(:,sampy);
%Split the first sample name to extract the tissue type we are playing
%with. - and to extract whether tumor or normal
for j=1:length(sample_names)
    broken=regexp(sample_names{j}, '_', 'split');
    tn_colors(j).Labels=sample_names{j};
    if broken{2}(1)=='T'
            tn_colors(j).Colors='r';
    else
        tn_colors(j).Colors='b';
        end
end
%Make the clustergram
cg=clustergram(sample_tn, 'ColumnLabels', sample_names, 'RowLabels', ...
    data.genes,'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2);
set(cg, 'ColumnLabelsColor', tn_colors);
clustergram_plot(['tn_charm_' modif '_' broken{1}]);

close all hidden;
clear tn_colors;


    
    