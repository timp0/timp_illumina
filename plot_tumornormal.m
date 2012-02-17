function [] = plot_tumornormal(data)
%Plot the tumor normal pairs

for i=1:5
    %Extract the tumor normal pairs using sample_class index
    sampy=data.sample_class==(i*2)|data.sample_class==(i*2+1);
    sample_names=data.samples(sampy);
    sample_tn=data.avg(:,sampy);
    %Split the first sample name to extract the tissue type we are playing
    %with.
    broken=regexp(sample_names{1}, '_', 'split');
    %Make the clustergram
    cg=clustergram(sample_tn, 'ColumnLabels', sample_names, 'RowLabels', ...
        data.genes,'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2);
    
    clustergram_plot(['tn_' broken{1}]);
    
    close all hidden;
end

    
    