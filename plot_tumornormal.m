function [] = plot_tumornormal2(data,modif)
%Plot the tumor normal pairs

for i=1:5
    %Extract the tumor normal pairs using sample_class index
    sampy=data.sample_class==(i*2)|data.sample_class==(i*2+1);
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
    clustergram_plot(['tn_' modif '_' broken{1}]);
    
    close all hidden;
    clear tn_colors;
end

    
    