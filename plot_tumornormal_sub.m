function [] = plot_tumornormal_sub1(data,subset,opty,modif)
%Plot the tumor normal pairs

for i=1:5
    %Extract the tumor normal pairs using sample_class index
    sampy=data.sample_class==(i*2)|data.sample_class==(i*2+1);
    do_plot(sampy,data,subset,opty,modif);
end

sampy=(data.sample_class<12)&(data.sample_class>1);
do_plot(sampy, data, subset, opty, ['total_' modif]);

end

    
    
function [] = do_plot(sampy,data,subset,opty,modif)
sample_names=data.samples(sampy);
sample_tn=data.avg(subset,sampy);
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
    data.genes(subset),'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2,'cluster',opty.cluster);
set(cg, 'ColumnLabelsColor', tn_colors);
clustergram_plot2(['tn_' modif '_' broken{1}], opty.huge);

close all hidden;
clear tn_colors;
end
