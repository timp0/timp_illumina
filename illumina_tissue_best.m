function [best_index, best_info] = illumina_tissue_best1(data)

i=2;

%Extract the tumor normal pairs using sample_class index
sampy=(data.sample_class==(i*2)|data.sample_class==(i*2+1))&data.sample_id<6;
normy=(data.sample_class==(i*2+1))&data.sample_id<6;
tumy=(data.sample_class==(i*2))&data.sample_id<6;

sample_names=data.samples(sampy);


tumor_data=data.avg(:,tumy);
normal_data=data.avg(:,normy);

for j=1:length(sample_names)
    broken=regexp(sample_names{j}, '_', 'split');
    tn_colors(j).Labels=sample_names{j};
    if broken{2}(1)=='T'
            tn_colors(j).Colors='r';
    else
        tn_colors(j).Colors='b';
        end
end

[best_index, best_info]=illumina_pvalue1(data.genes,['sub_' broken{1}],tumor_data,normal_data);

pvalue_table_csv(best_info, ['table_sub_' broken{1}]);

sample_tn=data.avg(best_index,sampy);
genie=data.genes(best_index);

%Make the clustergram
cg=clustergram(sample_tn, 'ColumnLabels', sample_names, 'RowLabels', ...
    genie,'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2);
set(cg, 'ColumnLabelsColor', tn_colors);
clustergram_plot(['tn_sub_' broken{1}]);

close all hidden;
clear tn_colors;

end
