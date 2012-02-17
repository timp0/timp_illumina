function [] = illumina_genes_sub(data, good_index, good_info, samp_index,namey)

sample_names=data.samples(samp_index);

for j=1:length(sample_names)
    broken=regexp(sample_names{j}, '_', 'split');
    tn_colors(j).Labels=sample_names{j};
    if broken{2}(1)=='T'
            tn_colors(j).Colors='r';
    else
        tn_colors(j).Colors='b';
    end
end

if isstruct(good_index)
pvalue_table_csv(good_info, ['table_sub_' namey]);
end

sample_tn=data.avg(good_index,samp_index);
genie=data.genes(good_index);

%Make the clustergram
cg=clustergram(sample_tn, 'ColumnLabels', sample_names, 'RowLabels', ...
    genie,'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2);
set(cg, 'ColumnLabelsColor', tn_colors);
clustergram_plot(['tn_sub_' namey]);

close all hidden;
clear tn_colors;
