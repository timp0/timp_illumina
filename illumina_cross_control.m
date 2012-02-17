function [] = illumina_cross_control1(data)

just=(ismember(data.sample_class,[6 7])&(data.sample_id==2)) | ...
    (ismember(data.sample_class,[8 9])&(data.sample_id==19)) | ...
    (ismember(data.sample_class,[10 11])&(data.sample_id==13)) | ...
    (data.sample_class==12);



%Make the clustergram
clustergram(data.avg(:,just), 'ColumnLabels', data.samples(just), 'RowLabels', ...
    data.genes.label,'Linkage', 'ward', 'Dendrogram', [10 50], 'Standardize', 2,'cluster',3);
clustergram_plot2('cross_plot', 1);

close all hidden;