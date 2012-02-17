function [best_index, best_info] = illumina_tissue_best1(data)

close all hidden;

i=2;

%Extract the tumor normal pairs using sample_class index
sampy=(data.sample_class==(i*2)|data.sample_class==(i*2+1))&data.sample_id<6;
normy=(data.sample_class==(i*2+1))&data.sample_id<6;
tumy=(data.sample_class==(i*2))&data.sample_id<6;

%Inelegant - but right now I don't care

sample_names=data.samples(sampy);

broken=regexp(sample_names{1}, '_', 'split');




tumor_data=data.avg(:,tumy);
normal_data=data.avg(:,normy);


[best_index, best_info]=illumina_pvalue1(data.genes,['sub_' broken{1}],tumor_data,normal_data);

illumina_genes_sub(data,best_index, best_info, sampy, broken{1});

end
