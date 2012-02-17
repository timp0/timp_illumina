function [tissue_info] = illumina_tissue_best2a(data,modif)

close all hidden;

for i=1:5

    %Extract the tumor normal pairs using sample_class index
    sampy=(data.sample_class==(i*2)|data.sample_class==(i*2+1));
    normy=(data.sample_class==(i*2+1));
    tumy=(data.sample_class==(i*2));

    %Inelegant - but right now I don't care
    
    sample_names=data.samples(sampy);
    
    broken=regexp(sample_names{1}, '_', 'split');
    
    
    
    
    tumor_data=data.avg(:,tumy);
    normal_data=data.avg(:,normy);
    
    tic
    [best_index, best_info]=illumina_pvalue1(data.genes,['sub_' modif '_' broken{1}],tumor_data,normal_data);
    toc
    illumina_genes_sub(data,best_index, best_info, sampy, [broken{1} '_' modif]);
    
    tissue_info(i).best_index=best_index;
    tissue_info(i).best_info=best_info;
    close all hidden;
    clear best_index;
    clear best_info;
end


end
