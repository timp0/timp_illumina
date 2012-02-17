function [diff_geneidx,goody] = illumina_pvalue1(genes,namey,tumor_data,normal_data)


nTumor=size(tumor_data,2);
nNormal=size(normal_data,2);


meanTumor=mean(tumor_data,2);
meanNormal=mean(normal_data,2);



pValues = mattest(tumor_data, normal_data, 'Permute', 1e5);


%The differential score of 13 corresponds to a p-value of 0.05, the differential score of 20
%corresponds to a p-value of 0.01, and the differential score of 30
%corresponds to a p-value of 0.001. A positive differential score
%represents up regulation, while a negative score represents down
%regulation. The differential score of 13 corresponds to a p-value of 0.05,
%the differential score of 20 corresponds to a p-value of 0.01, and the
%differential score of 30 corresponds to a p-value of 0.001. A positive
%differential score represents up regulation, while a negative score represents down regulation.

    
diffscore = @(p, avgTest, avgRef) -10*sign(avgTest - avgRef).*log10(p);

diffScores = diffscore(pValues, meanTumor, meanNormal);

[pFDR, qValues] = mafdr(pValues);

diffScoresFDRQ = diffscore(qValues, meanTumor, meanNormal);

diffStruct = mavolcanoplot(tumor_data, normal_data, qValues, 'pcutoff', 0.01, ...
    'LogTrans', false, 'Labels', genes, 'foldchange', 1.01,'PlotOnly', true);

volcano_plot(namey, 0);

close all hidden;

nDiffGenes=size(diffStruct.GeneLabels,1);
nGenes=size(genes,1);

diff_geneidx = false(nGenes, 1);

for i = 1:nDiffGenes
    index = find(strcmpi(diffStruct.GeneLabels{i},genes));
    
    diff_geneidx(index)=1;
    
    goody(i).gene=diffStruct.GeneLabels{i};
    goody(i).FDR=pFDR(i);
    goody(i).pValues=diffStruct.PValues(i);
    goody(i).Fold=diffStruct.FoldChanges(i);
end


    
end


