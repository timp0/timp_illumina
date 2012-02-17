function [diff_geneidx,goody] = illumina_pvalue2(genes,namey,tumor_data,normal_data)
%Get significant genes and p-values, fdrs for them

%Significant difference set
sig_dif=0.17;

%Sig p-value set
sig_p=.01;

nTumor=size(tumor_data,2);
nNormal=size(normal_data,2);


meanTumor=mean(tumor_data,2);
meanNormal=mean(normal_data,2);



pValues = mattest(tumor_data, normal_data, 'Permute', 1e5);

aDiff = meanTumor-meanNormal;

[pFDR, qValues] = mafdr(pValues);

close all
%Make volcano-like plot of methylation
plot(aDiff, -log10(pValues), '.');
%Get edges of axes
xLims=get(gca, 'XLim');
yLims=get(gca, 'YLim');
hold all
plot(xLims, [-log10(sig_p) -log10(sig_p)], '-.r', 'LineWidth', 1.5);
plot([-sig_dif -sig_dif], yLims, '-.r', 'LineWidth', 1.5);
plot([sig_dif sig_dif], yLims, '-.r', 'LineWidth', 1.5);

print('-dpdf', ['Movie/' namey '_volcano.pdf']);

close all;

diff_geneidx=((aDiff<=-.17)|(aDiff>=.17))&(pValues<=.01);

goody.gene=genes(diff_geneidx);
goody.FDR=pFDR(diff_geneidx);
goody.pValues=pValues(diff_geneidx);
goody.aDiff=aDiff(diff_geneidx);

    
end


