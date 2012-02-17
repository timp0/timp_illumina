function [diff_geneidx,goody] = illumina_pvalue3(data,opty)
%Get significant genes and p-values, fdrs for them

if nargin < 2
    opty.not=1;
end

if ~isfield(opty, 'modif')
    opty.modif='a';
end

if ~isfield(opty,'tumor_class')
    opty.tumor_class=[2 4 6 8 10];
end

if ~isfield(opty,'normal_class')
    opty.normal_class=[3 5 7 9 11];
end

if ~isfield(opty,'subset')
    opty.subset=true(length(data.genes.label),1);
end

if ~isfield(opty,'id_set')
    opty.id_set=unique(data.sample_id);
end





tumor_data=data.avg(opty.subset,(ismember(data.sample_class,opty.tumor_class)&ismember(data.sample_id,opty.id_set)));
normal_data=data.avg(opty.subset,(ismember(data.sample_class,opty.normal_class)&ismember(data.sample_id,opty.id_set)));



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
plot(aDiff, -log10(pValues), '.','MarkerSize', 14);
%Get edges of axes
xLims=get(gca, 'XLim');
yLims=get(gca, 'YLim');
hold all
plot(xLims, [-log10(sig_p) -log10(sig_p)], '-.r', 'LineWidth', 1.5);
plot([-sig_dif -sig_dif], yLims, '-.r', 'LineWidth', 1.5);
plot([sig_dif sig_dif], yLims, '-.r', 'LineWidth', 1.5);

print('-dpdf', ['Movie/' opty.modif '_volcano.pdf']);

close

subplot(2,1,1),cdfplot(-log10(pValues));
xlabel('-log(p)')
set(gca,'Xlim', [0 10])

subplot(2,1,2), hist(-log10(pValues))

hold all
[a,b]=ksdensity(-log10(pValues))
plot(b,a*384,'k', 'LineWidth', 2)
xlabel('-log(p)')
set(gca, 'Xlim', [0 10])

print('-dpdf', ['Movie/' opty.modif '_phisto.pdf']);

close all;

diff_geneidx=((aDiff<=-.17)|(aDiff>=.17))&(pValues<=.01);

goody.gene=data.genes.label(diff_geneidx);
goody.FDR=pFDR(diff_geneidx);
goody.pValues=pValues(diff_geneidx);
goody.aDiff=aDiff(diff_geneidx);

    
end


