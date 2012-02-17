function [] = pvalue_table_csv(good_info, namey)

fid=fopen(['Movie/' namey '.csv'], 'w');
fprintf(fid, 'Location_Name, FDR, pValue, Fold_Change\n');
for i=1:size(good_info,2)
    fprintf(fid, '%s,%g,%g,%g\n',good_info(i).gene, good_info(i).FDR, good_info(i).pValues, ...
        good_info(i).Fold);
end
