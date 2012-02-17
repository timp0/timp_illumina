function [] = plot_tumornormal3(data,modif)
%Plot the tumor normal pairs

pairnames={'Blah', 'Breast', 'Colon', 'Lung', 'Ovary', 'Wilms', 'Thyroid', 'Pancreas'};
opty.huge=1;

for i=2:8
    %Extract the tumor normal pairs using sample_class index
    opty.classes=[i];
    opty.modif=[modif '_' pairnames{i}];
    classcompare_sub3(data,opty);
    
end

end

