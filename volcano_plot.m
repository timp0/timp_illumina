function [] = volcano_plot(namey,huge)
%This function programmatically gets the handle for the Volcano and
%prints it.  

%Make all handles visible. This is necessary because clustergram
%objects are created with 'HandleVisibility' property set to 'off'.
set(0,'ShowHiddenHandles','on') 
%Get all handles from root
allhnds = get(0,'Children');
%Find the handles that correspond to clustergram objects
cgfigidx = strmatch('Volcano',get(allhnds,'Tag'));
cffighnd = allhnds(cgfigidx);
%Make the non-visible handles hidden again
set(0,'ShowhiddenHandles','off')

%Select only the last Volcano plot
if length(cffighnd)>1
    warning('More than one volcano open!!');
    cffighnd = cffighnd(end);
end
if huge
    set(cffighnd, 'paperSize', [20 20], 'paperPosition', [.25 .25 19.75 19.75]);
else
    set(cffighnd, 'paperSize', [8.5 11], 'paperPosition', [.25 .25 6.25 4.25]);
end

print(cffighnd, '-dpdf', ['Movie/' namey '_volcano.pdf']);







