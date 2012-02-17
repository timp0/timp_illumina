function [] = clustergram_plot2(namey,huge)
%This function programmatically gets the handle for the clustergram and
%prints it.  

%Make all handles visible. This is necessary because clustergram
%objects are created with 'HandleVisibility' property set to 'off'.
set(0,'ShowHiddenHandles','on') 
%Get all handles from root
allhnds = get(0,'Children');
%Find the handles that correspond to clustergram objects
cgfigidx = strmatch('Clustergram',get(allhnds,'Tag'));
cffighnd = allhnds(cgfigidx);
%Make the non-visible handles hidden again
set(0,'ShowhiddenHandles','off')

%Select only the last clustergram
if length(cffighnd)>1
    warning('More than one clustergram open!!');
    cffighnd = cffighnd(end);
end
if huge
    set(cffighnd, 'paperSize', [20 20], 'paperPosition', [.25 .25 19.75 19.75]);
end

print(cffighnd, '-dpdf', ['Movie\' namey '.pdf']);
print(cffighnd, '-painters','-dpdf', ['Movie\' namey '_ai.pdf']);







