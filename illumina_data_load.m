function [data] = illumina_data_load2()
%Ok = the point of this is to take the 3 perl generated files and load them
%into MATLAB structure data

%Base Filename
base='Feinberg_CuMeth_Background_Data_Table';


%First - samples
samp_file=fopen([base '_mat_sample.csv']);
sampy=textscan(samp_file,'%s %d %d %d', 'delimiter', ',');
data.samples=sampy{1};
data.samp_num=size(data.samples,1);
data.sample_class=sampy{2};
data.sample_id=sampy{3};
data.sample_other=sampy{4};
fclose(samp_file);

%Second - targets
targ_file=fopen([base '_mat_targer.csv']);
targy=textscan(targ_file,'%s %s','delimiter', ',');
data.genes=targy{1};
data.targets=targy{2};
fclose(targ_file);


data.raw=csvread([base '_mat_data.csv']);

for i=1:data.samp_num
    data.avg(:,i)=data.raw(:,-4+i*5);
    data.cy3(:,i)=data.raw(:,-3+i*5);
    data.cy5(:,i)=data.raw(:,-2+i*5);
end

end
