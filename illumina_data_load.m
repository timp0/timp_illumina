function [data] = illumina_data_load3(base)
%Ok = the point of this is to take the 3 perl generated files and load them
%into MATLAB structure data

if nargin < 1
    %Base Filename
    base='New_norm';
end


%First - samples
samp_file=fopen([base '_mat_sample.csv']);
sampy=textscan(samp_file,'%s %d %d %d', 'delimiter', ',','HeaderLines', 1);
data.samples=sampy{1};
data.samp_num=size(data.samples,1);
data.sample_class=sampy{2};
data.sample_id=sampy{3};
data.sample_other=sampy{4};
fclose(samp_file);

%Second - targets
targ_file=fopen([base '_mat_targer_ucscisl_hmmisl.csv']);
targy=textscan(targ_file,'%s %s %s %d %d %d %d %d %d %d %d %d %f %f %d %d %d %d %d %f %f','delimiter', ',','HeaderLines', 1);
data.genes.label=targy{1};
data.genes.illumina=targy{2};
data.genes.chromy=targy{3};
data.genes.loc=[targy{4} targy{5}];
data.genes.good=targy{6};
data.genes.region=targy{7};
data.genes.ucsc.rel=targy{8};
data.genes.ucsc.dist=targy{9};
data.genes.ucsc.loc=[targy{10} targy{11}];
data.genes.ucsc.numcpg=targy{12};
data.genes.ucsc.percpg=targy{13};
data.genes.ucsc.obsexp=targy{14};
data.genes.hmm.rel=targy{15};
data.genes.hmm.dist=targy{16};
data.genes.hmm.loc=[targy{17} targy{18}];
data.genes.hmm.numcpg=targy{19};
data.genes.hmm.percpg=targy{20};
data.genes.hmm.obsexp=targy{21};
fclose(targ_file);




raw=dlmread([base '_mat_data.csv'],',',1,1);

data.avg=raw(:,1:5:end-1);
data.cy3=raw(:,2:5:end-1);
data.cy5=raw(:,3:5:end-1);


end
