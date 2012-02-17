function [data] = illumina_methylation_control_samp1(data)
%Illumina methylation controls

close all;

%First - find the samples which correspond to the methylation controls
control_samples=data.sample_class==1;

%Extract Methylation percentages(real)
data.actual_meth=double(data.sample_other(control_samples));

%Extract Data from Control Samples

data.meas_meth=data.avg(:,control_samples);

for i=1:size(data.meas_meth,1)
  figure('Visible','off');
  [data.gene_fit{i} data.gof_fit{i}]=fit(data.actual_meth,data.meas_meth(i,:)', ...
'poly1');
  plot(data.actual_meth,data.meas_meth(i,:),'o','MarkerSize', 14);
  hold all;
  %For some reason, plotting a cfit object doesn't work if you put
  %LineWidth in there, so change after . . 
  h=plot(data.gene_fit{i});
  set(h, 'LineWidth', 1.5);
  print('-dpng', ['Movie/Controls/' data.genes{i} '.png']);
  close all;
end

                  

