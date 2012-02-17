function [data] = illumina_methylation_control_samp3(data)
%Illumina methylation controls

close all;

%First - find the samples which correspond to the methylation controls
control_samples=data.sample_class==1;

%Extract Methylation percentages(real)
actual_meth=double(data.sample_other(control_samples));

%Extract Data from Control Samples

data.adj=data.avg;

meas_meth=data.avg(:,control_samples);




for i=1:size(meas_meth,1)    
    figure('Visible','off');
    [data.genes.fit{i} data.genes.gof_fit{i}]=fit(actual_meth,meas_meth(i,:)', ...
        'poly1');
    plot(actual_meth,meas_meth(i,:),'o','MarkerSize', 14);
    hold all;
    %For some reason, plotting a cfit object doesn't work if you put
    %LineWidth in there, so change after . .
    h=plot(data.genes.fit{i});
    set(h, 'LineWidth', 1.5);
    
    coefs=coeffvalues(data.genes.fit{i});
    if coefs(1)<0
        data.neg_slope(i)=true;
    else
        data.neg_slope(i)=false;
    end
    
    print('-dpng', ['Movie/Controls/' data.genes.label{i} '.png']);    
    close all;
    
    data.adj(i,:)=data.genes.fit{i}(data.avg(i,:));
    
end

end
