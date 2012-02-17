function [bad_slope] = illumina_methylation_control_samp2(data)
%Illumina methylation controls

close all;

%First - find the samples which correspond to the methylation controls
control_samples=data.sample_class==1;

%Extract Methylation percentages(real)
data.actual_meth=double(data.sample_other(control_samples));

%Extract Data from Control Samples



data.meas_meth.avg=data.avg(:,control_samples);
data.meas_meth.cy5=data.cy5(:,control_samples);
data.meas_meth.cy3=data.cy3(:,control_samples);


for i=1:size(data.meas_meth.avg,1)    
    figure('Visible','off');
        [data.gene_fit{i} data.gof_fit{i}]=fit(data.actual_meth,data.meas_meth.avg(i,:)', ...
            'poly1');
        subplot(1,3,1);
        plot(data.actual_meth,data.meas_meth.avg(i,:),'o','MarkerSize', 14);
        hold all;
        %For some reason, plotting a cfit object doesn't work if you put
        %LineWidth in there, so change after . .
        h=plot(data.gene_fit{i});
        set(h, 'LineWidth', 1.5);
        
        coefs=coeffvalues(data.gene_fit{i});
        if coefs(1)<0
            bad_slope(i)=1;
        else
            bad_slope(i)=0;
        end
        
        
        cy3_fit{i}=fit(data.actual_meth, data.meas_meth.cy3(i,:)', 'poly1');
        subplot(1,3,2);
        plot(data.actual_meth, data.meas_meth.cy3(i,:),'o', 'MarkerSize', 14);
        hold all;
        h=plot(cy3_fit{i});
        set(h, 'LineWidth', 1.5);
        
        cy5_fit{i}=fit(data.actual_meth, data.meas_meth.cy5(i,:)', 'poly1');
        subplot(1,3,3);
        plot(data.actual_meth, data.meas_meth.cy5(i,:),'o', 'MarkerSize', 14);
        hold all;
        h=plot(cy5_fit{i});
        set(h, 'LineWidth', 1.5);
        
        
        print('-dpng', ['Movie/Controls/' data.genes{i} '.png']);
        close all;
    end
    
end

