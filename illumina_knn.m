function [tester] = illumina_knn1(data, opty)
%Function to do a k-nearest neighbor comparison of the classes and get
%class predictors


if nargin < 2
    opty.not=1;
end

if ~isfield(opty, 'cluster')
    opty.cluster=3;
end

if ~isfield(opty, 'modif')
    opty.modif='a';
end

if ~isfield(opty, 'huge')
    opty.huge=0;
end

if ~isfield(opty,'classes')
    opty.classes=[2 3; 4 5; 6 7; 8 9; 10 11];
end

if ~isfield(opty,'subset')
    opty.subset=true(length(data.genes.label),1);
end

if ~isfield(opty,'id_set')
    opty.id_set=unique(data.sample_id);
end




%Get samples

%Make the sample_class index
sampy=ismember(data.sample_class,opty.classes)&ismember(data.sample_id,opty.id_set);

sample_names=data.samples(sampy);
sample_dat=data.avg(opty.subset,sampy);
sample_num=length(sample_names);

%Define different classes
ori_classes=data.sample_class(sampy);
sample_classes=zeros(sample_num,1);

%Get different clusters of classes(multiple tumor types sometimes)
for i=1:size(opty.classes,2)
    sample_classes=sample_classes+i*ismember(ori_classes, opty.classes(:,i));
end




%Get indicices for cross validation - K-fold is the best answer as
%it will allow for leave-one-out in a organized, instead of random
%way with looping (disjointed data sets is the terminology used)
which=crossvalind('Kfold', sample_num,sample_num);


%Got to loop through and leave out each sample for thorough test
tester = classperf(sample_classes);
for i=1:sample_num
  test = (which == i); train= ~test;
  classed = knnclassify(sample_dat(:,test)',sample_dat(:,train)', ...
                        sample_classes(train),length(opty.classes));
  classperf(tester,classed,test);
end
