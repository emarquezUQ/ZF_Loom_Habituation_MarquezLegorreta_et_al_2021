
%%%%% this script is to collect the CNMF results of f60 group but this time
%%%%% adding the noise filtered out by the CNMF

%%% First it is taking the resutls from the CNMF (caiman). 
%%% Then it does a first exploration of the dataset with a kmeans to see
%%% possible types of respones and their distribution across fish and
%%% planes. 
%%% Finally we tested with a second kmeans  
%%% with a higher number of clusters to see
%%% if other respones were discovered. 


MatFiles=dir('*f60*analysis_matlab.mat'); %%to get the files
name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium

Noise=load(name, 'Noise');
Noise=Noise.Noise;

Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,:);  %%%to combine the Calcium and Fitness variables (need to ask Gilles what Fitness is). Fitness here is the variable were we take the good calcium responses from the HPC analysis and pairthem with their index number.

GoodNoise=Noise(Fitness,:);

GoodCalNoise=zeros(size(GoodNoise));
GoodCalNoise(:,:)=GoodCalcium+GoodNoise;

MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFilesCalcium called GoodNumber the size of Fitness.
for i = 2:length(MatFiles) %%%%to take the slices one by one starting by the second one cause we already did this with the first one
    %%%% we are going to do the same thing that before but for all the
    %%%% slices
name=strcat(MatFiles(i).name);%%%%to take the name of the slice in turn
C=load(name, 'DenoisedTraces');%%to load only the DenoisedTraces from the file
C=C.DenoisedTraces;%%%% <-- take the field called DenoisedTraces from the C structure and make it the new C

N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;%%%because indexing in python is from 0 and matlab is at 1

GC=C(F,:);

    GN=N(F,:);
    
    GCN=zeros(size(GN));
    GCN(:,:)=GC+GN;
    
    GoodCalNoise=vertcat(GoodCalNoise,GCN);
    
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;%%%to get rid of vairables we will not use anymore

clear Calcium Fitness Noise

% save('f60_CN_matfiles.mat','-v7.3');
%%

% GoodCalNoise=zeros(size(GoodNoise));
% GoodCalNoise(:,:)=GoodCalcium+GoodNoise;

clear   GoodCalcium GoodNoise 

  save('f60_CN_GoodCalNoise.mat','-v7.3');
 

ZS_CN=zscore(GoodCalNoise,1,2); %%%to normalize the data
ZS_CN=detrend(ZS_CN')';%%% to Remove a linear trend from ZS (why???)

clear GoodCalNoise

 save('f60_CN.mat','-v7.3');
 
%%
%%%here i am doing a kmeans to get some raw regressors
options = statset('UseParallel',1); [idxKmeans_ZS_CN Cmap_ZS_CN]=kmeans(ZS_CN,50,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS_CN2,GoodBetas_ZS_CN2]=Test_Regress(Cmap_ZS_CN,Cmap_ZS_CN,idxKmeans_ZS_CN,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

saveas(gcf,'50rawclusters_CN_f60','jpg');


% clear Calcium Fitness GoodCalcium
% 
 save('f60_postKmeans_CN.mat','-v7.3');
 
 %%
%%% i will be trying to get the raw clusters again by puting more clusters
 
[idxKmeans_ZS_CN_500CL Cmap_ZS_CN_500CL]=kmeans(ZS_CN,500,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
%[Model_ZS_CN_500CL,GoodBetas_ZS_CN2_500CL]=Test_Regress(Cmap_ZS_CN_500CL,Cmap_ZS_CN_500CL,idxKmeans_ZS_CN_500CL,0.1);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

%saveas(gcf,'500rawclusters_CN_f60','jpg');

counter2=0;
for i=1:10
    Fighandle=figure;counter=1;
    set(Fighandle, 'Position', [100, 100, 1500, 1000]);
    for j=1:50
        subplot(7,8,counter);plot(Cmap_ZS_CN_500CL(j+counter2,:));
        counter=counter+1;
    end
    counter2=counter2+50;
    print(Fighandle,strcat('500rawclusters_CN_f60_',num2str(i)),'-dpng','-r0');
    close all;
end

 
 save('f60_postKmeans_CN_500CL.mat','-v7.3');
 
