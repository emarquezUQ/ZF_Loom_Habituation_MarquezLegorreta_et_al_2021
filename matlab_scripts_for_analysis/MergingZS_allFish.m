
%%%%% this script is to get all the ROIs from the 4 wild type datasets and
%%%%% zscore them together. It is needed to get them raw and trim them as the have different lenghts.

%%%% s60 is too big and matlab protested. sowei had to trim it while
%%%% getting the raw data. 

load('All_More_BrainReg2.mat');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');


load('final_F20_step1.mat','idx_Fish_f20');

load('final_F60_step1_2.mat','idx_Fish_f60','ZS_short_F60');

load('All_means_maxResp_CL4_normalized.mat','S_trim');

load('final_S20_step1.mat','ZS_s20','idx_Fish_s20');

load('final_S60_step1.mat','idx_Fish_s60','ZS_short_S60');

%%%%% for s60. 
S60data_CN=load('s60_postKmeans_CN.mat','MatFiles');
MatFiles=S60data_CN.('MatFiles');



MatFiles=dir('*s60*analysis_matlab.mat'); %%to get the files
name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium
%MatFiles(1).number=size(Calcium,1);
%Spikes=load(name, 'Spikes');
%Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
%DF=load(name, 'dFonF');
%DF=DF.dFonF;
Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,ZS_short_S60(S_trim));  %%%to combine the Calcium and Fitness variables (need to ask Gilles what Fitness is). Fitness here is the variable were we take the good calcium responses from the HPC analysis and pairthem with their index number.
%GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,ZS_short_S60(S_trim));
%GoodDF=DF(Fitness,:);


MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFilesCalcium called GoodNumber the size of Fitness.
for i = 2:length(MatFiles) %%%%to take the slices one by one starting by the second one cause we already did this with the first one
    %%%% we are going to do the same thing that before but for all the
    %%%% slices
name=strcat(MatFiles(i).name);%%%%to take the name of the slice in turn
C=load(name, 'DenoisedTraces');%%to load only the DenoisedTraces from the file
C=C.DenoisedTraces;%%%% <-- take the field called DenoisedTraces from the C structure and make it the new C
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
% S=load(name, 'Spikes');
% S=S.Spikes;
N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;%%%because indexing in pythong is from 0 and matlab is at 1
%D=load(name, 'dFonF');
%D=D.dFonF;
GC=C(F,ZS_short_S60(S_trim));
%GS=S(F,:);
%GD=D(F,:);

    %Noise=vertcat(Noise,N);
    GN=N(F,ZS_short_S60(S_trim));
    %Calcium=vertcat(Calcium,C);
    %DF=vertcat(DF,D);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC); %The fish 20+ are longer
    GoodNoise=vertcat(GoodNoise,GN);
    %GoodDF=vertcat(GoodDF,GD);
    %GoodSpikes=vertcat(GoodSpikes,GS);

%MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N; %%%to get rid of vairables we will not use anymore

GoodCalNoise_s60=zeros(size(GoodNoise));
GoodCalNoise_s60(:,:)=GoodCalcium+GoodNoise;

  
clear Calcium Fitness GoodCalcium Noise GoodNoise 



%% triming


GoodCalNoise_s20=GoodCalNoise_s20_long(:,S_trim);

%%% to check it worked
figure;plot(mean(GoodCalNoise_s20(s20_cleaned_idxs.clust_s20_CL4_cleaned.clust_s20_CL4_7_cleaned,:)));

clear GoodCalNoise_s20_long


GoodCalNoise_f60=GoodCalNoise_f60_long(:,ZS_short_F60);
figure;plot(mean(GoodCalNoise_f60(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_7_cleaned,:)));

clear GoodCalNoise_f60_long

save('GoodCalNoise_all.mat','GoodCalNoise_f20','GoodCalNoise_f60','GoodCalNoise_s20','GoodCalNoise_s60','-v7.3');


GoodCalNoise=vertcat(GoodCalNoise_f20,GoodCalNoise_f60,GoodCalNoise_s20,GoodCalNoise_s60);

clear GoodCalNoise_f20 GoodCalNoise_f60 GoodCalNoise_s20 GoodCalNoise_s60

ZS_CN=zscore(GoodCalNoise,1,2); %%%to normalize the data
ZS_CN=detrend(ZS_CN')';%%% to Remove a linear trend from ZS 

clear GoodCalNoise

%%% puting the fish idx together too. 
idx_Fish_all=vertcat(idx_Fish_f20,idx_Fish_f60,idx_Fish_s20,idx_Fish_s60);



%% adjusting idx

%%% f20 doesnt need it cause it is at the begining. 
figure;plot(mean(ZS_CN(f20_cleaned_idxs.clust_f20_CL4_cleaned.clust_f20_CL4_7_cleaned,:)));



idx_f60_adjust=(length(idx_Fish_f20)+1:length(idx_Fish_f20)+length(idx_Fish_f60));

%%% to see if it works.
figure;plot(mean(ZS_CN(idx_f60_adjust(f60_cleaned_idxs.clust_f60_CL4_cleaned.clust_f60_CL4_7_cleaned),:)));


idx_s20_adjust=(length(idx_Fish_f20)+length(idx_Fish_f60)+1:length(idx_Fish_f20)+length(idx_Fish_f60)+length(idx_Fish_s20));

%%% to see if it works.
figure;plot(mean(ZS_CN(idx_s20_adjust(s20_cleaned_idxs.clust_s20_CL4_cleaned.clust_s20_CL4_7_cleaned),:)));


idx_s60_adjust=(length(idx_Fish_f20)+length(idx_Fish_f60)+length(idx_Fish_s20)+1:length(idx_Fish_f20)+length(idx_Fish_f60)+length(idx_Fish_s20)+length(idx_Fish_s60));

%%% to see if it works.
figure;plot(mean(ZS_CN(idx_s60_adjust(s60_cleaned_idxs.clust_s60_CL4_cleaned.clust_s60_CL4_7_cleaned),:)));



save('ZS_N_Fish_all.mat','ZS_CN','idx_Fish_all','idx_f60_adjust','idx_s20_adjust','idx_s60_adjust','-v7.3');



