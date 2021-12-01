


%%%%%%% this script is for the FMR1 loomhabituation
%%%%%%% data. it is specially to get the calcium traces from the CNMF 
%%%%%%% segmentation and then zscore them. Also getting the FishID and plane 
%%%%%%% of all the ROIs together including hets, fmr1 and wild types. 
 

MatFiles=dir('*analysis_matlab.mat'); %%to get the files
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
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but  +1 because python indexing starts at 0 and matlab at 1
GoodCalcium=Calcium(Fitness,:);  %%%to combine the Calcium and Fitness variables. Fitness here is the variable were we take the good calcium responses and pairthem with their index number.
%GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
%GoodDF=DF(Fitness,:);


MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFiles called GoodNumber the size of Fitness.
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
F=F.idx_components+1;%%%because indexing in python is from 0 and matlab is at 1
%D=load(name, 'dFonF');
%D=D.dFonF;
GC=C(F,:);
%GS=S(F,:);
%GD=D(F,:);

    Noise=vertcat(Noise,N);
    GN=N(F,:);
    %Calcium=vertcat(Calcium,C);
    %DF=vertcat(DF,D);
    %Spikes=vertcat(Spikes,S);
    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC); 
    GoodNoise=vertcat(GoodNoise,GN);
    %GoodDF=vertcat(GoodDF,GD);
    %GoodSpikes=vertcat(GoodSpikes,GS);

%MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;%%%to get rid of vairables we will not use anymore

%%

GoodCalNoise=zeros(size(GoodNoise));
GoodCalNoise(:,:)=GoodCalcium+GoodNoise;


ZS_CN=zscore(GoodCalNoise,1,2); %%%to normalize the data
ZS_CN=detrend(ZS_CN')';%%% to Remove a linear trend from ZS 


clear Calcium Fitness GoodCalcium Noise GoodCalNoise

 save('s20_fmr1_loomhab_CN.mat','-v7.3');
 
 
 %%
 
 %%%the following is to create 1 column variables with the number of the
%%%fish and the slice number
Numbers=[0 [MatFiles.GoodNumber]]; %%%to take make a vector with the GoodNumber field
counter=1;
idx_Plane=nan(length(ZS_CN),1);%%%% to make an empty (with nans) one column variable the size of ZS
idx_Fish=nan(length(ZS_CN),1);%%%% to make an empty (with nans) one column variable the size of ZS
name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1}); %%%to get the number of the plane    
    idx_Plane(Numbers(i)+1:Numbers(i+1))=Plane; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Plane   
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=str2num(Fish{1}{1}); %%%to get the number of the fish 
    idx_Fish(Numbers(i)+1:Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter %%%to get rid of vairables we will not use anymore

 %%%% NOTE: We corrected a small mistake on the generation of idx_Plane
 %%%% and idx_Fish that shifted 1 ROI per fish. So the correct idx_Fish is saved below.   
 
 save('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');
 
%% Exploring responses to looms

%%%% In this case due to the high number of ROIs I
%%%% will try a first linear regression with the s20 regressors
%%%% with a low threshold to filter noisy ROIs. 

 
load('rawregressS20.mat')

for i=1:size(rawregress,1)
    figure;
    plot(rawregress(i,:))
end

rawregress=rawregress(1:5,1:904); %% Because this movies were shorter

ModelResults=[];
parfor i=1:size(ZS_CN,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=stepwiselm(rawregress',ZS_CN(i,:),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
    
    %%%this is to put the results in the ModelResulsts variable
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom=[ModelResults.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them
idx_rsq1=find(rsquare_loom>0.1 & rsquare_loom<1); %%%then select the rsquare that are between 0.1 and 1
figure; %%%to plot them in a raster plot
imagesc(ZS_CN(idx_rsq1,:), [-0.5 4]);colormap hot
figure; histogram(rsquare_loom);

save('s20_fmr1_loomhab_CN_post_rsq01.mat','idx_rsq1','Numbers','idx_Fish','idx_Plane','ModelResults','rawregress','rsquare_loom','-v7.3');
 

%% kmeans to look for types of responses
%%%here i am doing a kmeans to all the ROIs together
options = statset('UseParallel',1); [idxKmeans_ZS_CN Cmap_ZS_CN]=kmeans(ZS_CN,50,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_ZS_CN1,GoodBetas_ZS_CN1]=Test_Regress(Cmap_ZS_CN,Cmap_ZS_CN,idxKmeans_ZS_CN,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3


% clear Calcium Fitness GoodCalcium
% 
 save('s20_postKmeans_CN.mat','idxKmeans_ZS_CN','Cmap_ZS_CN','Model_ZS_CN1','GoodBetas_ZS_CN1');
 
 %%
%%%here i am doing a kmeans to get some raw regressors but with the ROIs
%%%that passed the r2=0.1 filter
options = statset('UseParallel',1); [idxKmeans_ZS_CN_rsq01 Cmap_ZS_CN_rsq01]=kmeans(ZS_CN(idx_rsq1,:),50,'Options',options,'Distance','cityblock','Replicates',3,'MaxIter',1000,'Display','final');
[Model_ZS_CN2,GoodBetas_ZS_CN2]=Test_Regress(Cmap_ZS_CN_rsq01,Cmap_ZS_CN_rsq01,idxKmeans_ZS_CN_rsq01,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3


% clear Calcium Fitness GoodCalcium
% 
 save('s20_postKmeans_CN_rsq01.mat','idxKmeans_ZS_CN_rsq01','Cmap_ZS_CN_rsq01','Model_ZS_CN2','GoodBetas_ZS_CN2');
 
 
 %%
 
 %%%%% to get the corrected idx_Fish
 
 idx_Fish=nan(length(ZS_CN),1);
 name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one	
    name=strcat(MatFiles(i).name);
    
    [name2,~]=regexp(name,'loomhab_(\d+)_','tokens','match'); %%%to get the number of the fish 
    [name3,~]=regexp(name,'fish(\d+)_','tokens','match'); %%%to get the number of the fish
    
    Fish=strcat(name2{1}{1},name3{1}{1});Fish=str2double(Fish); %%%to get the number of the fish 
    idx_Fish(Numbers(i)+1:Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter %%%to get rid of vairables we will not use anymore

unique(idx_Fish)

  save('s20_good_idx_Fish.mat','idx_Fish');
 
 %%
 %%%%% finding the goodclusters to use
 
 %%%% getting the old inhib regressor
 load('inhib_s20_regress_200CL.mat','inhib_s20_regress');

 rawregress(6,:)=inhib_s20_regress(1,1:904);
 
 
 [Model_ZS_CN3,GoodBetas_ZS_CN3]=Test_Regress(rawregress,Cmap_ZS_CN_rsq01,idxKmeans_ZS_CN_rsq01,0.1);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

 
 %% to visualize the results
 
 idxKmeans_rawNrsq=idxKmeans_ZS_CN(idx_rsq1);
 
 
 GoodBetas_ZS_CN3=[13 21 29 32 35 36 38 40 43 44];
 
 
 GoodClust_goodmembers=[];Threshold=0.5;
for i=1:length(GoodBetas_ZS_CN3)    
    ZS_temp=ZS_CN(idx_rsq1(find(idxKmeans_rawNrsq==GoodBetas_ZS_CN3(i))),:);
    corr_temp=zeros(1,length(ZS_temp));
    parfor jj=1:size(ZS_temp,1)
        temp_corr=corrcoef(Cmap_ZS_CN(GoodBetas_ZS_CN3(i),:), ZS_temp(jj,:));
        corr_temp(jj)=temp_corr(1,2);
    end    
    GoodClust_goodmembers(i).ZS=ZS_temp(find(corr_temp>Threshold),:);
    idx_temp=idx_rsq1(find(idxKmeans_rawNrsq==GoodBetas_ZS_CN3(i)));
    GoodClust_goodmembers(i).idx=idx_temp(find(corr_temp>Threshold));
    GoodClust_goodmembers(i).mean=mean(GoodClust_goodmembers(i).ZS,1);
    GoodClust_goodmembers(i).STD=std(GoodClust_goodmembers(i).ZS,1,1);       
end
  

idx_Fish_cat=categorical(idx_Fish);


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1800, 900]);
rows=length(GoodBetas_ZS_CN3);counter=1;
for i=1:size(GoodBetas_ZS_CN3,2)
    
    subplot(rows,4,counter);plot(GoodClust_goodmembers(i).mean); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(GoodClust_goodmembers(i).ZS,[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(GoodClust_goodmembers(i).idx)); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_cat(GoodClust_goodmembers(i).idx));% ax=gca; ax.FontSize=4;%%% for the fish location
    counter=counter+4;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1800, 900]);
rows=length(GoodBetas_ZS_CN3);counter=1;
for i=1:size(GoodBetas_ZS_CN3,2)
    
    subplot(rows,2,counter);plot(GoodClust_goodmembers(i).mean); %%%to plot the mean
    %subplot(rows,3,counter+1);imagesc(GoodClust_goodmembers(i).ZS,[0 3]);%%%for the raster plot
    %subplot(rows,4,counter+2);histogram(idx_Plane(GoodClust_goodmembers(i).idx)); %%%for the plane location   
    subplot(rows,2,counter+1);histogram(idx_Fish_cat(GoodClust_goodmembers(i).idx)); ax=gca; ax.FontSize=4;%%% for the fish location
    counter=counter+2;
    
end
 




