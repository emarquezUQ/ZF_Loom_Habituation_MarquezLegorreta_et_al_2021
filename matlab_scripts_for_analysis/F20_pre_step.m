%%%%% this script is to collect the CNMF results of f20 group 

%%% First it is taking the resutls from the CNMF (caiman). 
%%% Then it does a first exploration of the dataset with a kmeans to see
%%% possible types of respones and their distribution across fish and
%%% planes. We check how they look, which slices and which fish
%%% they are in and then select the ones that are in most fish.
%%% We then used a linear regression to look for the ROIs with the selected
%%% proiles and clean them from noisy responses. 
%%% Finally we tested with a second kmeans to the filtered datasets to see
%%% if other respones were discovered. 

%%% Note: these are denoised traces, we found later that adding the noise
%%% allowed us to detect the inhibited responses and details on other
%%% profiles. But this first analysis was useful to know the main type of
%%% responses and generate regressors. 

MatFiles=dir('*analysis_matlab.mat'); %%to get CNMF output files
name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium

Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness. +1 Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,:); 

MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFilesCalcium called GoodNumber the size of Fitness.
for i = 2:length(MatFiles) %%%%to take the slices one by one starting by the second one cause we already did this with the first one
    %%%% we are going to do the same thing that before but for all the
    %%%% slices
name=strcat(MatFiles(i).name);%%%%to take the name of the slice in turn
C=load(name, 'DenoisedTraces');%%to load only the DenoisedTraces from the file
C=C.DenoisedTraces;%%%% <-- take the field called DenoisedTraces from the C structure and make it the new C

F=load(name, 'idx_components');
F=F.idx_components+1;%%%because indexing in python is from 0 and matlab is at 1

GC=C(F,:);

    Fitness=horzcat(Fitness,F);
    GoodCalcium=vertcat(GoodCalcium,GC); 
    
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N; 
ZS=zscore(GoodCalcium,1,2); %%%to normalize the data
ZS=detrend(ZS')';%%% to Remove a linear trend 



%%%here i am doing a kmeans to get some raw regressors
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,50,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_ZS,Cmap_ZS,idxKmeans_ZS,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

goodmaps=[14 17 31 47 50];
for i=1:length(goodmaps)
rawregress(i,:)=Cmap_ZS(goodmaps(i),:);
end

save ('rawregressF20_2.mat','rawregress');

for i=1:length(goodmaps)
figure;
plot(rawregress(i,:));
end


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(goodmaps);

for i=1:length(goodmaps)
   
    subplot(4,2,counter);plot(rawregress(i,:));
     counter=counter+1;
end


%%% to get the fishID and plane of each ROI 
Numbers=[1 [MatFiles.GoodNumber]]; 
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
for i=1:length(MatFiles) 	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1}); %%%to get the number of the plane    
    idx_Plane(Numbers(i):Numbers(i+1))=Plane; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Plane   
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=str2num(Fish{1}{1}); %%%to get the number of the fish 
    idx_Fish(Numbers(i):Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter 

%%%% to run a preliminary linear regression to the clusters found
ModelResults=[];
parfor i=1:size(ZS,1)  
    
    
    mdl=fitlm(rawregress',ZS(i,:));
    
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
end
rsquare_loom=[ModelResults.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them
idx_rsq=find(rsquare_loom>0.5 & rsquare_loom<1); %%%then select the rsquare that are between 0.5 and 1
figure; %%%to plot them in a raster plot
imagesc(ZS(idx_rsq,:), [-0.5 4]);colormap hot

ZSrsq05=ZS(idx_rsq,:);


%%%%% getting cofficients
coefficients={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(ModelResults)
    coef=[ModelResults(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    if ~isempty(temp)
        
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            if coef.pValue(coef_idx)<0.05%%%to select the coef that are bellow the p value we want, in this case 0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients); 
coefficients(idxempty)={0}; 
clearvars idxempty idx coef_idx coef  
coefficients=cell2mat(coefficients); 


coefficients_Nb=coefficients>0; 
coefficients_Nb=sum(coefficients_Nb,2);
idx_coef=find(coefficients_Nb>0);
idx_coef_rsq=intersect(idx_rsq,idx_coef); 

%%%now we run another kmeans with the filtered data set idx_coef_rsq. first
%%%asking for 5 clusters 
options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS(idx_coef_rsq,:),5,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_ZS_rsq,rawregress,idxKmeans_ZS_rsq,0.4);%%%here we do another linear regression and we select the ones with an r2 value above 0.4

idxKmeans1_coef_rsq=idxKmeans_ZS_rsq(idx_coef_rsq);%%% to make a variable with the common idx of the interesting clusters from the 1st kmeans and the traces hat passed the regression
rows=length(GoodBetas_ZS2);
counter=1;


GoodBetas=[GoodBetas_ZS2];
GoodBetas=GoodBetas_ZS2(GoodBetas);
rows=length(GoodBetas);
counter=1;

%%%this is to make a figure where we will plot the mean of clusters
%%%selected the raster plot, the levels and the fish where they are found.
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(GoodBetas);counter=1;
for i=GoodBetas
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans_ZS_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end


Numbers(1)=0; 
colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors
colors = colors*256;


%%%%%this is to do it together with the raster plots and histograms of
%%%%%localization of the clusters in fish and levels
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas
    idx_temp=find(idxKmeans_ZS_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1),'color',colors(counter2,:)/256); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end


clear Calcium Fitness GoodCalcium

 save('D:\Emmanuel\SnF_20vs60ISI_looms_1st_analysis\f20_r2050_CL5.mat','-v7.3');


