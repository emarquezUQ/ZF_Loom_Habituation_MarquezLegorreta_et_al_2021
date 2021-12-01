%%%%this script is to take a look at the responses of the F60 fish.

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


MatFiles=dir('*f60*analysis_matlab.mat'); %%to get the files
name=strcat(MatFiles(1).name); %%%to get the name of the files
Calcium=load(name, 'DenoisedTraces'); %%to load only the DenoisedTraces from the file, the raw data was denoised by the CNMF (The Cluster Analysis tool calculates clusters based on a Constrained non-negative matrix factorization (NMF) clustering method.)
Calcium=Calcium.DenoisedTraces; %%%% <-- take the field called DenoisedTraces from the Calcium structure and make it the new Calcium

Fitness=load(name, 'idx_components');%%to load only the idx_components from the file, they are based on what a Gcamp spike should be and they will filter the true spikes in our data
Fitness=Fitness.idx_components+1; %%%% <-- take the field called idx_components from the Fitness structure and make it the new Fitness but why +1?? Because python indexing starts at 0 ant matlab at 1
GoodCalcium=Calcium(Fitness,:);  

MatFiles(1).GoodNumber=length(Fitness); %%%% <-- Create a field inside MatFilesCalcium called GoodNumber the size of Fitness.
for i = 2:length(MatFiles) %%%%to take the slices one by one starting by the second one cause we already did this with the first one
    %%%% we are going to do the same thing that before but for all the
    %%%% slices
name=strcat(MatFiles(i).name);%%%%to take the name of the slice in turn
C=load(name, 'DenoisedTraces');%%to load only the DenoisedTraces from the file
C=C.DenoisedTraces;%%%% <-- take the field called DenoisedTraces from the C structure and make it the new C

F=load(name, 'idx_components');
F=F.idx_components+1;%%%because indexing in pythong is from 0 and matlab is at 1

GC=C(F,:);

MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N; 
ZS=zscore(GoodCalcium,1,2); %%%to normalize the data
ZS=detrend(ZS')';%%% to Remove a linear trend from ZS 



%%%here i am doing a kmeans to get some raw regressors
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,50,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_ZS,Cmap_ZS,idxKmeans_ZS,0.3);%%%here we do another linear regression and we select the ones with an r2 value above 0.3

goodmaps=[2 7 15 17 30 33 39 47 50];
for i=1:length(goodmaps)
rawregress(i,:)=Cmap_ZS(goodmaps(i),:);
end

save ('rawregressF60.mat','rawregress');

for i=1:length(goodmaps)
figure;
plot(rawregress(i,:));
end


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(goodmaps);

for i=1:length(goodmaps)
   
    subplot(3,3,counter);plot(rawregress(i,:));
     counter=counter+1;
end

%%%fish and the slice number
Numbers=[1 [MatFiles.GoodNumber]]; %%%to take make a vector with the GoodNumber field
counter=1;
idx_Plane=nan(length(ZS),1);%%%% to make an empty (with nans) one column variable the size of ZS
idx_Fish=nan(length(ZS),1);%%%% to make an empty (with nans) one column variable the size of ZS
name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1}); %%%to get the number of the plane    
    idx_Plane(Numbers(i):Numbers(i+1))=Plane; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Plane   
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=str2num(Fish{1}{1}); %%%to get the number of the fish 
    idx_Fish(Numbers(i):Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter %%%to get rid of vairables we will not use anymore




ModelResults=[];
parfor i=1:size(ZS,1)  %%%parfor is to do it in parallel
    
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

coefficients={}; %%%to make the coefficients variable that we will use. Regression coefficients represent the mean change in the response variable for one unit of change in the predictor variable while holding other predictors in the model constant.
for idx=1:length(ModelResults)%%% to make a variable the size of ModelResults
    coef=[ModelResults(idx).coef];%%%% and then put in another variable the coef field from ModelResults
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');%%%to take the name of the rows of the coef variable
    if ~isempty(temp)%%% if temp is not empty...
        
        for coef_idx=2:height(coef)%%%take the number of rows from coef, except the first one(i think because is the intercept)
            if coef.pValue(coef_idx)<0.05%%%to select the coef that are bellow the p value we want, in this case 0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx); %%%to make an array the size of idx,10 with the coefficient values that were significant
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients); %%%to make a variable with where we will aply in every cell the isempty function wich will help us find the empty places
coefficients(idxempty)={0}; %%% and put a 0 in the places where we found that there were empty cells
clearvars idxempty idx coef_idx coef  %%%clear variables
coefficients=cell2mat(coefficients); %%%to make a matrix of the coefficients array


coefficients_Nb=coefficients>0; %%%to take the values above 0 and make a variable
coefficients_Nb=sum(coefficients_Nb,2);%%% to sum the cells from the second dimension (the looms)
idx_coef=find(coefficients_Nb>0);%%%to take that have more than 2 responses (i might want to try to see what happens if i dont do this)
idx_coef_rsq=intersect(idx_rsq,idx_coef);%%% to make a variable with the common values from idx_rsq and idx_coef, so with the r2 values between 0.3-1 and the coefficients that were significant

%%%now we run another kmeans with the filtered data set idx_coef_rsq. first
%%%asking for 10 clusters and then for 5
options = statset('UseParallel',1); [idxKmeans_ZS_rsq Cmap_ZS_rsq]=kmeans(ZS(idx_coef_rsq,:),5,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS2,GoodBetas_ZS2]=Test_Regress(Cmap_ZS_rsq,rawregress,idxKmeans_ZS_rsq,0.4);%%%here we do another linear regression and we select the ones with an r2 value above 0.4


idxKmeans1_coef_rsq=idxKmeans_ZS(idx_coef_rsq);%%% to make a variable with the common idx of the interesting clusters from the 1st kmeans and the traces hat passed the regression
rows=length(goodmaps);
counter=1;


%%%this is to make a figure where we will plot the mean of clusters
%%%selected the raster plot, the levels and the fish where they are found.
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=length(goodmaps);counter=1;
for i=goodmaps
    %idx_temp=find(idxKmeans1_coef_rsq==i);
    idx_temp=find(idxKmeans1_coef_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    counter=counter+4;
end



%%%% picking up less clusters and doing everything again
%%% i chose the clusters that were in most fish 

goodmaps=[2 15 33 39]; %%% the first kmeans didnt detect a sound response but i am sure it is there. but I will need to sue a regressor
for i=1:length(goodmaps)
rawregress(i,:)=Cmap_ZS(goodmaps(i),:);
end

save ('rawregressF60.mat','rawregress');


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(goodmaps);

for i=1:length(goodmaps)
   
    subplot(3,3,counter);plot(rawregress(i,:));
     counter=counter+1;
end


%%%fish and the slice number
Numbers=[1 [MatFiles.GoodNumber]]; %%%to take make a vector with the GoodNumber field
counter=1;
idx_Plane=nan(length(ZS),1);%%%% to make an empty (with nans) one column variable the size of ZS
idx_Fish=nan(length(ZS),1);%%%% to make an empty (with nans) one column variable the size of ZS
name=strcat(MatFiles(1).name);%%%to get the name of the files (is actually to create the variable name before the loop)
for i=1:length(MatFiles) %%%%to take slices one by one	
    name=strcat(MatFiles(i).name);
    [Plane,~]=regexp(name,'Slice(\d+)_','tokens','match');Plane=str2num(Plane{1}{1}); %%%to get the number of the plane    
    idx_Plane(Numbers(i):Numbers(i+1))=Plane; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Plane   
    [Fish,~]=regexp(name,'fish(\d+)_','tokens','match');Fish=str2num(Fish{1}{1}); %%%to get the number of the fish 
    idx_Fish(Numbers(i):Numbers(i+1))=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter %%%to get rid of vairables we will not use anymore


Test=[];i=1;
name=strcat(MatFiles(i).name);
Rs=load(name, 'ROIs');%%% to load the ROIs
Rs=Rs.ROIs;%%% and put them on a variable
F=load(name, 'idx_components');%%% to load the ROIs idx_components again
F=F.idx_components+1;%%% and put them on a variable
Rs=Rs(:,F);%%% to combine the ROis and the idx_components
MatFiles(i).ROI=Rs;%%% to make an ROI field in MatFiles
Test(i)=size(Rs,2)==MatFiles(i).GoodNumber; %%% to check of the size of the columns of the ROIs is the same as the first (cause i=1) goodnumber 
for i = 2:length(MatFiles)%%% to take the planes (except the first one as we already did it... ) and do the same as above
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    MatFiles(i).ROI=Rs;
    Test(i)=size(Rs,2)==(MatFiles(i).GoodNumber-MatFiles(i-1).GoodNumber); %%%to to check of the size of the columns of the ROIs is the same as the goodnumber in turn. when they are you get a 1 in a column in Test until you get 25 columns. but i dont understand why we do it...
end
clearvars GC C S F N name i;

idxKmeans_final=zeros(size(ZS,1),1);%%%to make a variable with 0s the size of ZS
idxKmeans_final(idx_coef_rsq)=idxKmeans1_coef_rsq;%%% and put the filtered responses by r2 and filters

temp=[];
counter=1;
for i=goodmaps %%%to take the clusters we want
    idx_temp=find(idxKmeans_final==i);   %%%and put the filtered values in a new variable in different columns
    temp{counter}=idx_temp;    
    counter=counter+1;    
end

Numbers(1)=0; 

colors = distinguishable_colors(length(goodmaps),[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors

colors = colors*256;



%%%%%this is to do it together with the raster plots and histograms of
%%%%%localization of the clusters in fish and levels
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;counter2=1;xplot=floor(sqrt(length(goodmaps)));yplot=ceil(length(goodmaps)/xplot);
for i=goodmaps
    idx_temp=find(idxKmeans1_coef_rsq==i);
    subplot(rows,4,counter);plot(mean(ZS(idx_coef_rsq(idx_temp),:),1),'color',colors(counter2,:)/256); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS(idx_coef_rsq(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane(idx_coef_rsq(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish(idx_coef_rsq(idx_temp)));%%% for the fish location
    %counter=counter+1;
    counter2=counter2+1
     counter=counter+4;
end


clear Calcium Fitness GoodCalcium

 save('f60_r2050_CL4.mat','-v7.3');
 
 