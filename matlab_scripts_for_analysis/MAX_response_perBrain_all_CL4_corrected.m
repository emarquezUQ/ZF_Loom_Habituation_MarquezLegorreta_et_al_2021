
%%%% this script is to find the maximum responses per fish, per cluster and
%%%% per brain region. Doing it for the 4 datasets (f20, f60, s20, s60). 

%%% for CL4


load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx', 'S_trim');

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);

load('Zbrain_Masks.mat');
load('All_More_BrainReg.mat','PerBrainRegions');


RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Cerebellum','Tegmentum','Habenula',};

%%
%%%for f20

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');


%%

%%% this is to get the means of each cluster per fish and per brain region of interest. 



clust_f20_CL4_cleaned=f20_cleaned_idxs.clust_f20_CL4_cleaned;

clust_f20_CL4_cleaned_cell={};
clust=fieldnames(clust_f20_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL4_cleaned_cell.(clustersF{j,1})=clust_f20_CL4_cleaned.(clust{j});   
end    

%test_idx_CL4_per_fishNbrain_f20=struct;
Means_CL4_per_fishNbrain_f20=struct;
fish=unique(idx_Fish_f20);


for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f20_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.f20.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_mean=NaN(1,1344); 
else
temp_mean=mean(ZS_f20(temp_idx,:));
end
Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(tempfish,:)=temp_mean;


end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_f20{1,clust}))

end
end


%%% now ploting the means.


for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f20_CL4_cleaned_cell))
subplot(3,3,counter);
if 2<sum(~isnan(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(:,1)))
plot(nanmean(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})));
else
end
counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f20_perfishNbrain2=struct;

for brain=1:length(RegionList)

for clust=1:length(clustersF)

if 2<sum(~isnan(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(:,1)))
    
    for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1));

        for k=1:30

        temp_Max_resp(1,k+1)= max(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_f20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint zero
        end


    Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2);%%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
    end
else 
end    
    

end
end


%%% to check it worked

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
temp_regionClust=fieldnames(Max_resp_f20_perfishNbrain2.(RegionList{brain}));
for clust=1:length(temp_regionClust)
subplot(3,3,counter);

 
for i=1:length(fish)
   plot(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(temp_regionClust{clust})(i,:)); 
    hold on;
end


counter=counter+1;
end
end


%%% now ploting the means.
%%% 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
temp_regionClust=fieldnames(Max_resp_f20_perfishNbrain2.(RegionList{brain}));
for clust=1:length(temp_regionClust)
subplot(3,3,counter);
plot(nanmean(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(temp_regionClust{clust})));


counter=counter+1;


end
end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain2




%%
%%%for f60

load('final_F60_step1_2.mat','ZS_f60','idx_Fish_f60','ZS_short_F60');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');


%%

%%% this is to get the means of each cluster per fish and per brain region of interest. 



clust_f60_CL4_cleaned=f60_cleaned_idxs.clust_f60_CL4_cleaned;

clust_f60_CL4_cleaned_cell={};
clust=fieldnames(clust_f60_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f60_CL4_cleaned_cell.(clustersF{j,1})=clust_f60_CL4_cleaned.(clust{j});   
end    


Means_CL4_per_fishNbrain_f60=struct;
fish=unique(idx_Fish_f60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47


for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f60==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f60_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.f60.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_mean=NaN(1,1344); 
else
temp_mean=mean(ZS_f60(temp_idx,ZS_short_F60));
end
Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(tempfish,:)=temp_mean;


end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_f20{1,clust}))

end
end


%%% now ploting the means.

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_f60_CL4_cleaned_cell))
subplot(3,3,counter);
if 2<sum(~isnan(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(:,1)))
plot(nanmean(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})));
else
end
counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f60_perfishNbrain2=struct;

for brain=1:length(RegionList)

for clust=1:length(clustersF)

if 2<sum(~isnan(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(:,1)))
    
    for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1));

        for k=1:30

        temp_Max_resp(1,k+1)= max(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_f60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint zero
        end


    Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2);%%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
    end
else 
end    
    

end
end


%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
temp_regionClust=fieldnames(Max_resp_f60_perfishNbrain2.(RegionList{brain}));
for clust=1:length(temp_regionClust)
subplot(3,3,counter);

 
for i=1:length(fish)
   plot(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(temp_regionClust{clust})(i,:)); 
    hold on;
end


counter=counter+1;
end
end


%%% now ploting the means.
%%% 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
temp_regionClust=fieldnames(Max_resp_f60_perfishNbrain2.(RegionList{brain}));
for clust=1:length(temp_regionClust)
subplot(3,3,counter);
plot(nanmean(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(temp_regionClust{clust})));


counter=counter+1;


end
end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain2 Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain2



%%

%%% now for s20

load('final_S20_step1.mat','ZS_s20','idx_Fish_s20');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');


%%

%%% this is to get the means of each cluster per fish and per brain region of interest. 



clust_s20_CL4_cleaned=s20_cleaned_idxs.clust_s20_CL4_cleaned;

clust_s20_CL4_cleaned_cell={};
clust=fieldnames(clust_s20_CL4_cleaned);
for j=1:size(clustersS,1)
 clust_s20_CL4_cleaned_cell.(clustersS{j,1})=clust_s20_CL4_cleaned.(clust{j});   
end    


Means_CL4_per_fishNbrain_s20=struct;
fish=unique(idx_Fish_s20);



for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s20_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.s20.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_mean=NaN(1,1344); 
else
temp_mean=mean(ZS_s20(temp_idx,S_trim));
end
Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(tempfish,:)=temp_mean;


end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_f20{1,clust}))

end
end


%%% now ploting the means.

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_s20_CL4_cleaned_cell))
subplot(3,3,counter);
if 2<sum(~isnan(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(:,1)))
plot(nanmean(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})));
else
end
counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_s20_perfishNbrain2=struct;

for brain=1:length(RegionList)

for clust=1:length(clustersF)

if 2<sum(~isnan(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(:,1)))
    
    for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1));

        for k=1:30

        temp_Max_resp(1,k+1)= max(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_s20.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint zero
        end


    Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2);%%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
    end
else 
end    
    

end
end


%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
temp_regionClust=fieldnames(Max_resp_s20_perfishNbrain2.(RegionList{brain}));
for clust=1:length(temp_regionClust)
subplot(3,3,counter);

 
for i=1:length(fish)
   plot(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(temp_regionClust{clust})(i,:)); 
    hold on;
end


counter=counter+1;
end
end


%%% now ploting the means.
%%% 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
temp_regionClust=fieldnames(Max_resp_s20_perfishNbrain2.(RegionList{brain}));
for clust=1:length(temp_regionClust)
subplot(3,3,counter);
plot(nanmean(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(temp_regionClust{clust})));


counter=counter+1;


end
end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain2 Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain2 Means_CL4_per_fishNbrain_s20 Max_resp_s20_perfishNbrain2

%%


%%% now for s60

load('final_S60_step1.mat','ZS_s60','idx_Fish_s60','ZS_short_S60');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');



%%

%%% this is to get the means of each cluster per fish and per brain region of interest. 



clust_s60_CL4_cleaned=s60_cleaned_idxs.clust_s60_CL4_cleaned;

clust_s60_CL4_cleaned_cell={};
clust=fieldnames(clust_s60_CL4_cleaned);
for j=1:size(clustersS,1)
 clust_s60_CL4_cleaned_cell.(clustersS{j,1})=clust_s60_CL4_cleaned.(clust{j});   
end    

Means_CL4_per_fishNbrain_s60=struct;
fish=unique(idx_Fish_s60);


for brain=1:length(RegionList)

for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s60==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s60_CL4_cleaned_cell.(clustersF{clust}));
temp_idx=intersect(temp_idx_clust,PerBrainRegions.s60.(RegionList{brain}).idx);

if length(temp_idx)<5
    
    temp_mean=NaN(1,1344); 
else
temp_mean=mean(ZS_s60(temp_idx,ZS_short_S60(S_trim)));
end
Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(tempfish,:)=temp_mean;


end
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 
for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
subplot(3,3,counter);
for i=1:length(fish)
   plot(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(i,:)); 
    hold on;
    
end
counter=counter+1;
%figure;plot(mean(Means_CL4_per_fishNbrain_f20{1,clust}))

end
end


%%% now ploting the means.

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
for clust=1:length(fieldnames(clust_s60_CL4_cleaned_cell))
subplot(3,3,counter);
if 2<sum(~isnan(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(:,1)))
plot(nanmean(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})));
else
end
counter=counter+1;


end
end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_s60_perfishNbrain2=struct;

for brain=1:length(RegionList)

for clust=1:length(clustersF)

if 2<sum(~isnan(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(:,1)))
    
    for  i=1:length(fish)  
    temp_Max_resp(1,1)= max(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1))-min(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(1)+1));

        for k=1:30

        temp_Max_resp(1,k+1)= max(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(i,loom_moments{1,k}))-min(Means_CL4_per_fishNbrain_s60.(RegionList{brain}).(clustersF{clust})(i,Loomf20_onset_idx(k)+1)); %%% it is k+1 cause i am adding a timepoint zero
        end


    Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,2);%%% i changed to temp_Max_resp(1,2) cause in the first one I will have my timpoint 0
    end
else 
end    
    

end
end


%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
temp_regionClust=fieldnames(Max_resp_s60_perfishNbrain2.(RegionList{brain}));
for clust=1:length(temp_regionClust)
subplot(3,3,counter);

 
for i=1:length(fish)
   plot(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(temp_regionClust{clust})(i,:)); 
    hold on;
end


counter=counter+1;
end
end


%%% now ploting the means.
%%% 

for brain=1:length(RegionList)

figure('Position',[100 0 900 900]);
counter=1;
temp_regionClust=fieldnames(Max_resp_s60_perfishNbrain2.(RegionList{brain}));
for clust=1:length(temp_regionClust)
subplot(3,3,counter);
plot(nanmean(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(temp_regionClust{clust})));


counter=counter+1;


end
end

clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Zbrain_Masks S_trim RegionList PerBrainRegions datasets clustersF clustersS Means_CL4_per_fishNbrain_f20 Max_resp_f20_perfishNbrain2 Means_CL4_per_fishNbrain_f60 Max_resp_f60_perfishNbrain2 Means_CL4_per_fishNbrain_s20 Max_resp_s20_perfishNbrain2 Means_CL4_per_fishNbrain_s60 Max_resp_s60_perfishNbrain2


%%



save('Max_response_perBrain_all_CL4_corrected.mat');


