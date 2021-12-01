
%%%% this script is to get the proportion of ROIs for the four datasets with
%%%% the CL4 and per brain region

%%% for CL4

%%% this script is to look at the number/proportion of ROIs across fish and
%%% datasets and to be able to look also at in across brain regions. 

%%%% first you need to load what you need


load('All_More_BrainReg.mat')


datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);

load('Zbrain_Masks.mat');


RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Hindbrain','Cerebellum'};



%%
%%% for f20


load('final_F20_step1.mat','idx_Fish_f20');
f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
%BrainReg_F20=load('BrainReg_F20.mat');

PerBrainRegions_f20=PerBrainRegions.f20;

idx_rsq_test_f20short_cleaned=f20_cleaned_idxs.idx_rsq_test_f20short_cleaned;
idx_f20_multisense_cleaned=f20_cleaned_idxs.idx_multisense_cleaned;
idx_rsq_Mov_cleaned=f20_cleaned_idxs.idx_rsq_Mov_cleaned;
clust_f20_CL4_cleaned=f20_cleaned_idxs.clust_f20_CL4_cleaned;

clust_f20_CL4_cleaned_cell={};
clust=fieldnames(clust_f20_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL4_cleaned_cell.(clustersF{j,1})=clust_f20_CL4_cleaned.(clust{j});   
end    



%%% to look at the proportion per brain region to see if it works
%%% NOTE: i am basing the proportion on the response rsq filtered ROIs, not in
%%% the total number of ROIs in that brian region. it could be intersting
%%% to look at both. 

counter=1;
figure;
for j=1:size(clustersF,1)
    
    idx_temp=clust_f20_CL4_cleaned_cell.(clustersF{j,1});
    
    for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_f20.(regionName).idx,idx_temp))/sum(ismember(PerBrainRegions_f20.(regionName).idx,idx_rsq_test_f20short_cleaned));
end
subplot(3,3,counter);bar(ratioperRegion);set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{j,1})); ylim([0 1])
    
   %hold on; 
   counter=counter+1;
end

%%% now i need to do it an average per fish

Prop_CL4_per_fish_f20={};
fish=unique(idx_Fish_f20);
for clust=1:size(clustersF,1)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_f20_CL4_cleaned_cell.(clustersF{clust,1}));
for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_f20.(regionName).idx,temp_idx))/sum(ismember(PerBrainRegions_f20.(regionName).idx,intersect(idx_rsq_test_f20short_cleaned,temp_idx_fish)));
end

Prop_CL4_per_fish_f20{clust}(tempfish,:)=ratioperRegion;
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 
counter=1;
figure;
for clust=1:size(clustersF,1)

for i=1:length(fish)
   
   subplot(3,3,counter); plot(Prop_CL4_per_fish_f20{clust}(i,:),'.');set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
     
    hold on;
    %subplot(3,3,counter);bar(mean(Prop_CL4_per_fish_f20{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])

end
 
counter=counter+1;
end


counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);bar(mean(Prop_CL4_per_fish_f20{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;
end

counter=1; 
figure;
for clust=1:size(clustersF,1)


subplot(3,3,counter);boxplot(Prop_CL4_per_fish_f20{clust});set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;

end

%%
All_Prop_CL4_per_fish={};

for j=1:size(clustersF,1)
    All_Prop_CL4_per_fish.f20.(clustersF{j,1})=Prop_CL4_per_fish_f20{j};  
end

clearvars -except All_Prop_CL4_per_fish PerBrainRegions clustersF clustersS datasets RegionList

%%% the box plot is very nice.
counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);boxplot(All_Prop_CL4_per_fish.f20.(clustersF{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;
end

%%
%%% for f60

load('final_F60_step1_2.mat','idx_Fish_f60');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
%BrainReg_F60=load('BrainReg_F60.mat');

PerBrainRegions_f60=PerBrainRegions.f60;

idx_rsq_test_f60short_cleaned=f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned;
idx_f60_multisense_cleaned=f60_cleaned_idxs.idx_multisense_cleaned;
idx_rsq_Mov_cleaned=f60_cleaned_idxs.idx_rsq_Mov_cleaned;
clust_f60_CL4_cleaned=f60_cleaned_idxs.clust_f60_CL4_cleaned;

clust_f60_CL4_cleaned_cell={};
clust=fieldnames(clust_f60_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f60_CL4_cleaned_cell.(clustersF{j,1})=clust_f60_CL4_cleaned.(clust{j});   
end    



%%% to look at the proportion per brain region to see if it works
%%% NOTE: i am basing the proportion on the response rsq filtered ROIs, not in
%%% the total number of ROIs in that brian region. 

counter=1;
figure;
for j=1:size(clustersF,1)
    
    idx_temp=clust_f60_CL4_cleaned_cell.(clustersF{j,1});
    
    for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_f60.(regionName).idx,idx_temp))/sum(ismember(PerBrainRegions_f60.(regionName).idx,idx_rsq_test_f60short_cleaned));
end
subplot(3,3,counter);bar(ratioperRegion);set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{j,1})); ylim([0 1])
    
   %hold on; 
   counter=counter+1;
end

%%% now i need to do it an average per fish

Prop_CL4_per_fish_f60={};
fish=unique(idx_Fish_f60);
fish(find(fish==47))=[]; %%% cause I also took out fish 47...
for clust=1:size(clustersF,1)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f60==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_f60_CL4_cleaned_cell.(clustersF{clust,1}));
for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_f60.(regionName).idx,temp_idx))/sum(ismember(PerBrainRegions_f60.(regionName).idx,intersect(idx_rsq_test_f60short_cleaned,temp_idx_fish)));
end

Prop_CL4_per_fish_f60{clust}(tempfish,:)=ratioperRegion;
end
end

%%

for j=1:size(clustersF,1)
    All_Prop_CL4_per_fish.f60.(clustersF{j,1})=Prop_CL4_per_fish_f60{j};  
end

clearvars -except All_Prop_CL4_per_fish PerBrainRegions clustersF clustersS datasets RegionList

%%% the box plot is very nice.
counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);boxplot(All_Prop_CL4_per_fish.f60.(clustersF{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
 counter=counter+1;
end







%%
%%% for s20

load('final_S20_step1.mat','idx_Fish_s20');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
%BrainReg_S20=load('BrainReg_S20.mat');

PerBrainRegions_s20=PerBrainRegions.s20;

idx_rsq_test_s20short_cleaned=s20_cleaned_idxs.idx_rsq_test_s20short_cleaned;
idx_s20_multisense_cleaned=s20_cleaned_idxs.idx_multisense_cleaned;
idx_rsq_Mov_cleaned=s20_cleaned_idxs.idx_rsq_Mov_cleaned;
clust_s20_CL4_cleaned=s20_cleaned_idxs.clust_s20_CL4_cleaned;

clust_s20_CL4_cleaned_cell={};
clust=fieldnames(clust_s20_CL4_cleaned);
for j=1:size(clustersS,1)
 clust_s20_CL4_cleaned_cell.(clustersS{j,1})=clust_s20_CL4_cleaned.(clust{j});   
end    



%%% to look at the proportion per brain region to see if it works
%%% NOTE: i am basing the proportion on the response rsq filtered ROIs, not in
%%% the total number of ROIs in that brian region. 

counter=1;
figure;
for j=1:size(clustersS,1)
    
    idx_temp=clust_s20_CL4_cleaned_cell.(clustersS{j,1});
    
    for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_s20.(regionName).idx,idx_temp))/sum(ismember(PerBrainRegions_s20.(regionName).idx,idx_rsq_test_s20short_cleaned));
end
subplot(3,3,counter);bar(ratioperRegion);set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersS{j,1})); ylim([0 1])
    
   %hold on; 
   counter=counter+1;
end

%%% now i need to do it an average per fish

Prop_CL4_per_fish_s20={};
fish=unique(idx_Fish_s20);
for clust=1:size(clustersS,1)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s20==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_s20_CL4_cleaned_cell.(clustersS{clust,1}));
for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_s20.(regionName).idx,temp_idx))/sum(ismember(PerBrainRegions_s20.(regionName).idx,intersect(idx_rsq_test_s20short_cleaned,temp_idx_fish)));
end

Prop_CL4_per_fish_s20{clust}(tempfish,:)=ratioperRegion;
end
end

%%

for j=1:size(clustersS,1)
    All_Prop_CL4_per_fish.s20.(clustersS{j,1})=Prop_CL4_per_fish_s20{j};  
end

clearvars -except All_Prop_CL4_per_fish PerBrainRegions clustersS clustersF datasets RegionList

%%% the box plot is very nice.
counter=1; 
figure;
for clust=1:size(clustersS,1)
subplot(3,3,counter);boxplot(All_Prop_CL4_per_fish.s20.(clustersS{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersS{clust,1})); ylim([0 1])
 counter=counter+1;
end



%%


%%
%%% for s60

load('final_S60_step1.mat','idx_Fish_s60');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');
%BrainReg_S60=load('BrainReg_S60.mat');

PerBrainRegions_s60=PerBrainRegions.s60;

idx_rsq_test_s60short_cleaned=s60_cleaned_idxs.idx_rsq_test_s60short_cleaned;
idx_s60_multisense_cleaned=s60_cleaned_idxs.idx_multisense_cleaned;
idx_rsq_Mov_cleaned=s60_cleaned_idxs.idx_rsq_Mov_cleaned;
clust_s60_CL4_cleaned=s60_cleaned_idxs.clust_s60_CL4_cleaned;

clust_s60_CL4_cleaned_cell={};
clust=fieldnames(clust_s60_CL4_cleaned);
for j=1:size(clustersS,1)
 clust_s60_CL4_cleaned_cell.(clustersS{j,1})=clust_s60_CL4_cleaned.(clust{j});   
end    



%%% to look at the proportion per brain region to see if it works
%%% NOTE: i am basing the proportion on the response rsq filtered ROIs, not in
%%% the total number of ROIs in that brian region. 

counter=1;
figure;
for j=1:size(clustersS,1)
    
    idx_temp=clust_s60_CL4_cleaned_cell.(clustersS{j,1});
    
    for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_s60.(regionName).idx,idx_temp))/sum(ismember(PerBrainRegions_s60.(regionName).idx,idx_rsq_test_s60short_cleaned));
end
subplot(3,3,counter);bar(ratioperRegion);set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersS{j,1})); ylim([0 1])
    
   %hold on; 
   counter=counter+1;
end

%%% now i need to do it an average per fish

Prop_CL4_per_fish_s60={};
fish=unique(idx_Fish_s60);
for clust=1:size(clustersS,1)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s60==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_s60_CL4_cleaned_cell.(clustersS{clust,1}));
for i=1:length(RegionList)
    regionName=RegionList{i};
    ratioperRegion(:,i)=sum(ismember(PerBrainRegions_s60.(regionName).idx,temp_idx))/sum(ismember(PerBrainRegions_s60.(regionName).idx,intersect(idx_rsq_test_s60short_cleaned,temp_idx_fish)));
end

Prop_CL4_per_fish_s60{clust}(tempfish,:)=ratioperRegion;
end
end

%%

for j=1:size(clustersS,1)
    All_Prop_CL4_per_fish.s60.(clustersS{j,1})=Prop_CL4_per_fish_s60{j};  
end

clearvars -except All_Prop_CL4_per_fish PerBrainRegions clustersS clustersF datasets RegionList

%%% the box plot is very nice.
counter=1; 
figure;
for clust=1:size(clustersS,1)
subplot(3,3,counter);boxplot(All_Prop_CL4_per_fish.s60.(clustersS{clust}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersS{clust,1})); ylim([0 1])
 counter=counter+1;
end


save('All_Prop_CL4_per_fish_rsq_2.mat','All_Prop_CL4_per_fish');

%%
%%%% to plot them all together
colors={'r'; 'm'; 'g'; 'b'};
counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);
for data=1:size(datasets,1)
boxplot(All_Prop_CL4_per_fish.(datasets(data,:)).(clustersF{clust}),'Colors',colors{data});set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
hold on;
end
%boxplot(All_Prop_CL4_per_fish.s20.(clustersF{clust}),'Colors','r');set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
counter=counter+1;
end

%%

