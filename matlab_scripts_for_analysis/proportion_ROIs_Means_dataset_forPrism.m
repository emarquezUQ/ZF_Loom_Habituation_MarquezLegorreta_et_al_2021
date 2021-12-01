
%%%%% this script is to sort out the proportion analysis to be able to plot
%%%%% the 4 datasets and individual data points instead of bar graphs (for
%%%%% figure 3 b and c). 

%%

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);


RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Hindbrain','Cerebellum'};


%%

%%% merging things together with the proportions according to the brain
%%% areas and the rsq

All_Prop_CL4_per_fish_rsq=load('All_Prop_CL4_per_fish_rsq_2.mat');

All_Prop_CL4_per_fish_rsq=All_Prop_CL4_per_fish_rsq.All_Prop_CL4_per_fish;


All_Prop_CL4_per_fish_rsq_perCluster=load('All_Prop_CL4_per_fish_rsq_perCluster2.mat');

All_Prop_CL4_per_fish_rsq_perCluster=All_Prop_CL4_per_fish_rsq_perCluster.All_Prop_CL4_per_fish;


All_Prop_CL4_per_fish_nonloom_perRegion=load('All_Prop_CL4_per_fish_brainRegAllROIs_2.mat');

All_Prop_CL4_per_fish_nonloom_perRegion=All_Prop_CL4_per_fish_nonloom_perRegion.All_Prop_CL4_per_fish;

%%
%%%%% checking

colors={'r'; 'm'; 'g'; 'b'};
counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);
for data=1:size(datasets,1)
plot(mean(All_Prop_CL4_per_fish_rsq.(datasets(data,:)).(clustersF{clust})),strcat('.',colors{data}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
hold on;
end
counter=counter+1;
end


colors={'r'; 'm'; 'g'; 'b'};
counter=1; 
figure;
for clust=1:size(clustersF,1)
subplot(3,3,counter);
for data=1:size(datasets,1)
plot(mean(All_Prop_CL4_per_fish_rsq_perCluster.(datasets(data,:)).(clustersF{clust})),strcat('.',colors{data}));set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
hold on;
end
counter=counter+1;
end

%% Making the matrices to copy to Prism

clusters={'fasthab';'slopehab';'nonhab'};


%%%% for proportion of cluster ID of loom responsive ROIs within brain regions 
%%% for figure 3c
All_prop_brainReg_rsq=[];
counter=1; 

for clust=1:size(clusters,1)

for data=1:size(datasets,1)
temp_mean=nanmean(All_Prop_CL4_per_fish_rsq.(datasets(data,:)).(clusters{clust}));

All_prop_brainReg_rsq(:,counter)=temp_mean;
counter=counter+1;
end
end


%%%% for proportion of cluster ID of loom responsive ROIs in all brain 
%%%% for figure 3b
All_prop_brainReg_rsq_allbrain=[];
counter=1; 

for clust=1:size(clusters,1)

for data=1:size(datasets,1)
temp_mean=nanmean(All_Prop_CL4_per_fish_rsq_perCluster.(datasets(data,:)).(clusters{clust}));

All_prop_brainReg_rsq_allbrain(:,counter)=temp_mean;
counter=counter+1;
end
end


%%%% for proportion of cluster ID of loom responsive ROIs withinbrain regions over all ROIs of that region 
%%%% for possible suppl figure
All_prop_brainReg=[];
counter=1; 

for clust=1:size(clusters,1)

for data=1:size(datasets,1)
temp_mean=nanmean(All_Prop_CL4_per_fish_nonloom_perRegion.(datasets(data,:)).(clusters{clust}));

All_prop_brainReg(:,counter)=temp_mean;
counter=counter+1;
end
end



save('All_Prop_CL4_Means_perDataset_forPrism.mat','All_prop_brainReg_rsq','All_prop_brainReg_rsq_allbrain','All_prop_brainReg');


%% 
%%% to see box plots too

% colors={'r'; 'm'; 'g'; 'b'};
% counter=1; 
% figure;
% for clust=1:size(clustersF,1)
% subplot(3,3,counter);
% for data=1:size(datasets,1)
% boxplot(All_Prop_CL4_per_fish_rsq.(datasets(data,:)).(clustersF{clust}),'Colors',colors{data});set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
% hold on;
% end
% %boxplot(All_Prop_CL4_per_fish.s20.(clustersF{clust}),'Colors','r');set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
% counter=counter+1;
% end
% 
% 
% 
% colors={'r'; 'm'; 'g'; 'b'};
% counter=1; 
% figure;
% for clust=1:size(clustersF,1)
% subplot(3,3,counter);
% for data=1:size(datasets,1)
% boxplot(All_Prop_CL4_per_fish_rsq_perCluster.(datasets(data,:)).(clustersF{clust}),'Colors',colors{data});set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
% hold on;
% end
% %boxplot(All_Prop_CL4_per_fish.s20.(clustersF{clust}),'Colors','r');set(gca,'xticklabel',RegionList),xtickangle(45),title((clustersF{clust,1})); ylim([0 1])
% counter=counter+1;
% end
