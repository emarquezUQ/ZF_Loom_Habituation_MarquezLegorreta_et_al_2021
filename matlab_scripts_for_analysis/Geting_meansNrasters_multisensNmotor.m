

%%% this script is to get the rasters of the sound, multisensory and
%%% movements raster plots. 
%%% based on the "Geting_meansNrasters.m" script, so i can keep the same scale. 

%%%  doing it for f20

load('f20_cleaned_idxs.mat');

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');

load('All_More_BrainReg2.mat');

load('Tail_mov_F20.mat','idx_rsq_Mov');

load('means_F20_CL4n7.mat')

%%%% for the rasterplot of sound responsive ROIs
figure;imagesc(ZS_f20(clust_f20_CL7_cleaned.clust_f20_CL7_3_cleaned,:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_sound.tif');


%%%% for the rasterplot of multisensory (visual and auditory) responsive ROIs
figure;imagesc(ZS_f20(idx_multisense_cleaned,:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_multisense.tif');


%%% Note: after cleaning the movement ROIs
%%% the idxs got ordered in a different way. it's needed to correct this
%%% so it can be seen by fish. 

idx_rsq_Mov_cleaned3=zeros(size(idx_rsq_Mov));
idx_rsq_Mov_cleaned3=find(ismember(idx_rsq_Mov,idx_rsq_Mov_cleaned));
idx_rsq_Mov_cleaned3=idx_rsq_Mov(idx_rsq_Mov_cleaned3);


figure;scatter(ROI_temp2.f20(idx_rsq_Mov_cleaned,1),ROI_temp2.f20(idx_rsq_Mov_cleaned,2),'filled');
hold on;
scatter(ROI_temp2.f20(idx_rsq_Mov_cleaned3,1),ROI_temp2.f20(idx_rsq_Mov_cleaned3,2),'filled');
hold on;


figure;imagesc(ZS_f20(idx_rsq_Mov_cleaned3,:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_mov_corrected.tif');

%%% now the mean of one of the fish

unique(idx_Fish_f20)

%%% with fish 44

%figure;plot(mean(ZS_f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==44)),:)));

%%% to get the values to put in prism

mean(ZS_f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==44)),:));

%%% to see distribution
figure;
scatter(ROI_temp2.f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==44)),1),ROI_temp2.f20(intersect(idx_rsq_Mov_cleaned3,find(idx_Fish_f20==44)),2),'filled');



