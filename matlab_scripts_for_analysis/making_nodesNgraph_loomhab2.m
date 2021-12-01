


%%%% this script is to generate the representative nodes of each main cluster
%%%% from all datasets to be used to extract means of each fish in a way
%%%% that i can compare them.  

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

%load('All_More_BrainReg.mat','PerBrainRegions');
 load('All_More_BrainReg2.mat')

load('Nodes_N_means_alldatasets.mat','ROI_temp2_all','Zbrain_brainMask2D');

load('ZS_N_Fish_all.mat','idx_Fish_all','idx_f60_adjust','idx_s20_adjust','idx_s60_adjust');

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};


f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');




%%

%%% to order the clusters per name
clust_f20_CL7_cleaned=f20_cleaned_idxs.clust_f20_CL7_cleaned;

clust_f20_CL7_cleaned_cell={};
clust=fieldnames(clust_f20_CL7_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL7_cleaned_cell.(clustersF{j,1})=clust_f20_CL7_cleaned.(clust{j});   
end   


clust_f60_CL7_cleaned=f60_cleaned_idxs.clust_f60_CL7_cleaned;

clust_f60_CL7_cleaned_cell={};
clust=fieldnames(clust_f60_CL7_cleaned);
for j=1:size(clustersF,1)
 clust_f60_CL7_cleaned_cell.(clustersF{j,1})=clust_f60_CL7_cleaned.(clust{j});   
end   

clust_s20_CL7_cleaned=s20_cleaned_idxs.clust_s20_CL7_cleaned;

clust_s20_CL7_cleaned_cell={};
clust=fieldnames(clust_s20_CL7_cleaned);
for j=1:size(clustersS,1)
 clust_s20_CL7_cleaned_cell.(clustersS{j,1})=clust_s20_CL7_cleaned.(clust{j});   
end   


clust_s60_CL7_cleaned=s60_cleaned_idxs.clust_s60_CL7_cleaned;

clust_s60_CL7_cleaned_cell={};
clust=fieldnames(clust_s60_CL7_cleaned);
for j=1:size(clustersS,1)
 clust_s60_CL7_cleaned_cell.(clustersS{j,1})=clust_s60_CL7_cleaned.(clust{j});   
end   

%%

 %%% I need to order them. based on ClustersF.
 %%% fasthabs first and no sound
 goodorder_clust=[4 2 5 6 1 7];
 

Nodes2=struct;
Nodes2.Mod_loc=[];
Nodes2.Mod_clust=[];
Nodes2.Mod_brain=[];
Nodes2.Mod_KmeansID={};
 counter=1;
for clust=goodorder_clust
          
 idx_temp=vertcat(clust_f20_CL7_cleaned_cell.(clustersF{clust}),(idx_f60_adjust(clust_f60_CL7_cleaned_cell.(clustersF{clust})))',(idx_s20_adjust(clust_s20_CL7_cleaned_cell.(clustersF{clust})))',(idx_s60_adjust(clust_s60_CL7_cleaned_cell.(clustersF{clust})))');   
  
 
 
 for brain=1:length(RegionList)
     
     idx_brain_temp=vertcat(PerBrainRegions.f20.(RegionList{brain}).idx,(idx_f60_adjust(PerBrainRegions.f60.(RegionList{brain}).idx))',(idx_s20_adjust(PerBrainRegions.s20.(RegionList{brain}).idx))',(idx_s60_adjust(PerBrainRegions.s60.(RegionList{brain}).idx))');
     
     brain_clust_idx=intersect(idx_temp,idx_brain_temp);
     
     if length(brain_clust_idx)<200
       continue 
     elseif length(brain_clust_idx)<500 & length(brain_clust_idx)>200
         moduleN=1;
     elseif length(brain_clust_idx)<1000 & length(brain_clust_idx)>500
     moduleN=2;
     elseif length(brain_clust_idx)<3000 & length(brain_clust_idx)>1000
     moduleN=3;
      elseif length(brain_clust_idx)>3000 
     moduleN=4;
     end   
 
options = statset('UseParallel',1); [idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2_all(brain_clust_idx,:),moduleN,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

Nodes2.Mod_loc=vertcat(Nodes2.Mod_loc,Cmap_ROIs);
Nodes2.Mod_clust=vertcat(Nodes2.Mod_clust,ones(size(Cmap_ROIs,1),1)*clust);
Nodes2.Mod_brain=vertcat(Nodes2.Mod_brain,ones(size(Cmap_ROIs,1),1)*brain);


KID=unique(idxKmeans_ROIs);
for ID=1:length(KID)
 Nodes2.Mod_KmeansID{counter,1}=brain_clust_idx(find(idxKmeans_ROIs==ID));   
counter=counter+1;
end

 end

end

%%% to check if it works
clust=6;
idx_temp=vertcat(clust_f20_CL7_cleaned_cell.(clustersF{clust}),(idx_f60_adjust(clust_f60_CL7_cleaned_cell.(clustersF{clust})))',(idx_s20_adjust(clust_s20_CL7_cleaned_cell.(clustersF{clust})))',(idx_s60_adjust(clust_s60_CL7_cleaned_cell.(clustersF{clust})))');    
figure;scatter(ROI_temp2_all(brain_clust_idx,1),ROI_temp2_all(brain_clust_idx,2));

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;scatter(ROI_temp2_all(temp_KiD{4,1},1),ROI_temp2_all(temp_KiD{4,1},2),15,'filled');



figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust,'rgggbm');
view(-90,90);

%%% it seems that it worked. now to get ROIs of each fish per
%%% dataset.

%% for f20
%load('final_F20_step1.mat','idx_Fish_f20');
%%% OR (for all datasets at once)
load('ZS_N_Fish_all.mat')


%% getting the ROIs of each node per fish

fish=unique(idx_Fish_f20);
for tempfish=1:length(fish)
for clust=goodorder_clust
      
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f20_CL7_cleaned_cell.(clustersF{clust}));
 
    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes2.Mod_clust==clust);
    for node=(find(Nodes2.Mod_clust==clust))'
    
        temp_idx=intersect(temp_idx_clust,Nodes2.Mod_KmeansID{node,1});
    %figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.   
%     figure;
%     plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%     hold on;
%     scatter(ROI_temp2_all(temp_idx,1),ROI_temp2_all(temp_idx,2),'filled');
%     view(-90,90);
      
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx; 
                          
    end
  
    Nodes2.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust})=idx_ROIs_Node;
         
end
end



%%% to get the means

for tempfish=1:length(fish)
for clust=goodorder_clust
        

    for node=(find(Nodes2.Mod_clust==clust))'
    
   
    temp_idx=Nodes2.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        %plot(temp_mean);
    end
    
    Nodes2.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(102,1344);
    for clust=goodorder_clust
        
        for node=(find(Nodes2.Mod_clust==clust))'
            
        temp_mean=Nodes2.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes2.f20.mean_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes2.f20.corr_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
        
        Nodes2.f20.NaNtest_K.(strcat('fish_',num2str(fish(tempfish))))=isnan(R_temp);
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.f20.NaNtest_K.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);
max(max(Matrix_mean));

%% now for f60


%load('final_F60_step1_2.mat','idx_Fish_f60');


%% getting the ROIs of each node per fish
%Nodes_f60=struct;
fish=unique(idx_Fish_all(idx_f60_adjust));

fish(find(fish==47))=[]; %%%  fish 47 was discarded

for tempfish=1:length(fish)
for clust=goodorder_clust
   
   
temp_idx_fish=find(idx_Fish_all(idx_f60_adjust)==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f60_CL7_cleaned_cell.(clustersF{clust}));


    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes2.Mod_clust==clust);
    for node=(find(Nodes2.Mod_clust==clust))'
    
        temp_idx=intersect(idx_f60_adjust(temp_idx_clust),Nodes2.Mod_KmeansID{node,1});
%     figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.   
%     figure;
%     plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%     hold on;
%     scatter(ROI_temp2_all(temp_idx,1),ROI_temp2_all(temp_idx,2),'filled');
%     view(-90,90);
      
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx; 
                          
    end
  
    Nodes2.f60.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust})=idx_ROIs_Node;
         
    
end
end


%%% to get the means
for tempfish=1:length(fish)
for clust=goodorder_clust
        

    for node=(find(Nodes2.Mod_clust==clust))'
    
   
    temp_idx=Nodes2.f60.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        plot(temp_mean);
    end
    
    Nodes2.f60.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(102,1344);
    for clust=goodorder_clust
        
        for node=(find(Nodes2.Mod_clust==clust))'
            
        temp_mean=Nodes2.f60.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes2.f60.mean_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes2.f60.corr_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
        
        Nodes2.f60.NaNtest_K.(strcat('fish_',num2str(fish(tempfish))))=isnan(R_temp);
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.f60.NaNtest_K.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);
max(max(Matrix_mean));



%% now for s20




%% getting the ROIs of each node per fish
%Nodes_s20=struct;
fish=unique(idx_Fish_all(idx_s20_adjust));


for tempfish=1:length(fish)
for clust=goodorder_clust
   
   
temp_idx_fish=find(idx_Fish_all(idx_s20_adjust)==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s20_CL7_cleaned_cell.(clustersF{clust}));


    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes2.Mod_clust==clust);
    for node=(find(Nodes2.Mod_clust==clust))'
    
        temp_idx=intersect(idx_s20_adjust(temp_idx_clust),Nodes2.Mod_KmeansID{node,1});
%     figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.   
%     figure;
%     plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%     hold on;
%     scatter(ROI_temp2_all(temp_idx,1),ROI_temp2_all(temp_idx,2),'filled');
%     view(-90,90);
      
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx; 
                          
    end
  
    Nodes2.s20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust})=idx_ROIs_Node;
         
    
end
end


%%% to get the means
for tempfish=1:length(fish)
for clust=goodorder_clust
        

    for node=(find(Nodes2.Mod_clust==clust))'
    
   
    temp_idx=Nodes2.s20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        plot(temp_mean);
    end
    
    Nodes2.s20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(102,1344);
    for clust=goodorder_clust
        
        for node=(find(Nodes2.Mod_clust==clust))'
            
        temp_mean=Nodes2.s20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes2.s20.mean_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes2.s20.corr_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
        
        Nodes2.s20.NaNtest_K.(strcat('fish_',num2str(fish(tempfish))))=isnan(R_temp);
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.s20.NaNtest_K.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);
max(max(Matrix_mean));




%% now for s60



%% getting the ROIs of each node per fish
%Nodes_s60=struct;
fish=unique(idx_Fish_all(idx_s60_adjust));


for tempfish=1:length(fish)
for clust=goodorder_clust
   
   
temp_idx_fish=find(idx_Fish_all(idx_s60_adjust)==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s60_CL7_cleaned_cell.(clustersF{clust}));


    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes2.Mod_clust==clust);
    for node=(find(Nodes2.Mod_clust==clust))'
    
        temp_idx=intersect(idx_s60_adjust(temp_idx_clust),Nodes2.Mod_KmeansID{node,1});
%     figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.   
%     figure;
%     plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
%     hold on;
%     scatter(ROI_temp2_all(temp_idx,1),ROI_temp2_all(temp_idx,2),'filled');
%     view(-90,90);
      
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx; 
                          
    end
  
    Nodes2.s60.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust})=idx_ROIs_Node;
         
    
end
end


%%% to get the means
for tempfish=1:length(fish)
for clust=goodorder_clust
        

    for node=(find(Nodes2.Mod_clust==clust))'
    
   
    temp_idx=Nodes2.s60.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        plot(temp_mean);
    end
    
    Nodes2.s60.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(102,1344);
    for clust=goodorder_clust
        
        for node=(find(Nodes2.Mod_clust==clust))'
            
        temp_mean=Nodes2.s60.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes2.s60.mean_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes2.s60.corr_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
        
        Nodes2.s60.NaNtest_K.(strcat('fish_',num2str(fish(tempfish))))=isnan(R_temp);
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.s60.NaNtest_K.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);
max(max(Matrix_mean));





%% saving


save('Nodes_N_means_alldatasets2.mat','Nodes2','ROI_temp2_all','Zbrain_brainMask2D');

%% next step

%%%% this script is to make a correlation in each fish for each of the loom
%%%% resonses based on the nodes. Then to make an
%%%% average of these correlation matrices per dataset to be able to
%%%% compare them. 

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

load('Nodes_N_means_alldatasets2.mat','Nodes2','Zbrain_brainMask2D');

load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx');

%% getting the corrmatrix for each loom and each fish
Data_corrMat2=struct;

for data=1:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes2.(datatemp).mean_matrix_K);
    
    for tempfish=1:length(fish)
        temp_mean=Nodes2.(datatemp).mean_matrix_K.(fish{tempfish});
        
        temp_R={};
        
        for k=1:31
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 ||k==31
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat2.(datatemp).(fish{tempfish}).loomsR=temp_R;
        
    end
             
end


%% making means of each loom per dataset

for data=1:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Data_corrMat2.(datatemp));
    
    Mean_corrMat={};
    for k=1:31
        
        temp_Mean_corrMat=[];
        
        for tempfish=1:length(fish)
            
        temp_mat=Data_corrMat2.(datatemp).(fish{tempfish}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMat2.(datatemp).Mean_corrMat=Mean_corrMat;
    
end




%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust);view(-90,90);xlim([400 1350]);
title('Model cluster Nodes2');
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_brain);view(-90,90);xlim([400 1350]);
title('Model brain Nodes2');


counter=1;
figure;
for data=1:4
    datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1});caxis([-1 1]);colormap('parula'); %% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2});caxis([-1 1]); colormap('parula');%% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3});caxis([-1 1]);colormap('parula');
subplot(4,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4});caxis([-1 1]);colormap('parula');
subplot(4,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5});caxis([-1 1]);colormap('parula');
subplot(4,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6});caxis([-1 1]);colormap('parula');
subplot(4,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11});caxis([-1 1]);colormap('parula');%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12});caxis([-1 1]);colormap('parula'); %% for 11th loom
title(datatemp);

counter=counter+8;
end


%% now to plot some graphs. 

for data=1:4
    figure;
    count=1;
    datatemp=datasets(data,:);
    sgtitle(datatemp);
   
    for k=[2 3 4 5 6 11 12]
    
    R=Data_corrMat2.(datatemp).Mean_corrMat{1,k};
        
        
n=length(Nodes2.Mod_loc);

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

[~,~,weights] = find(tril(R,-1));

% create the graph object:
G = graph(s,t,weights,n);

% mark the lines to remove from the graph:
threshold = 0.75; %  minimum correlation to plot
line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>

% plot it:
subplot(1,7,count);
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];
% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*1;
axis off

% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes2.Mod_loc(:,1);
y = Nodes2.Mod_loc(:,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust,'rygcbm','.',20,'off');

view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold off;

count=count+1

    end
end


save('graph_loomNdataset2.mat','Data_corrMat2');


%% checking a few things

for data=1:4
     datatemp=datasets(data,:);
        
     fish=fieldnames(Nodes2.(datatemp).corr_matrix);
     
    for tempfish=1:length(fish)
           
    Nodes2.NaNtest.(datatemp).(fish{tempfish})=isnan(Nodes2.(datatemp).corr_matrix.(fish{tempfish}));
    end
end


%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
for data=1:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes2.(datatemp).NaNtest_K);
     
    
    Matrix_mean=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.(datatemp).NaNtest_K.(fish{f}));     
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,4,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled');colormap('jet'); caxis([0 1]);% colorbar;
view(-90,90);
title(datatemp);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));


end

