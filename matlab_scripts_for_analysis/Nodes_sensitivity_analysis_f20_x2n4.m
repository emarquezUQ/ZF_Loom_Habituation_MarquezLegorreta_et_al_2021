

%%%%% this is the script to test number of node sensitivity to the results. 
%%%%% to do some graph analysis with different amount of nodes than the
%%%%% 99 nodes and see if the results hold. With the f20 dataset.

%%
%%%% Loading:

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

load('All_More_BrainReg2.mat')

load('Nodes_N_means_alldatasets.mat','ROI_temp2_all','Zbrain_brainMask2D');

load('ZS_N_Fish_all.mat','idx_Fish_all','idx_f60_adjust','idx_s20_adjust','idx_s60_adjust');

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');

load('Max_response_perBrain_all_CL4_corrected_distribution.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx');



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


 %%% to order them correctly. based on ClustersF.
 %%% fasthabs first and no sound
 goodorder_clust=[4 2 5 6 1 7];
 
 
%% 
%%%%% Same method as for the 99 nodes
Nodes1=struct;
Nodes1.Mod_loc=[];
Nodes1.Mod_clust=[];
Nodes1.Mod_brain=[];
Nodes1.Mod_KmeansID={};
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

Nodes1.Mod_loc=vertcat(Nodes1.Mod_loc,Cmap_ROIs);
Nodes1.Mod_clust=vertcat(Nodes1.Mod_clust,ones(size(Cmap_ROIs,1),1)*clust);
Nodes1.Mod_brain=vertcat(Nodes1.Mod_brain,ones(size(Cmap_ROIs,1),1)*brain);

KID=unique(idxKmeans_ROIs);
for ID=1:length(KID)
 Nodes1.Mod_KmeansID{counter,1}=brain_clust_idx(find(idxKmeans_ROIs==ID));   
counter=counter+1;
end

 end

end

%%% to check if it works

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
gscatter(Nodes1.Mod_loc(:,1),Nodes1.Mod_loc(:,2),Nodes1.Mod_clust,'rgggbm');
view(-90,90);

%%

%%%%% to get the response traces for f20

%% for f20
load('ZS_N_idx_Fish_all.mat','idx_Fish_f20');
%%% and 
load('ZS_N_Fish_all.mat')


%% getting the ROIs of each node per fish

fish=unique(idx_Fish_f20);
for tempfish=1:length(fish)
for clust=goodorder_clust
      
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f20_CL7_cleaned_cell.(clustersF{clust}));
 
    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes1.Mod_clust==clust);
    for node=(find(Nodes1.Mod_clust==clust))'
    
        temp_idx=intersect(temp_idx_clust,Nodes1.Mod_KmeansID{node,1});
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx; 
                          
    end
  
    Nodes1.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust})=idx_ROIs_Node;
         
end
end



%%% to get the means

for tempfish=1:length(fish)
for clust=goodorder_clust
        

    for node=(find(Nodes1.Mod_clust==clust))'
    
   
    temp_idx=Nodes1.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        %plot(temp_mean);
    end
    
    Nodes1.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(102,1344);
    for clust=goodorder_clust
        
        for node=(find(Nodes1.Mod_clust==clust))'
            
        temp_mean=Nodes1.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes1.f20.mean_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes1.f20.corr_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
        
        Nodes1.f20.NaNtest_K.(strcat('fish_',num2str(fish(tempfish))))=isnan(R_temp);
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes1.f20.NaNtest_K.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);
max(max(Matrix_mean));

%%

%% getting the corrmatrix for each loom and each fish
Data_corrMat1=struct;

for data=1%:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes1.(datatemp).mean_matrix_K);
    
    for tempfish=1:length(fish)
        temp_mean=Nodes1.(datatemp).mean_matrix_K.(fish{tempfish});
        
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
        
    Data_corrMat1.(datatemp).(fish{tempfish}).loomsR=temp_R;
        
    end
             
end


%% making means of each loom per dataset

for data=1%:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Data_corrMat1.(datatemp));
    
    Mean_corrMat={};
    for k=1:31
        
        temp_Mean_corrMat=[];
        
        for tempfish=1:length(fish)
            
        temp_mat=Data_corrMat1.(datatemp).(fish{tempfish}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMat1.(datatemp).Mean_corrMat=Mean_corrMat;
    
end




%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes1.Mod_loc(:,1),Nodes1.Mod_loc(:,2),Nodes1.Mod_clust);view(-90,90);xlim([400 1350]);
title('Model cluster Nodes1');
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes1.Mod_loc(:,1),Nodes1.Mod_loc(:,2),Nodes1.Mod_brain);view(-90,90);xlim([400 1350]);
title('Model brain Nodes1');


counter=1;
figure;
for data=1%:4
    datatemp=datasets(data,:);
subplot(1,8,counter);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,1});caxis([-1 1]);colormap('parula'); %% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,2});caxis([-1 1]); colormap('parula');%% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,3});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+3);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,4});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+4);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,5});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+5);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,6});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+6);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,11});caxis([-1 1]);colormap('parula');%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,12});caxis([-1 1]);colormap('parula'); %% for 11th loom
title(datatemp);

counter=counter+8;
end



%% 
%%%%% checking if the nodes are well represented. 

%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
NaN_nodes={};
for data=1%:4
     datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes1.(datatemp).NaNtest_K);  
     
    
    Matrix_mean=[];
    temp_NaN_nodes=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes1.(datatemp).NaNtest_K.(fish{f}));  
 temp=double(diag(Nodes1.(datatemp).NaNtest_K.(fish{f})));
 temp_NaN_nodes=horzcat(temp_NaN_nodes,temp);
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,4,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes1.Mod_loc(:,1),Nodes1.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colormap('jet'); caxis([0 1]);%colorbar;
view(-90,90);
title(datatemp);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));

NaN_nodes{data}=temp_NaN_nodes;
end

fish_perNode1=[];
meanProp_good1=[];
for  data=1%:4
    fish_perNode1(:,data)=abs(sum(NaN_nodes{data},2)-size(NaN_nodes{data},2));
    meanProp_good1(:,data)=1-(sum(NaN_nodes{data},2)/size(NaN_nodes{data},2));
    
end

discard1=find(min(meanProp_good1,[],2)<0.25); %%%  
keep1=find(ismember([1:102],discard1)==0);


figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes1.Mod_loc(keep1,1),Nodes1.Mod_loc(keep1,2),Nodes1.Mod_brain(keep1));view(-90,90);
title('Model Nodes');
%%% adding numbers
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes1.Mod_loc(:,1),Nodes1.Mod_loc(:,2),Nodes1.Mod_brain(:));view(-90,90);
title('Model Nodes');

counter=1;
figure;
for data=1%:4
     datatemp=datasets(data,:);
subplot(1,8,counter);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,1}(keep1,keep1)); caxis([-1 1])%% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,2}(keep1,keep1));caxis([-1 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,3}(keep1,keep1));caxis([-1 1])
subplot(1,8,counter+3);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,4}(keep1,keep1));caxis([-1 1])
subplot(1,8,counter+4);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,5}(keep1,keep1));caxis([-1 1])
subplot(1,8,counter+5);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,6}(keep1,keep1));caxis([-1 1])
subplot(1,8,counter+6);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,11}(keep1,keep1));caxis([-1 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,12}(keep1,keep1));caxis([-1 1]) %% for 11th loom
title(datatemp);

counter=counter+8;
end

%%

%% testing_graph_analysis_with_BCT3

%%%% I will also compare it with the original analysis. 

OrginalMatAll_corrected=load('graph_analysis_loomhab3.mat');

%%
%%%% To look at density, the number of degrees and strength.

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected1=struct;
for data=1%:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat1.(datatemp).Mean_corrMat{1,moment(m)}(keep1,keep1)),0.75);
     MatAll_corrected1.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%

%kden = density_und(CIJ);

for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected1.(datatemp).(loom{i}).Mat);

MatAll_corrected1.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected1.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end


%%%% degrees and Strength

for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)

deg=degrees_und(MatAll_corrected1.(datatemp).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected1.(datatemp).(loom{i}).Mat));

MatAll_corrected1.(datatemp).(loom{i}).deg=deg;
MatAll_corrected1.(datatemp).(loom{i}).str=str;
end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected1.(datatemp).(loom{i}).deg;
end
plot(nanmean(temp));
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).deg;
end
plot(nanmean(temp));
hold on;

end


%%%% participation and cluster coeficient


for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)

temp_mat=MatAll_corrected1.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes1.Mod_clust(keep1));
Ccoef=clustering_coef_wu(temp_mat);

MatAll_corrected1.(datatemp).(loom{i}).Ccoef=Ccoef;
MatAll_corrected1.(datatemp).(loom{i}).P=P;

end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected1.(datatemp).(loom{i}).Ccoef;
end
plot(nanmean(temp));
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).Ccoef;
end
plot(nanmean(temp));
hold on;

end

%% now for 2 times more nodes

%%%%% this is to try to have 2 times more nodes
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
         moduleN=2;
     elseif length(brain_clust_idx)<1000 & length(brain_clust_idx)>500
     moduleN=4;
     elseif length(brain_clust_idx)<3000 & length(brain_clust_idx)>1000
     moduleN=6;
      elseif length(brain_clust_idx)>3000 
     moduleN=8;
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

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_clust,'rgggbm');
view(-90,90);


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

    temp_matrix=NaN(length(Nodes2.Mod_brain),1344);
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


%%


%% getting the corrmatrix for each loom and each fish
Data_corrMat2=struct;

for data=1%:4
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

for data=1%:4
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
title('Model cluster Nodes1');
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_brain);view(-90,90);xlim([400 1350]);
title('Model brain Nodes1');


counter=1;
figure;
for data=1%:4
    datatemp=datasets(data,:);
subplot(1,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1});caxis([-1 1]);colormap('parula'); %% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2});caxis([-1 1]); colormap('parula');%% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11});caxis([-1 1]);colormap('parula');%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12});caxis([-1 1]);colormap('parula'); %% for 11th loom
title(datatemp);

counter=counter+8;
end



%% 

%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
NaN_nodes={};
for data=1%:4
     datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes2.(datatemp).NaNtest_K);  
     
    
    Matrix_mean=[];
    temp_NaN_nodes=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes2.(datatemp).NaNtest_K.(fish{f}));  
 temp=double(diag(Nodes2.(datatemp).NaNtest_K.(fish{f})));
 temp_NaN_nodes=horzcat(temp_NaN_nodes,temp);
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,4,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colormap('jet'); caxis([0 1]);%colorbar;
view(-90,90);
title(datatemp);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));

NaN_nodes{data}=temp_NaN_nodes;
end

fish_perNode2=[];
meanProp_good2=[];
for  data=1%:4
    fish_perNode2(:,data)=abs(sum(NaN_nodes{data},2)-size(NaN_nodes{data},2));
    meanProp_good2(:,data)=1-(sum(NaN_nodes{data},2)/size(NaN_nodes{data},2));
    
end

discard2=find(min(meanProp_good2,[],2)<0.25); %%%
keep2=find(ismember([1:length(Nodes2.Mod_brain)],discard2)==0);

length(discard2)
length(keep2)

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(keep2,1),Nodes2.Mod_loc(keep2,2),Nodes2.Mod_brain(keep2));view(-90,90);
title('Model Nodes');
%%% adding numbers
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_brain(:));view(-90,90);
title('Model Nodes');

counter=1;
figure;
for data=1%:4
     datatemp=datasets(data,:);
subplot(1,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1}(keep2,keep2)); caxis([-1 1])%% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep2,keep2));caxis([-1 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3}(keep2,keep2));caxis([-1 1])
subplot(1,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4}(keep2,keep2));caxis([-1 1])
subplot(1,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep2,keep2));caxis([-1 1])
subplot(1,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6}(keep2,keep2));caxis([-1 1])
subplot(1,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep2,keep2));caxis([-1 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep2,keep2));caxis([-1 1]) %% for 11th loom
title(datatemp);

counter=counter+8;
end

%%
%% testing_graph_analysis_with_BCT3

%%
%%%% To look at density, the number of degrees and strength.

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected2=struct;
for data=1%:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep2,keep2)),0.75);
     MatAll_corrected2.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%


for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected2.(datatemp).(loom{i}).Mat);

MatAll_corrected2.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected2.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end




%%%% degrees and Strength

for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)

deg=degrees_und(MatAll_corrected2.(datatemp).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected2.(datatemp).(loom{i}).Mat));

MatAll_corrected2.(datatemp).(loom{i}).deg=deg;
MatAll_corrected2.(datatemp).(loom{i}).str=str;
end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected2.(datatemp).(loom{i}).deg;
end
plot(nanmean(temp));
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).deg;
end
plot(nanmean(temp));
hold on;

end


%%%% participation and cluster coeficient

for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)

temp_mat=MatAll_corrected2.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes2.Mod_clust(keep2));
Ccoef=clustering_coef_wu(temp_mat);

MatAll_corrected2.(datatemp).(loom{i}).Ccoef=Ccoef;
MatAll_corrected2.(datatemp).(loom{i}).P=P;

end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected2.(datatemp).(loom{i}).P;
end
plot(nanmean(temp));
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).P;
end
plot(nanmean(temp));
hold on;

end

%% now for 4 times more nodes

%%%%% this is to try to have 4 times more nodes
Nodes4=struct;
Nodes4.Mod_loc=[];
Nodes4.Mod_clust=[];
Nodes4.Mod_brain=[];
Nodes4.Mod_KmeansID={};
 counter=1;
for clust=goodorder_clust
          
 idx_temp=vertcat(clust_f20_CL7_cleaned_cell.(clustersF{clust}),(idx_f60_adjust(clust_f60_CL7_cleaned_cell.(clustersF{clust})))',(idx_s20_adjust(clust_s20_CL7_cleaned_cell.(clustersF{clust})))',(idx_s60_adjust(clust_s60_CL7_cleaned_cell.(clustersF{clust})))');   
  
 
 
 for brain=1:length(RegionList)
     
     idx_brain_temp=vertcat(PerBrainRegions.f20.(RegionList{brain}).idx,(idx_f60_adjust(PerBrainRegions.f60.(RegionList{brain}).idx))',(idx_s20_adjust(PerBrainRegions.s20.(RegionList{brain}).idx))',(idx_s60_adjust(PerBrainRegions.s60.(RegionList{brain}).idx))');
     
     brain_clust_idx=intersect(idx_temp,idx_brain_temp);
     
     if length(brain_clust_idx)<200
       continue 
     elseif length(brain_clust_idx)<500 & length(brain_clust_idx)>200
         moduleN=4;
     elseif length(brain_clust_idx)<1000 & length(brain_clust_idx)>500
     moduleN=8;
     elseif length(brain_clust_idx)<3000 & length(brain_clust_idx)>1000
     moduleN=12;
      elseif length(brain_clust_idx)>3000 
     moduleN=16;
     end   
 
options = statset('UseParallel',1); [idxKmeans_ROIs Cmap_ROIs]=kmeans(ROI_temp2_all(brain_clust_idx,:),moduleN,'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

Nodes4.Mod_loc=vertcat(Nodes4.Mod_loc,Cmap_ROIs);
Nodes4.Mod_clust=vertcat(Nodes4.Mod_clust,ones(size(Cmap_ROIs,1),1)*clust);
Nodes4.Mod_brain=vertcat(Nodes4.Mod_brain,ones(size(Cmap_ROIs,1),1)*brain);


KID=unique(idxKmeans_ROIs);
for ID=1:length(KID)
 Nodes4.Mod_KmeansID{counter,1}=brain_clust_idx(find(idxKmeans_ROIs==ID));   
counter=counter+1;
end

 end

end

%%% to check if it works

figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
gscatter(Nodes4.Mod_loc(:,1),Nodes4.Mod_loc(:,2),Nodes4.Mod_clust,'rgggbm');
view(-90,90);


%% getting the ROIs of each node per fish

fish=unique(idx_Fish_f20);
for tempfish=1:length(fish)
for clust=goodorder_clust
      
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f20_CL7_cleaned_cell.(clustersF{clust}));
 
    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes4.Mod_clust==clust);
    for node=(find(Nodes4.Mod_clust==clust))'
    
        temp_idx=intersect(temp_idx_clust,Nodes4.Mod_KmeansID{node,1});

      
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx; 
                          
    end
  
    Nodes4.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust})=idx_ROIs_Node;
         
end
end



%%% to get the means

for tempfish=1:length(fish)
for clust=goodorder_clust
        

    for node=(find(Nodes4.Mod_clust==clust))'
    
   
    temp_idx=Nodes4.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(temp_idx,:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(temp_idx,:));
        %plot(temp_mean);
    end
    
    Nodes4.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(length(Nodes4.Mod_brain),1344);
    for clust=goodorder_clust
        
        for node=(find(Nodes4.Mod_clust==clust))'
            
        temp_mean=Nodes4.f20.ROIs_BrainNK_idx.(strcat('fish_',num2str(fish(tempfish)))).(clustersF{clust}).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes4.f20.mean_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes4.f20.corr_matrix_K.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
        
        Nodes4.f20.NaNtest_K.(strcat('fish_',num2str(fish(tempfish))))=isnan(R_temp);
end


%%


%% getting the corrmatrix for each loom and each fish
Data_corrMat4=struct;

for data=1%:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes4.(datatemp).mean_matrix_K);
    
    for tempfish=1:length(fish)
        temp_mean=Nodes4.(datatemp).mean_matrix_K.(fish{tempfish});
        
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
        
    Data_corrMat4.(datatemp).(fish{tempfish}).loomsR=temp_R;
        
    end
             
end


%% making means of each loom per dataset

for data=1%:4
    datatemp=datasets(data,:);
    
    fish=fieldnames(Data_corrMat4.(datatemp));
    
    Mean_corrMat={};
    for k=1:31
        
        temp_Mean_corrMat=[];
        
        for tempfish=1:length(fish)
            
        temp_mat=Data_corrMat4.(datatemp).(fish{tempfish}).loomsR{1,k};
        
        temp_Mean_corrMat=cat(3,temp_Mean_corrMat,temp_mat);
        end
        
        Mean_corrMat{1,k}=nanmean(temp_Mean_corrMat,3);
    end
    Data_corrMat4.(datatemp).Mean_corrMat=Mean_corrMat;
    
end




%%
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes4.Mod_loc(:,1),Nodes4.Mod_loc(:,2),Nodes4.Mod_clust);view(-90,90);xlim([400 1350]);
title('Model cluster Nodes1');
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes4.Mod_loc(:,1),Nodes4.Mod_loc(:,2),Nodes4.Mod_brain);view(-90,90);xlim([400 1350]);
title('Model brain Nodes1');


counter=1;
figure;
for data=1%:4
    datatemp=datasets(data,:);
subplot(1,8,counter);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,1});caxis([-1 1]);colormap('parula'); %% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,2});caxis([-1 1]); colormap('parula');%% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,3});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+3);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,4});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+4);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,5});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+5);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,6});caxis([-1 1]);colormap('parula');
subplot(1,8,counter+6);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,11});caxis([-1 1]);colormap('parula');%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,12});caxis([-1 1]);colormap('parula'); %% for 11th loom
title(datatemp);

counter=counter+8;
end



%% 

%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
NaN_nodes={};
for data=1%:4
     datatemp=datasets(data,:);
    
    fish=fieldnames(Nodes4.(datatemp).NaNtest_K);  
     
    
    Matrix_mean=[];
    temp_NaN_nodes=[];
for f=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes4.(datatemp).NaNtest_K.(fish{f}));  
 temp=double(diag(Nodes4.(datatemp).NaNtest_K.(fish{f})));
 temp_NaN_nodes=horzcat(temp_NaN_nodes,temp);
end

%Matrix_mean=(sum(Matrix_mean,3))/length(fish); %% is the same
Matrix_mean=nanmean(Matrix_mean,3);  


subplot(1,4,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(Nodes4.Mod_loc(:,1),Nodes4.Mod_loc(:,2),20,(1-mean(Matrix_mean)),'filled'); colormap('jet'); caxis([0 1]);%colorbar;
view(-90,90);
title(datatemp);
hold off;

counter=counter+1;

meanProp=horzcat(meanProp,(mean(Matrix_mean)'));

NaN_nodes{data}=temp_NaN_nodes;
end

fish_perNode4=[];
meanProp_good4=[];
for  data=1%:4
    fish_perNode4(:,data)=abs(sum(NaN_nodes{data},2)-size(NaN_nodes{data},2));
    meanProp_good4(:,data)=1-(sum(NaN_nodes{data},2)/size(NaN_nodes{data},2));
    
end

discard4=find(min(meanProp_good4,[],2)<0.25); %%%  
keep4=find(ismember([1:length(Nodes4.Mod_brain)],discard4)==0);

length(discard4)
length(keep4)

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes4.Mod_loc(keep4,1),Nodes4.Mod_loc(keep4,2),Nodes4.Mod_brain(keep4));view(-90,90);
title('Model Nodes');
%%% adding numbers
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes4.Mod_loc(:,1),Nodes4.Mod_loc(:,2),Nodes4.Mod_brain(:));view(-90,90);
title('Model Nodes');



counter=1;
figure;
for data=1%:4
     datatemp=datasets(data,:);
subplot(1,8,counter);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,1}(keep4,keep4)); caxis([-1 1])%% for pre loom
subplot(1,8,counter+1);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,2}(keep4,keep4));caxis([-1 1]) %% for 1st loom
subplot(1,8,counter+2);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,3}(keep4,keep4));caxis([-1 1])
subplot(1,8,counter+3);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,4}(keep4,keep4));caxis([-1 1])
subplot(1,8,counter+4);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,5}(keep4,keep4));caxis([-1 1])
subplot(1,8,counter+5);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,6}(keep4,keep4));caxis([-1 1])
subplot(1,8,counter+6);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,11}(keep4,keep4));caxis([-1 1])%% for 10th loom
subplot(1,8,counter+7);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,12}(keep4,keep4));caxis([-1 1]) %% for 11th loom
title(datatemp);

counter=counter+8;
end

%%
%% testing_graph_analysis_with_BCT3

%%
%%%% To look at density, the number of degrees and strength.

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected4=struct;
for data=1%:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat4.(datatemp).Mean_corrMat{1,moment(m)}(keep4,keep4)),0.75);
     MatAll_corrected4.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%

for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected4.(datatemp).(loom{i}).Mat);

MatAll_corrected4.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected4.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end




%%%% degrees and Strength

for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)

deg=degrees_und(MatAll_corrected4.(datatemp).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected4.(datatemp).(loom{i}).Mat));

MatAll_corrected4.(datatemp).(loom{i}).deg=deg;
MatAll_corrected4.(datatemp).(loom{i}).str=str;
end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected4.(datatemp).(loom{i}).deg;
end
plot(nanmean(temp));
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).deg;
end
plot(nanmean(temp));
hold on;

end


%%%% participation and cluster coeficient

for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)

temp_mat=MatAll_corrected4.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes4.Mod_clust(keep4));
Ccoef=clustering_coef_wu(temp_mat);

MatAll_corrected4.(datatemp).(loom{i}).Ccoef=Ccoef;
MatAll_corrected4.(datatemp).(loom{i}).P=P;

end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected4.(datatemp).(loom{i}).P;
end
plot(nanmean(temp));
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).P;
end
plot(nanmean(temp));
hold on;

end

%%

save('nodes_sensitivity_analysis_f20.mat','Nodes1','Nodes2','Nodes4','Data_corrMat1','Data_corrMat2','Data_corrMat4','keep1','keep2','keep4','discard1','discard2','discard4','MatAll_corrected1','MatAll_corrected2','MatAll_corrected4','OrginalMatAll_corrected','Zbrain_brainMask2D');

%%



%%%% making some figures to compare

OriginalData_corrMat=load('graphs_FnS_all.mat');

%%%%% matrices

counter=1;
figure;
for data=1%:4
     datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,1}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,2}(OriginalData_corrMat.keep,OriginalData_corrMat.keep));caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,3}(OriginalData_corrMat.keep,OriginalData_corrMat.keep));caxis([-1 1])
subplot(4,8,counter+3);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,4}(OriginalData_corrMat.keep,OriginalData_corrMat.keep));caxis([-1 1])
subplot(4,8,counter+4);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,5}(OriginalData_corrMat.keep,OriginalData_corrMat.keep));caxis([-1 1])
subplot(4,8,counter+5);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,6}(OriginalData_corrMat.keep,OriginalData_corrMat.keep));caxis([-1 1])
subplot(4,8,counter+6);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,11}(OriginalData_corrMat.keep,OriginalData_corrMat.keep));caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,12}(OriginalData_corrMat.keep,OriginalData_corrMat.keep));caxis([-1 1]) %% for 11th loom
title(datatemp);

counter=counter+8;
end


%counter=1;
for data=1%:4
     datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,1}(keep1,keep1)); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,2}(keep1,keep1));caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,3}(keep1,keep1));caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,4}(keep1,keep1));caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,5}(keep1,keep1));caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,6}(keep1,keep1));caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,11}(keep1,keep1));caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat1.(datatemp).Mean_corrMat{1,12}(keep1,keep1));caxis([-1 1]) %% for 11th loom
title(strcat(datatemp,'/',num2str(length(keep1))));

counter=counter+8;
end


%counter=1;
for data=1%:4
     datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1}(keep2,keep2)); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep2,keep2));caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3}(keep2,keep2));caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4}(keep2,keep2));caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep2,keep2));caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6}(keep2,keep2));caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep2,keep2));caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep2,keep2));caxis([-1 1]) %% for 11th loom
title(strcat(datatemp,'/',num2str(length(keep2))));

counter=counter+8;
end


%counter=1;
for data=1%:4
     datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,1}(keep4,keep4)); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,2}(keep4,keep4));caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,3}(keep4,keep4));caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,4}(keep4,keep4));caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,5}(keep4,keep4));caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,6}(keep4,keep4));caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,11}(keep4,keep4));caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,12}(keep4,keep4));caxis([-1 1]) %% for 11th loom
title(strcat(datatemp,'/',num2str(length(keep4))));

counter=counter+8;
end


%%%% density, mean participation and mean cluster coef


figure;
subplot(1,3,1);
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected1.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end


temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected2.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end


temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected4.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end
title('density');




subplot(1,3,2);
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).P;
end
plot(nanmean(temp));
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected1.(datatemp).(loom{i}).P;
end
plot(nanmean(temp));
hold on;

end


temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected2.(datatemp).(loom{i}).P;
end
plot(nanmean(temp));
hold on;

end


temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected4.(datatemp).(loom{i}).P;
end
plot(nanmean(temp));
hold on;

end
title('Participation');


subplot(1,3,3);
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).Ccoef;
end
plot(nanmean(temp));
hold on;

end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected1.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected1.(datatemp).(loom{i}).Ccoef;
end
plot(nanmean(temp));
hold on;

end


temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected2.(datatemp).(loom{i}).Ccoef;
end
plot(nanmean(temp));
hold on;

end


temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(:,i)=MatAll_corrected4.(datatemp).(loom{i}).Ccoef;
end
plot(nanmean(temp));
hold on;

end
title('clustering coeficient');

%% cross-validation and comparing variability


%%%% for original data
for data=1%:4

names = fieldnames(OriginalData_corrMat.Data_corrMat2.(datasets(data,:)));
CorrMatrices_mean_original=zeros(31,length(names)-1,99,99);
for loom=1:31        
    for fish_rem_nb=1:length(names)-1
        temp=nan(length(names)-2,99,99);    
        counter=1;
        for fish_nb=1:length(names)-1
            if fish_nb ~= fish_rem_nb
                fish_name=names(fish_nb);
                temp(counter,:,:)=OriginalData_corrMat.Data_corrMat2.(datasets(data,:)).(fish_name{1}).loomsR{1,loom}(OriginalData_corrMat.keep,OriginalData_corrMat.keep);
                counter=counter+1;            
            end
        end
        CorrMatrices_mean_original(loom,fish_rem_nb,:,:)=squeeze(nanmean(temp,1));
    end       
end

end


%%%% for ~200 nodes
for data=1%:4

names = fieldnames(Data_corrMat2.(datasets(data,:)));
Nnodes=length(keep2);
CorrMatrices_mean_2=zeros(31,length(names)-1,Nnodes,Nnodes);
for loom=1:31        
    for fish_rem_nb=1:length(names)-1
        temp=nan(length(names)-2,Nnodes,Nnodes);    
        counter=1;
        for fish_nb=1:length(names)-1
            if fish_nb ~= fish_rem_nb
                fish_name=names(fish_nb);
                temp(counter,:,:)=Data_corrMat2.(datasets(data,:)).(fish_name{1}).loomsR{1,loom}(keep2,keep2);
                counter=counter+1;            
            end
        end
        CorrMatrices_mean_2(loom,fish_rem_nb,:,:)=squeeze(nanmean(temp,1));
    end       
end

end



%%%% for ~400 nodes
for data=1%:4

names = fieldnames(Data_corrMat4.(datasets(data,:)));
Nnodes=length(keep4);
CorrMatrices_mean_4=zeros(31,length(names)-1,Nnodes,Nnodes);
for loom=1:31        
    for fish_rem_nb=1:length(names)-1
        temp=nan(length(names)-2,Nnodes,Nnodes);    
        counter=1;
        for fish_nb=1:length(names)-1
            if fish_nb ~= fish_rem_nb
                fish_name=names(fish_nb);
                temp(counter,:,:)=Data_corrMat4.(datasets(data,:)).(fish_name{1}).loomsR{1,loom}(keep4,keep4);
                counter=counter+1;            
            end
        end
        CorrMatrices_mean_4(loom,fish_rem_nb,:,:)=squeeze(nanmean(temp,1));
    end       
end

end


%% to show it

%%%% 5 random group of matrices without a fish at random

figure;
counter=1;
for fish=randperm(length(names)-1,5)
for loom=[1 2 3 4 5 6 11 12]

    subplot(5,8,counter);imagesc(squeeze(squeeze(CorrMatrices_mean_2(loom,fish,:,:))));caxis([-1 1])
    
    counter=counter+1;
end
end



%%%% changes in density. plot the density for f20 and then as scatter plot
%%%% the individual densities of the substracted means

%%%%% getting the thresholed matrices (I need the BCT toolbox)


%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)),0.75);
     MatAll_corrected.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%


for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected.(datatemp).(loom{i}).Mat);

MatAll_corrected.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end



%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_crossval_original=[];

for fish=1:length(names)-1

    moment=[1 2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    %loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(squeeze(squeeze(CorrMatrices_mean_original(moment(m),fish,:,:)))),0.75);
     MatAll_corrected_crossval_original(m,fish,:,:)=Mat;
     end

end

%%%% calculating the density


crossval_density=[];
for fish=1:length(names)-1
for m=1:length(moment)

temp_kden=density_und(squeeze(squeeze(MatAll_corrected_crossval_original(m,fish,:,:))));

crossval_density(m,fish)=temp_kden;
end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

for fish=1:length(names)-1
for m=1:length(moment)-1    
scatter(m,crossval_density(m+1,fish),'filled');
hold on;
end
end
hold off;

%% now for ~200



%%%% changes in density. plot the density for f20 and then as scatter plot
%%%% the individual densities of the substracted means

%%%%% getting the thresholed matrices (I need the BCT toolbox)



%% density
%

for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected2.(datatemp).(loom{i}).Mat);

MatAll_corrected2.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected2.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end



%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_crossval_2=[];

for fish=1:length(names)-1

    moment=[1 2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    %loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(squeeze(squeeze(CorrMatrices_mean_2(moment(m),fish,:,:)))),0.75);
     MatAll_corrected_crossval_2(m,fish,:,:)=Mat;
     end

end

%%%% calculating the density


crossval_density=[];
for fish=1:length(names)-1
for m=1:length(moment)

temp_kden=density_und(squeeze(squeeze(MatAll_corrected_crossval_2(m,fish,:,:))));

crossval_density(m,fish)=temp_kden;
end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected2.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

for fish=1:length(names)-1
for m=1:length(moment)-1    
scatter(m,crossval_density(m+1,fish),'filled');
hold on;
end
end
hold off;


%% now for ~400


for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected4.(datatemp).(loom{i}).Mat);

MatAll_corrected4.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected4.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end



%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_crossval_4=[];

for fish=1:length(names)-1

    moment=[1 2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    %loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(squeeze(squeeze(CorrMatrices_mean_2(moment(m),fish,:,:)))),0.75);
     MatAll_corrected_crossval_4(m,fish,:,:)=Mat;
     end

end

%%%% calculating the density


crossval_density=[];
for fish=1:length(names)-1
for m=1:length(moment)

temp_kden=density_und(squeeze(squeeze(MatAll_corrected_crossval_4(m,fish,:,:))));

crossval_density(m,fish)=temp_kden;
end
end


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected4.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

for fish=1:length(names)-1
for m=1:length(moment)-1    
scatter(m,crossval_density(m+1,fish),'filled');
hold on;
end
end
hold off;


