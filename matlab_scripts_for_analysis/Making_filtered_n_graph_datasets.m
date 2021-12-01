

%%%% This script is to get the filtered datasets to be able to perform
%%%% analysis more easily. 

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};


load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);


%%% getting the nodes data 
load('Nodes_N_means_alldatasets2.mat')

Nodes.Nod_coor=Nodes2.Mod_loc;
Nodes.Nod_clustID=Nodes2.Mod_clust;
Nodes.Nod_brainID=Nodes2.Mod_brain;

Nodes.f20.NodeMats=Nodes2.f20.mean_matrix_K;
Nodes.f60.NodeMats=Nodes2.f60.mean_matrix_K;
Nodes.s20.NodeMats=Nodes2.s20.mean_matrix_K;
Nodes.s60.NodeMats=Nodes2.s60.mean_matrix_K;




%%% getting the coords, ZS, ClustID and FishID of the rsq ROIs 

%load('ZS_N_Fish_all.mat','idx_Fish_all','idx_f60_adjust','idx_s20_adjust','idx_s60_adjust');
load('ZS_N_Fish_all.mat')

f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');

%%% cluster names
clusters=fieldnames(Nodes2.f20.ROIs_BrainNK_idx.fish_17);


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


%%% I need to order them correctly. based on ClustersF.
 %%% fasthabs first and no sound
goodorder_clust=[4 2 5 6 1 7];


%% getting the ROIs 
%%%% f60
unique(idx_Fish_all(idx_f60_adjust)) %%%%

unique(idx_Fish_all(idx_f60_adjust(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned))) 


ZS_rsq_f60=ZS_CN((idx_f60_adjust(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned)),:);

FishID_rsq_f60=idx_Fish_all((idx_f60_adjust(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned)));

ROIsCoord_rsq_f60=ROI_temp2_all((idx_f60_adjust(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned)),:);

ClustID_rsq_f60=zeros(size(idx_Fish_all));
for clust=goodorder_clust
   

temp_idx_clust=clust_f60_CL7_cleaned_cell.(clustersF{clust});

    
        temp_idx=idx_f60_adjust(temp_idx_clust);
    figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.   
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROI_temp2_all(temp_idx,1),ROI_temp2_all(temp_idx,2),'filled');
    view(-90,90);
      
    ClustID_rsq_f60(temp_idx)=clust;
                          
 end
  
    
ClustID_rsq_f60=ClustID_rsq_f60(idx_f60_adjust(f60_cleaned_idxs.idx_rsq_test_f60short3_cleaned));


%%% to check if it worked

for clust=goodorder_clust
     
        temp_idx=find(ClustID_rsq_f60==clust);
    figure;plot(mean(ZS_rsq_f60(temp_idx,:))); %% to check if it works.   
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROIsCoord_rsq_f60(temp_idx,1),ROIsCoord_rsq_f60(temp_idx,2),'filled');
    view(-90,90);
                             
end
 
%%% 
F60.ZS_rsq_f60=ZS_rsq_f60;
F60.FishID_rsq_f60=FishID_rsq_f60;
F60.ROIsCoord_rsq_f60=ROIsCoord_rsq_f60;
F60.ClustID_rsq_f60=ClustID_rsq_f60;

%% now for s20



ZS_rsq_s20=ZS_CN((idx_s20_adjust(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned)),:);

FishID_rsq_s20=idx_Fish_all((idx_s20_adjust(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned)));

ROIsCoord_rsq_s20=ROI_temp2_all((idx_s20_adjust(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned)),:);

ClustID_rsq_s20=zeros(size(idx_Fish_all));
for clust=goodorder_clust
   
temp_idx_clust=clust_s20_CL7_cleaned_cell.(clustersF{clust});

        temp_idx=idx_s20_adjust(temp_idx_clust);
    figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.   
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROI_temp2_all(temp_idx,1),ROI_temp2_all(temp_idx,2),'filled');
    view(-90,90);
      
    ClustID_rsq_s20(temp_idx)=clust;
                          
 end
      
ClustID_rsq_s20=ClustID_rsq_s20(idx_s20_adjust(s20_cleaned_idxs.idx_rsq_test_s20short_cleaned));


%%% to check if it worked

for clust=goodorder_clust
     
        temp_idx=find(ClustID_rsq_s20==clust);
    figure;plot(mean(ZS_rsq_s20(temp_idx,:))); %% to check if it works.   
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROIsCoord_rsq_s20(temp_idx,1),ROIsCoord_rsq_s20(temp_idx,2),'filled');
    view(-90,90);
                             
end
 
%%% 
S20.ZS_rsq_s20=ZS_rsq_s20;
S20.FishID_rsq_s20=FishID_rsq_s20;
S20.ROIsCoord_rsq_s20=ROIsCoord_rsq_s20;
S20.ClustID_rsq_s20=ClustID_rsq_s20;


%% now for f20



ZS_rsq_f20=ZS_CN((f20_cleaned_idxs.idx_rsq_test_f20short_cleaned),:);

FishID_rsq_f20=idx_Fish_all((f20_cleaned_idxs.idx_rsq_test_f20short_cleaned));

ROIsCoord_rsq_f20=ROI_temp2_all((f20_cleaned_idxs.idx_rsq_test_f20short_cleaned),:);

ClustID_rsq_f20=zeros(size(idx_Fish_all));
for clust=goodorder_clust
   
temp_idx_clust=clust_f20_CL7_cleaned_cell.(clustersF{clust});

        temp_idx=temp_idx_clust;
    figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.   
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROI_temp2_all(temp_idx,1),ROI_temp2_all(temp_idx,2),'filled');
    view(-90,90);
      
    ClustID_rsq_f20(temp_idx)=clust;
                          
 end
      
ClustID_rsq_f20=ClustID_rsq_f20(f20_cleaned_idxs.idx_rsq_test_f20short_cleaned);


%%% to check if it worked

for clust=goodorder_clust
     
        temp_idx=find(ClustID_rsq_f20==clust);
    figure;plot(mean(ZS_rsq_f20(temp_idx,:))); %% to check if it works.   
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROIsCoord_rsq_f20(temp_idx,1),ROIsCoord_rsq_f20(temp_idx,2),'filled');
    view(-90,90);
                             
end
 
%%% 
F20.ZS_rsq_f20=ZS_rsq_f20;
F20.FishID_rsq_f20=FishID_rsq_f20;
F20.ROIsCoord_rsq_f20=ROIsCoord_rsq_f20;
F20.ClustID_rsq_f20=ClustID_rsq_f20;


%% now for s60



ZS_rsq_s60=ZS_CN((idx_s60_adjust(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned)),:);

FishID_rsq_s60=idx_Fish_all((idx_s60_adjust(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned)));

ROIsCoord_rsq_s60=ROI_temp2_all((idx_s60_adjust(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned)),:);

ClustID_rsq_s60=zeros(size(idx_Fish_all));
for clust=goodorder_clust
   
temp_idx_clust=clust_s60_CL7_cleaned_cell.(clustersF{clust});

        temp_idx=idx_s60_adjust(temp_idx_clust);
    figure;plot(mean(ZS_CN(temp_idx,:))); %% to check if it works.   
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROI_temp2_all(temp_idx,1),ROI_temp2_all(temp_idx,2),'filled');
    view(-90,90);
      
    ClustID_rsq_s60(temp_idx)=clust;
                          
 end
      
ClustID_rsq_s60=ClustID_rsq_s60(idx_s60_adjust(s60_cleaned_idxs.idx_rsq_test_s60short_cleaned));


%%% to check if it worked

for clust=goodorder_clust
     
        temp_idx=find(ClustID_rsq_s60==clust);
    figure;plot(mean(ZS_rsq_s60(temp_idx,:))); %% to check if it works.   
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROIsCoord_rsq_s60(temp_idx,1),ROIsCoord_rsq_s60(temp_idx,2),'filled');
    view(-90,90);
                             
end
 
%%% 
S60.ZS_rsq_s60=ZS_rsq_s60;
S60.FishID_rsq_s60=FishID_rsq_s60;
S60.ROIsCoord_rsq_s60=ROIsCoord_rsq_s60;
S60.ClustID_rsq_s60=ClustID_rsq_s60;


%%

for clust=goodorder_clust
     
        temp_idx1=find(ClustID_rsq_s20==clust);
        temp_idx2=find(ClustID_rsq_f60==clust);
        
    figure;plot(mean(ZS_rsq_s20(temp_idx1,:))); %% to check if it works.   
     hold on;plot(mean(ZS_rsq_f60(temp_idx2,:)));
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROIsCoord_rsq_s20(temp_idx1,1),ROIsCoord_rsq_s20(temp_idx1,2),'filled');
    hold on;
    scatter(ROIsCoord_rsq_f60(temp_idx2,1),ROIsCoord_rsq_f60(temp_idx2,2),'filled');
    view(-90,90);
                             
end


for clust=goodorder_clust
      
    temp_idx1=find(S20.ClustID_rsq_s20==clust);
    temp_idx2=find(F60.ClustID_rsq_f60==clust);
    
    figure;plot(mean(S20.ZS_rsq_s20(temp_idx1,:)));   
     hold on;plot(mean(F60.ZS_rsq_f60(temp_idx2,:)));
    figure;
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(S20.ROIsCoord_rsq_s20(temp_idx1,1),S20.ROIsCoord_rsq_s20(temp_idx1,2),'filled');
    hold on;
    scatter(F60.ROIsCoord_rsq_f60(temp_idx2,1),F60.ROIsCoord_rsq_f60(temp_idx2,2),'filled');
    view(-90,90);
                             
end

save('FnS_all_samples.mat','F20','F60','S20','S60','Zbrain_brainMask2D');

save('Nodes_FnS_all_samples.mat','Nodes','Zbrain_brainMask2D','RegionList');


%%%% for the matrices

%%%% is all in there. but i need the "keep" list. 
load('graph_loomNdataset2.mat','Data_corrMat2');

datasets=['f20'; 'f60'; 's20'; 's60'];
%%%%% to filter nodes as I did for the paper: 


%%% to see which ROIs have less representation
counter=1;
figure;
meanProp=[];
NaN_nodes={};
for data=1:4
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

fish_perNode=[];
meanProp_good=[];
for  data=1:4
    fish_perNode(:,data)=abs(sum(NaN_nodes{data},2)-size(NaN_nodes{data},2));
    meanProp_good(:,data)=1-(sum(NaN_nodes{data},2)/size(NaN_nodes{data},2));
    
end

discard=find(min(fish_perNode,[],2)<3); %%% 
%discard=find(min(meanProp_good,[],2)<0.25); 
keep=find(ismember([1:102],discard)==0);


figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_brain(keep));view(-90,90);
title('Model Nodes');
%%% adding numbers
figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');hold on;gscatter(Nodes2.Mod_loc(:,1),Nodes2.Mod_loc(:,2),Nodes2.Mod_brain(:));view(-90,90);
title('Model Nodes');
a = [1:102]'; b = num2str(a); c = cellstr(b);
dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
text(x+dx, y+dy, c);


counter=1;
figure;
for data=1:4
     datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1}(keep,keep)); caxis([-1 1])%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep));caxis([-1 1]) %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3}(keep,keep));caxis([-1 1])
subplot(4,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4}(keep,keep));caxis([-1 1])
subplot(4,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep,keep));caxis([-1 1])
subplot(4,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6}(keep,keep));caxis([-1 1])
subplot(4,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep));caxis([-1 1])%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep));caxis([-1 1]) %% for 11th loom
title(datatemp);

counter=counter+8;
end

save('graphs_FnS_all.mat','Data_corrMat2','keep','datasets');


%%% now the matrices of the fmr1 and WT. note that the keep list will be different 

load('NodesNgraphFmr1Loomhab4.mat')

groupnames=fieldnames(Nodes4.ROIs_idx);

counter=1;
figure;
for g=1:3
     group=groupnames{g,1};
subplot(3,8,counter);imagesc(Data_corrMat4.(group).Mean_corrMat{1,1}(keep,keep)); %% for pre loom
subplot(3,8,counter+1);imagesc(Data_corrMat4.(group).Mean_corrMat{1,2}(keep,keep)); %% for 1st loom
subplot(3,8,counter+2);imagesc(Data_corrMat4.(group).Mean_corrMat{1,3}(keep,keep));
subplot(3,8,counter+3);imagesc(Data_corrMat4.(group).Mean_corrMat{1,4}(keep,keep));
subplot(3,8,counter+4);imagesc(Data_corrMat4.(group).Mean_corrMat{1,5}(keep,keep));
subplot(3,8,counter+5);imagesc(Data_corrMat4.(group).Mean_corrMat{1,6}(keep,keep));
subplot(3,8,counter+6);imagesc(Data_corrMat4.(group).Mean_corrMat{1,11}(keep,keep));%% for 10th loom
subplot(3,8,counter+7);imagesc(Data_corrMat4.(group).Mean_corrMat{1,12}(keep,keep)); %% for 11th loom
title(group);

counter=counter+8;
end


save('graphs_fmr1Exp.mat','Data_corrMat4','keep','groupnames');

