
%%%% this script is to generate the representative nodes of each main cluster
%%%% from all datasets to be used to extract means of each fish in a way
%%%% that i can compare them. 

%%%% this code was the first atempt and has some variables saved that we
%%%% used later but we ended up using a slighly different approach to
%%%% generate the nodes. See "making_nodesNgraph_loomhab2.m"


datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

%load('All_More_BrainReg.mat','PerBrainRegions');
 load('All_More_BrainReg2.mat')

load('ZS_N_Fish_all.mat','idx_Fish_all','idx_f60_adjust','idx_s20_adjust','idx_s60_adjust');

%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain'};

%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain','Cerebellum','Tegmentum','Habenula',};


f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');


%% making the shape of the brain
load('Zbrain_Masks.mat');

%%% this is to make a matrix with the location of the main structures 
Zbrain_brainMask=vertcat(Zbrain_Masks{[76 113 259 274 294],3}); %%% 78 is the eyes
 
Zbrain_brainMask=unique(Zbrain_brainMask,'rows');

figure;scatter(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),'.');

Zbrain_brainMask2D_bound=boundary(Zbrain_brainMask(:,1),Zbrain_brainMask(:,2),1);

figure;plot(Zbrain_brainMask(Zbrain_brainMask2D_bound,1),Zbrain_brainMask(Zbrain_brainMask2D_bound,2));

Zbrain_brainMask2D=[];
Zbrain_brainMask2D(:,1)=Zbrain_brainMask(Zbrain_brainMask2D_bound,1);
Zbrain_brainMask2D(:,2)=Zbrain_brainMask(Zbrain_brainMask2D_bound,2);

figure;plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');

%%

%%% to get all the ROIs together
ROI_temp2_all=[];
for data=1:4
 ROI_temp2_all=vertcat(ROI_temp2_all,ROI_temp2.(datasets(data,:)));   
     
end

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


%%%
moduleN=[13 22 2 26 10 26 3]; %%% for 100 modules. without sound
%moduleN=[13 22 2 26 10 26 1]; %%% for 100 modules
%moduleN=[6 10 2 13 5 13 1];  %%% for 50 modules
%moduleN=[4 5 2 6 3 5 1];  %%% for 26 modules
Nodes=struct;
Nodes.Mod_loc=[];
Nodes.Mod_clust=[];
for clust=[1 2 4 5 6 7]%1:length(fieldnames(clust_f20_CL7_cleaned_cell))

    
        
 idx_temp=vertcat(clust_f20_CL7_cleaned_cell.(clustersF{clust}),(idx_f60_adjust(clust_f60_CL7_cleaned_cell.(clustersF{clust})))',(idx_s20_adjust(clust_s20_CL7_cleaned_cell.(clustersF{clust})))',(idx_s60_adjust(clust_s60_CL7_cleaned_cell.(clustersF{clust})))');   
    
options = statset('UseParallel',1); [idxKmeans_ROIs_f20 Cmap_ROIs_f20]=kmeans(ROI_temp2_all(idx_temp,:),moduleN(clust),'Options',options,'Distance','cityblock','Replicates',5,'MaxIter',1000,'Display','final');

Nodes.Mod_loc=vertcat(Nodes.Mod_loc,Cmap_ROIs_f20);
Nodes.Mod_clust=vertcat(Nodes.Mod_clust,ones(size(Cmap_ROIs_f20,1),1)*clust);

% temp_module_means=[];
% for i=1:size(Cmap_ROIs_f20,1)
%     temp_idx=find(idxKmeans_ROIs_f20==i);
%     temp_module_means(i,:)=mean(ZS_f20(idx_temp(temp_idx),:));
%     
% end
% Nodes.Mod_mean=vertcat(Nodes.Mod_mean,temp_module_means);

end

%%% to check it worked
clust=6;
idx_temp=vertcat(clust_f20_CL7_cleaned_cell.(clustersF{clust}),(idx_f60_adjust(clust_f60_CL7_cleaned_cell.(clustersF{clust})))',(idx_s20_adjust(clust_s20_CL7_cleaned_cell.(clustersF{clust})))',(idx_s60_adjust(clust_s60_CL7_cleaned_cell.(clustersF{clust})))');    
figure;scatter(ROI_temp2_all(idx_temp,1),ROI_temp2_all(idx_temp,2));


figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_loc(:,3));
gscatter(Nodes.Mod_loc(:,1),Nodes.Mod_loc(:,2),Nodes.Mod_clust);
view(-90,90);


%%% need to get ROIs of each fish per dataset


%% for f20
%load('final_F20_step1.mat','idx_Fish_f20','ZS_f20');
%%% OR 
load('ZS_N_Fish_all.mat')

%%% to quickly check what are the distances among nodes
for clust=[1 2 4 5 6 7]

d=pdist(Nodes.Mod_loc(find(Nodes.Mod_clust==clust),:));
quantile(d,[0.025 0.25 0.50 0.75 0.975])
(min(d)/2)-1

end

  
for clust=[1 2 4 5 6 7]
idx_temp=clust_f20_CL7_cleaned_cell.(clustersF{clust});
D = pdist2(ROI_temp2_all(idx_temp,:),Nodes.Mod_loc); 
d=pdist(Nodes.Mod_loc(find(Nodes.Mod_clust==clust),:));
%temp_min=(min(d)/2)-1;
temp_min=quantile(d,[0.025]);

row=[];
for node=find(Nodes.Mod_clust==clust)
[temp_row,~]=find(D(:,node)<temp_min);

%%% in this bit I am tryting to fix the overlaping. 
 
if ~isempty(find(ismember(temp_row,row)))
 temp=find(ismember(temp_row,row));
 temp_row(temp)=[];
 row=vertcat(row,temp_row);
 
else
 row=vertcat(row,temp_row);
end
 
 %%% to check they dont overlap. 
B = unique(row); % which will give you the unique elements of A in array B
Ncount = histc(row, B); % this willgive the number of occurences of each unique element
unique(Ncount)
 
 
end

%figure;scatter3(ROI_temp2_all(idx_temp(row),1),ROI_temp2_all(idx_temp(row),2),ROI_temp2_all(idx_temp(row),3),'filled');
   
%  B = unique(row); % which will give you the unique elements of A in array B
%  Ncount = histc(row, B); % this willgive the number of occurences of each unique element
%  unique(Ncount)

end

%% getting the ROIs of each node per fish
Nodes_f20=struct;
fish=unique(idx_Fish_f20);

for clust=[1 2 4 5 6 7]
   d=pdist(Nodes.Mod_loc(find(Nodes.Mod_clust==clust),:));
    %temp_min=(min(d)/2)-1;
    temp_min=quantile(d,[0.025]);

    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f20_CL7_cleaned_cell.(clustersF{clust}));

D = pdist2(ROI_temp2_all(temp_idx_clust,:),Nodes.Mod_loc); 
    

    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    [temp_row,~]=find(D(:,node)<temp_min);
    
    %%% in this bit I am tryting to fix the overlaping. 
        if ~isempty(find(ismember(temp_row,row)))
        temp=find(ismember(temp_row,row));
        temp_row(temp)=[];
        row=vertcat(row,temp_row);
 
        else
        row=vertcat(row,temp_row);
        end
    
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx_clust(temp_row);
    
    end

   
    Nodes_f20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish))))=idx_ROIs_Node;
      
    
end
end

%%% to check if i fixed the overlaping. it seemed that it worked. 
for clust=[1 2 4 5 6 7]
   
for tempfish=1:length(fish)
      

    idx_ROIs_Node=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    
    temp_idx=Nodes_f20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).idx;
    
    idx_ROIs_Node=vertcat(idx_ROIs_Node,temp_idx);
    end
     
    B = unique(idx_ROIs_Node); % which will give you the unique elements of A in array B
    Ncount = histc(idx_ROIs_Node, B); % this willgive the number of occurences of each unique element
    unique(Ncount)
    
    
end
end



%%% to get the means
for clust=[1 2 4 5 6 7]
   
    
for tempfish=1:length(fish)
   

    for node=(find(Nodes.Mod_clust==clust))'
    
   
    temp_idx=Nodes_f20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_f20(temp_idx,:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_f20(temp_idx,:));
        %plot(temp_mean);
    end
    
    Nodes_f20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(100,1344);
    for clust=[1 2 4 5 6 7]
        
        for node=(find(Nodes.Mod_clust==clust))'
            
        temp_mean=Nodes_f20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes_f20.mean_matrix.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes_f20.corr_matrix.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes_f20.corr_matrix.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);

%% now for f60


%load('final_F60_step1_2.mat','idx_Fish_f60');


%% getting the ROIs of each node per fish
Nodes_f60=struct;
fish=unique(idx_Fish_all(idx_f60_adjust));

fish(find(fish==47))=[]; %%% cause I also took out fish 47...


for clust=[1 2 4 5 6 7]
   d=pdist(Nodes.Mod_loc(find(Nodes.Mod_clust==clust),:));
    %temp_min=(min(d)/2)-1;
    temp_min=quantile(d,[0.025]);

    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_all(idx_f60_adjust)==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_f60_CL7_cleaned_cell.(clustersF{clust}));

D = pdist2(ROI_temp2_all(temp_idx_clust,:),Nodes.Mod_loc); 
    

    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    [temp_row,~]=find(D(:,node)<temp_min);
    
    %%% in this bit I am tryting to fix the overlaping. 
        if ~isempty(find(ismember(temp_row,row)))
        temp=find(ismember(temp_row,row));
        temp_row(temp)=[];
        row=vertcat(row,temp_row);
 
        else
        row=vertcat(row,temp_row);
        end
    
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx_clust(temp_row);
    
    end

   
    Nodes_f60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish))))=idx_ROIs_Node;
      
    
end
end

%%% to check if i fixed the overlaping. it seemed that it worked. 
for clust=[1 2 4 5 6 7]
   
for tempfish=1:length(fish)
      

    idx_ROIs_Node=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    
    temp_idx=Nodes_f60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).idx;
    
    idx_ROIs_Node=vertcat(idx_ROIs_Node,temp_idx);
    end
     
    B = unique(idx_ROIs_Node); % which will give you the unique elements of A in array B
    Ncount = histc(idx_ROIs_Node, B); % this willgive the number of occurences of each unique element
    unique(Ncount)
    
    
end
end



%%% to get the means
for clust=[1 2 4 5 6 7]
   
    
for tempfish=1:length(fish)
   

    for node=(find(Nodes.Mod_clust==clust))'
    
   
    temp_idx=Nodes_f60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(idx_f60_adjust(temp_idx),:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(idx_f60_adjust(temp_idx),:));
        %plot(temp_mean);
    end
    
    Nodes_f60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(100,1344);
    for clust=[1 2 4 5 6 7]
        
        for node=(find(Nodes.Mod_clust==clust))'
            
        temp_mean=Nodes_f60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes_f60.mean_matrix.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes_f60.corr_matrix.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes_f60.corr_matrix.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);




%% now for s20


%% getting the ROIs of each node per fish
Nodes_s20=struct;
fish=unique(idx_Fish_all(idx_s20_adjust));



for clust=[1 2 4 5 6 7]
   d=pdist(Nodes.Mod_loc(find(Nodes.Mod_clust==clust),:));
    %temp_min=(min(d)/2)-1;
    temp_min=quantile(d,[0.025]);

    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_all(idx_s20_adjust)==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s20_CL7_cleaned_cell.(clustersF{clust}));

D = pdist2(ROI_temp2_all(temp_idx_clust,:),Nodes.Mod_loc); 
    

    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    [temp_row,~]=find(D(:,node)<temp_min);
    
    %%% in this bit I am tryting to fix the overlaping. 
        if ~isempty(find(ismember(temp_row,row)))
        temp=find(ismember(temp_row,row));
        temp_row(temp)=[];
        row=vertcat(row,temp_row);
 
        else
        row=vertcat(row,temp_row);
        end
    
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx_clust(temp_row);
    
    end

   
    Nodes_s20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish))))=idx_ROIs_Node;
      
    
end
end

%%% to check if i fixed the overlaping. it seemed that it worked. 
for clust=[1 2 4 5 6 7]
   
for tempfish=1:length(fish)
      

    idx_ROIs_Node=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    
    temp_idx=Nodes_s20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).idx;
    
    idx_ROIs_Node=vertcat(idx_ROIs_Node,temp_idx);
    end
     
    B = unique(idx_ROIs_Node); % which will give you the unique elements of A in array B
    Ncount = histc(idx_ROIs_Node, B); % this willgive the number of occurences of each unique element
    unique(Ncount)
    
    
end
end



%%% to get the means
for clust=[1 2 4 5 6 7]
   
    
for tempfish=1:length(fish)
   

    for node=(find(Nodes.Mod_clust==clust))'
    
   
    temp_idx=Nodes_s20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(idx_s20_adjust(temp_idx),:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(idx_s20_adjust(temp_idx),:));
        %plot(temp_mean);
    end
    
    Nodes_s20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(100,1344);
    for clust=[1 2 4 5 6 7]
        
        for node=(find(Nodes.Mod_clust==clust))'
            
        temp_mean=Nodes_s20.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes_s20.mean_matrix.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes_s20.corr_matrix.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes_s20.corr_matrix.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);

%% now for s60

%% getting the ROIs of each node per fish
Nodes_s60=struct;
fish=unique(idx_Fish_all(idx_s60_adjust));



for clust=[1 2 4 5 6 7]
   d=pdist(Nodes.Mod_loc(find(Nodes.Mod_clust==clust),:));
    %temp_min=(min(d)/2)-1;
    temp_min=quantile(d,[0.025]);

    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_all(idx_s60_adjust)==fish(tempfish));
temp_idx_clust=intersect(temp_idx_fish,clust_s60_CL7_cleaned_cell.(clustersF{clust}));

D = pdist2(ROI_temp2_all(temp_idx_clust,:),Nodes.Mod_loc); 
    

    idx_ROIs_Node=struct;
    row=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    [temp_row,~]=find(D(:,node)<temp_min);
    
    %%% in this bit I am tryting to fix the overlaping. 
        if ~isempty(find(ismember(temp_row,row)))
        temp=find(ismember(temp_row,row));
        temp_row(temp)=[];
        row=vertcat(row,temp_row);
 
        else
        row=vertcat(row,temp_row);
        end
    
    
    idx_ROIs_Node.(strcat('node_',num2str(node))).idx=temp_idx_clust(temp_row);
    
    end

   
    Nodes_s60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish))))=idx_ROIs_Node;
      
    
end
end

%%% to check if i fixed the overlaping. it seemed that it worked. 
for clust=[1 2 4 5 6 7]
   
for tempfish=1:length(fish)
      

    idx_ROIs_Node=[];
    node=find(Nodes.Mod_clust==clust);
    for node=(find(Nodes.Mod_clust==clust))'
    
    temp_idx=Nodes_s60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).idx;
    
    idx_ROIs_Node=vertcat(idx_ROIs_Node,temp_idx);
    end
     
    B = unique(idx_ROIs_Node); % which will give you the unique elements of A in array B
    Ncount = histc(idx_ROIs_Node, B); % this willgive the number of occurences of each unique element
    unique(Ncount)
    
    
end
end



%%% to get the means
for clust=[1 2 4 5 6 7]
   
    
for tempfish=1:length(fish)
   

    for node=(find(Nodes.Mod_clust==clust))'
    
   
    temp_idx=Nodes_s60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).idx;
    if  isempty(temp_idx)
         temp_mean=NaN;
        
    elseif size(temp_idx,1)==1
        
        temp_mean=ZS_CN(idx_s60_adjust(temp_idx),:);
        %plot(temp_mean);
    else
        temp_mean=mean(ZS_CN(idx_s60_adjust(temp_idx),:));
        %plot(temp_mean);
    end
    
    Nodes_s60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).mean=temp_mean;
    
    
   end
    
end
end

%%% making means matrices per fish and testing a correlation
for tempfish=1:length(fish)

    temp_matrix=NaN(100,1344);
    for clust=[1 2 4 5 6 7]
        
        for node=(find(Nodes.Mod_clust==clust))'
            
        temp_mean=Nodes_s60.(clustersF{clust}).(strcat('fish_',num2str(fish(tempfish)))).(strcat('node_',num2str(node))).mean;
        
        if isnan(temp_mean) 
        else
          temp_matrix(node,:)= temp_mean; 
        end
        end
    end
        Nodes_s60.mean_matrix.(strcat('fish_',num2str(fish(tempfish))))=temp_matrix;
        
        R_temp=corrcoef(temp_matrix');
        
        Nodes_s60.corr_matrix.(strcat('fish_',num2str(fish(tempfish))))=R_temp;
        
end

Matrix_mean=[];
for tempfish=1:length(fish)
 Matrix_mean=cat(3,Matrix_mean,Nodes_s60.corr_matrix.(strcat('fish_',num2str(fish(tempfish)))));   
    
end

Matrix_mean=nanmean(Matrix_mean,3);


%% saving

Nodes.f20=Nodes_f20;
Nodes.f60=Nodes_f60;
Nodes.s20=Nodes_s20;
Nodes.s60=Nodes_s60;



save('Nodes_N_means_alldatasets.mat','Nodes','ROI_temp2_all','Zbrain_brainMask2D');


