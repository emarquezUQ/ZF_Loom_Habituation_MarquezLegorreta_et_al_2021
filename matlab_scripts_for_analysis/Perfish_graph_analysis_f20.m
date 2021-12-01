
%%%%% trying to do a perfish graph analysis with the f20 dataset


%%%% this script is to make a graph with the ROIs of one single fish for
%%%% the first, 10th and 11th looms and compare them. I will separate the clusters.


%%
%%%% some things I need

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

%load('All_More_BrainReg.mat','PerBrainRegions');
 load('All_More_BrainReg2.mat')
%load('zbrain3D.mat')
 
load('Nodes_N_means_alldatasets.mat','ROI_temp2_all','Zbrain_brainMask2D');

load('ZS_N_Fish_all.mat','idx_Fish_all','idx_f60_adjust','idx_s20_adjust','idx_s60_adjust');

%RegionList={'Pallium','Subpallium','Thalamus','Pretectum','Tectum','Hindbrain'};

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};


f20_cleaned_idxs=load('f20_cleaned_idxs.mat');
f60_cleaned_idxs=load('f60_cleaned_idxs.mat');
s20_cleaned_idxs=load('s20_cleaned_idxs.mat');
s60_cleaned_idxs=load('s60_cleaned_idxs.mat');

load('Max_response_perBrain_all_CL4_corrected_distribution.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx');


%%%% I will also compare it with the original analysis. 

OrginalMatAll_corrected=load('graph_analysis_loomhab3.mat');

OriginalData_corrMat=load('graphs_FnS_all.mat');



%% for f20
%%%for f20. fish1
load('ZS_N_idx_Fish_all.mat','idx_Fish_f20');
%%% and 
load('ZS_N_Fish_all.mat')

fish_list=unique(idx_Fish_f20);

fish=fieldnames(OriginalData_corrMat.Data_corrMat2.f20);
fish(length(fish))=[];


%%%I am using the f20 dataset to test
%%% for CL7

%%% to order the clusters per name
clust_f20_CL7_cleaned=f20_cleaned_idxs.clust_f20_CL7_cleaned;

clust_f20_CL7_cleaned_cell={};
clust=fieldnames(clust_f20_CL7_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL7_cleaned_cell.(clustersF{j,1})=clust_f20_CL7_cleaned.(clust{j});   
end   




 %%% I need to order them. based on ClustersF.
 %%% fasthabs first and no sound
 goodorder_clust=[4 2 5 6 1 7];

ClustID_f20_SingleFish=struct;

for f=1:length(unique(idx_Fish_f20))

ClustID_f20_SingleFish.(fish{f}).idx=[];
ClustID_f20_SingleFish.(fish{f}).clustID_CL7=[];
ClustID_f20_SingleFish.(fish{f}).clustID_CL4=[];
for clust=1:length(goodorder_clust)
    idx_temp=clust_f20_CL7_cleaned_cell.(clustersF{goodorder_clust(clust)});
    idx_temp=intersect(find(idx_Fish_f20==fish_list(f)),idx_temp);
    temp_id=ones(size(idx_temp))*clust;
    
    if clust==1 || clust==2 || clust==3
       temp_id=ones(size(idx_temp))*1; 
    else
        temp_id=ones(size(idx_temp))*clust;
    end
    
    ClustID_f20_SingleFish.(fish{f}).idx=vertcat(ClustID_f20_SingleFish.(fish{f}).idx,idx_temp);
    ClustID_f20_SingleFish.(fish{f}).clustID_CL7=vertcat(ClustID_f20_SingleFish.(fish{f}).clustID_CL7,temp_id);
    ClustID_f20_SingleFish.(fish{f}).clustID_CL4=vertcat(ClustID_f20_SingleFish.(fish{f}).clustID_CL4,temp_id);
end
    
    
%%
figure;
%patch(brain3D,'EdgeColor','none','FaceAlpha',0.1);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
%figure; 
%scatter3(Modules.Mod_loc(:,1),Modules.Mod_loc(:,2),Modules.Mod_loc(:,3));
gscatter(ROI_temp2_all(ClustID_f20_SingleFish.(fish{f}).idx,1),ROI_temp2_all(ClustID_f20_SingleFish.(fish{f}).idx,2),ClustID_f20_SingleFish.(fish{f}).clustID_CL4,'gbrm');
view(-90,90);

end
%%

%% getting the corrmatrix for each loom and each fish

Data_corrMat_singleFish=struct;


for data=1%:4
    datatemp=datasets(data,:);
    
    
           
    for tempfish=1:length(fish)
        temp_mean=ZS_CN(ClustID_f20_SingleFish.(fish{tempfish}).idx,:);
        
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
        
    Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR=temp_R;
        
    end
             
end


counter=1;
figure;
for tempfish=randperm(length(fish),5)


for data=1%:4
    datatemp=datasets(data,:);
       
subplot(5,8,counter);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,1});caxis([-1 1]);colormap('parula'); %% for pre loom
subplot(5,8,counter+1);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,2});caxis([-1 1]); colormap('parula');%% for 1st loom
subplot(5,8,counter+2);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,3});caxis([-1 1]);colormap('parula');
subplot(5,8,counter+3);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,4});caxis([-1 1]);colormap('parula');
subplot(5,8,counter+4);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,5});caxis([-1 1]);colormap('parula');
subplot(5,8,counter+5);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,6});caxis([-1 1]);colormap('parula');
subplot(5,8,counter+6);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,11});caxis([-1 1]);colormap('parula');%% for 10th loom
subplot(5,8,counter+7);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,12});caxis([-1 1]);colormap('parula'); %% for 11th loom
title(datatemp);

counter=counter+8;
end
end

%%

%% testing_graph_analysis_with_BCT3



%%
%%%% I will look at density, the number of degrees and strength.

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_singleFish=struct;


for data=1%:4
    datatemp=datasets(data,:);
    
     
           
    for tempfish=1:length(fish)
    
    
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,moment(m)}),0.75);
     MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
    end

end
    
%% density
%


%kden = density_und(CIJ);

for data=1%:4
    datatemp=datasets(data,:);
    
    
    for tempfish=1:length(fish)
        
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp_kden=density_und(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).Mat);
            
            MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).kden=temp_kden;
        end
    end
end

figure;
temp=[];
temp_mean=[];
for data=1%:4
    datatemp=datasets(data,:);
    
    for tempfish=1:length(fish)
        
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp(1,i)=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).kden;
        end
        
        
        temp_mean(tempfish,:)=temp;
        
        plot(temp);
        hold on;
        
    end
    
    plot(nanmean(temp_mean),'r','LineWidth',3);
    hold on;
    
end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp,'k','LineWidth',3);
hold on;

end


%%%% degrees and Strength

for data=1%:4
    datatemp=datasets(data,:);
    
    
    
    for tempfish=1:length(fish)
        
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            deg=degrees_und(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).Mat);
            str=strengths_und_sign(abs(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).Mat));
            
            MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).deg=deg;
            MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).str=str;
        end
    end
end

figure;
temp_mean=[];
for data=1%:4
    datatemp=datasets(data,:);
    
    
    for tempfish=1:length(fish)
        temp=[];
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp(:,i)=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).deg;
        end
        
        temp_mean(tempfish,:)=nanmean(temp);
        
        plot(nanmean(temp));
        hold on;
        
    end
    
    plot(nanmean(temp_mean),'r','LineWidth',3);
        hold on;
end

temp=[];
for data=1%:4
    datatemp=datasets(data,:);
    loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));
    
    for i=1:length(loom)
        
        temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).deg;
    end
    plot(nanmean(temp),'k','LineWidth',3);
    hold on;
    
end


%%%% participation and cluster coeficient


for data=1%:4
    datatemp=datasets(data,:);
    
    
    for tempfish=1:length(fish)
        
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp_mat=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).Mat;
            temp_mat(isnan(temp_mat))=0;
            P=participation_coef(temp_mat,ClustID_f20_SingleFish.(fish{tempfish}).clustID_CL4);
            Ccoef=clustering_coef_wu(temp_mat);
            
            MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).Ccoef=Ccoef;
            MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).P=P;
            
        end
    end
end

figure;
temp_mean=[];
for data=1%:4
    datatemp=datasets(data,:);
    
    
    for tempfish=1:length(fish)
        temp=[];
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp(:,i)=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).Ccoef;
        end
        
        temp_mean(tempfish,:)=nanmean(temp);
        
        plot(nanmean(temp));
        hold on;
        
    end
    plot(nanmean(temp_mean),'r','LineWidth',3);
        hold on;
    
end

temp=[];
for data=1%:4
    datatemp=datasets(data,:);
    loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));
    
    for i=1:length(loom)
        
        temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).Ccoef;
    end
    plot(nanmean(temp),'k','LineWidth',3);
    hold on;
        
end

%%

save('Perfish_Graph_analysis_f20.mat','ClustID_f20_SingleFish','Data_corrMat_singleFish','MatAll_corrected_singleFish');



