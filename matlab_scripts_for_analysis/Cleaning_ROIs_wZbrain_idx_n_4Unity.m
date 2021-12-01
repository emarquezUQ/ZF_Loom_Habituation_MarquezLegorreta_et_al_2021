
%%%% this script is to clean up that are not in the brain. 
%%%% First to have a cleaned datasets for rasters and means and further analysis.
%%%% a second part is to get the coordiantes and write cvs file to plot in
%%%% Unity


%%% i first make a mask with all the brain regions
load('Zbrain_Masks.mat');


Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
 
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

figure;scatter(Zbrain_AllMask(:,1),Zbrain_AllMask(:,2),'.');

%scatter3(Zbrain_AllMask(:,1),Zbrain_AllMask(:,2),Zbrain_AllMask(:,3),'.');



%%% i need to load all the idx and ROI coordinates

%%

%%% first for f20


%%%% first you need to load what you need

load('BrainReg_F20.mat','ROI_temp2')
load('final_F20_step1.mat','gooodmaps','rawregressF20');
load('final_F20_step1.mat','idx_Fish_f20','idx_Plane_f20','idx_rsq_test_f20short','High_corr_Nb_f20','High_corr_Nb_f20_short');

% figure;
% hold on;
scatter(ROI_temp2(idx_rsq_test_f20short,1),ROI_temp2(idx_rsq_test_f20short,2));

%Removing the eyes

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

 
idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_test_f20short),1),ROI_temp2(intersect(idx_brain2,idx_rsq_test_f20short),2));

idx_rsq_test_f20short_cleaned=intersect(idx_brain2,idx_rsq_test_f20short);




%%% To get the ROIs of each cluster. First with the raw 7.
clust_f20_CL7_cleaned={};
for i=1:length(unique(High_corr_Nb_f20))

    idx_temp1=find(High_corr_Nb_f20==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_f20short(idx_temp1);
    
    filename=strcat('clust_f20_CL7_',num2str(i),'_cleaned');
    
    clust_f20_CL7_cleaned.(filename)=intersect(idx_brain2,idx_temp2);
  

end

%%% to check 
%figure;plot(mean(ZS_f20(idx_temp2,ZS_short_F20)));


%%% To get the ROIs of each cluster. now with the 4 main ones.
clust_f20_CL4_cleaned={};
for i=gooodmaps

    idx_temp1=find(High_corr_Nb_f20_short==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_f20short(idx_temp1);
    
    filename=strcat('clust_f20_CL4_',num2str(i),'_cleaned');
    
    clust_f20_CL4_cleaned.(filename)=intersect(idx_brain2,idx_temp2);


end

%%% to check 
% idx_temp1=find(High_corr_Nb_f20_short==2);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)
% idx_temp2=idx_rsq_test_f20short(idx_temp1);
% figure;plot(mean(ZS_f20(idx_temp2,ZS_short_F20)));
% figure;scatter(ROI_temp2(idx_temp2,1),ROI_temp2(idx_temp2,2));


%%% now for the tail movements


load('Tail_mov_F20.mat','idx_rsq_Mov');


% hold on;
scatter(ROI_temp2(idx_rsq_Mov,1),ROI_temp2(idx_rsq_Mov,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),1),ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),2));



%%% To get the ROIs of of the movements 

    
    idx_rsq_Mov_cleaned=intersect(idx_brain2,idx_rsq_Mov);

    
    
    
    %%% and for the multisense
    
    
    load('multisense_F20.mat','idx_multisense','High_corr_multisense_Nb_f20');

% hold on;
scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_multisense),1),ROI_temp2(intersect(idx_brain2,idx_multisense),2));
    
    
    idx_multisense_cleaned=intersect(idx_brain2,idx_multisense);
    
    
    %%% To get the ROIs of multisense subtype.
    multisense_f20_clusters={};
for i=1:length(unique(High_corr_multisense_Nb_f20))

    idx_temp1=find(High_corr_multisense_Nb_f20==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_multisense(idx_temp1);
    
    filename=strcat('multisense_f20_cluster_',num2str(i),'_cleaned');
    
    multisense_f20_clusters.(filename)=intersect(idx_brain2,idx_temp2);

 

end
    
save('f20_cleaned_idxs.mat','idx_rsq_test_f20short_cleaned','clust_f20_CL7_cleaned','clust_f20_CL4_cleaned','idx_rsq_Mov_cleaned','idx_multisense_cleaned','multisense_f20_clusters');

clear all


%%

%%% for f60



load('final_F60_step1_2.mat','idx_Fish_f60','idx_Plane_f60','ZS_short_F60');
load('final_F60_step1_3_short_correction.mat');
load('final_F20_step1.mat','gooodmaps','rawregressF20');
load('BrainReg_F60.mat','ROI_temp2');




% figure;
% hold on;
scatter(ROI_temp2(idx_rsq_test_f60short3,1),ROI_temp2(idx_rsq_test_f60short3,2));

%Removing the eyes

 

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_test_f60short3),1),ROI_temp2(intersect(idx_brain2,idx_rsq_test_f60short3),2));


idx_rsq_test_f60short3_cleaned=intersect(idx_brain2,idx_rsq_test_f60short3);


%%% To get the ROIs of each cluster. First with the raw 7.
clust_f60_CL7_cleaned={};
for i=1:length(unique(High_corr_Nb_f60_3))

    idx_temp1=find(High_corr_Nb_f60_3==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_f60short3(idx_temp1);
    
    filename=strcat('clust_f60_CL7_',num2str(i),'_cleaned');
    
    clust_f60_CL7_cleaned.(filename)=intersect(idx_brain2,idx_temp2);
  

end



%%% to check 
%figure;plot(mean(ZS_f60(idx_temp2,ZS_short_F60)));


%%% To get the ROIs of each cluster. now iwht the 4 main ones.
clust_f60_CL4_cleaned={};
for i=gooodmaps

    idx_temp1=find(High_corr_Nb_f60_short==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_f60short3(idx_temp1);
    
    filename=strcat('clust_f60_CL4_',num2str(i),'_cleaned');
    
    clust_f60_CL4_cleaned.(filename)=intersect(idx_brain2,idx_temp2);


end
%%% to check 
% idx_temp1=find(High_corr_Nb_f60_short==2);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)
% idx_temp2=idx_rsq_test_f60short3(idx_temp1);
% figure;plot(mean(ZS_f60(idx_temp2,ZS_short_F60)));
% figure;scatter(ROI_temp2(idx_temp2,1),ROI_temp2(idx_temp2,2));


%%% now for the tail movements


load('Tail_mov_F60.mat','idx_rsq_Mov');


% hold on;
scatter(ROI_temp2(idx_rsq_Mov,1),ROI_temp2(idx_rsq_Mov,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),1),ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),2));



%%% To get the ROIs of of the movements 

    
    idx_rsq_Mov_cleaned=intersect(idx_brain2,idx_rsq_Mov);
    
    
    
    %%% and for the multisense
    
    
    load('multisense_F60.mat','idx_multisense','High_corr_multisense_Nb_f60');

    
    % hold on;
scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_multisense),1),ROI_temp2(intersect(idx_brain2,idx_multisense),2));
    
    

idx_multisense_cleaned=intersect(idx_brain2,idx_multisense);
    
    
    %%% To get the ROIs of multisense subtype.
    multisense_f60_clusters={};
for i=1:length(unique(High_corr_multisense_Nb_f60))

    idx_temp1=find(High_corr_multisense_Nb_f60==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_multisense(idx_temp1);
    
    filename=strcat('multisense_f60_cluster_',num2str(i),'_cleaned');
    
    multisense_f60_clusters.(filename)=intersect(idx_brain2,idx_temp2);

 

end
    
save('f60_cleaned_idxs.mat','idx_rsq_test_f60short3_cleaned','clust_f60_CL7_cleaned','clust_f60_CL4_cleaned','idx_rsq_Mov_cleaned','idx_multisense_cleaned','multisense_f60_clusters');

clear all



%%

%%% first for s20




load('BrainReg_S20.mat','ROI_temp2')
load('final_S20_step1.mat','gooodmaps','rawregressS20');
load('final_S20_step1.mat','idx_Fish_s20','idx_Plane_s20','idx_rsq_test_s20short','High_corr_Nb_s20','High_corr_Nb_s20_short');



% figure;
% hold on;
scatter(ROI_temp2(idx_rsq_test_s20short,1),ROI_temp2(idx_rsq_test_s20short,2));

%Removing the eyes

 

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_test_s20short),1),ROI_temp2(intersect(idx_brain2,idx_rsq_test_s20short),2));

idx_rsq_test_s20short_cleaned=intersect(idx_brain2,idx_rsq_test_s20short);




%%% To get the ROIs of each cluster. First with the raw 7.
clust_s20_CL7_cleaned={};
for i=1:length(unique(High_corr_Nb_s20))

    idx_temp1=find(High_corr_Nb_s20==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_s20short(idx_temp1);
    
    filename=strcat('clust_s20_CL7_',num2str(i),'_cleaned');
    
    clust_s20_CL7_cleaned.(filename)=intersect(idx_brain2,idx_temp2);
  

end

%%% to check 
%figure;plot(mean(ZS_s20(idx_temp2,ZS_short_S20)));


%%% To get the ROIs of each cluster. now iwht the 4 main ones.
clust_s20_CL4_cleaned={};
for i=gooodmaps

    idx_temp1=find(High_corr_Nb_s20_short==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_s20short(idx_temp1);
    
    filename=strcat('clust_s20_CL4_',num2str(i),'_cleaned');
    
    clust_s20_CL4_cleaned.(filename)=intersect(idx_brain2,idx_temp2);


end

%%% to check 
% idx_temp1=find(High_corr_Nb_s20_short==2);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)
% idx_temp2=idx_rsq_test_s20short(idx_temp1);
% figure;plot(mean(ZS_s20(idx_temp2,ZS_short_S20)));
% figure;scatter(ROI_temp2(idx_temp2,1),ROI_temp2(idx_temp2,2));


%%% now for the tail movements


load('Tail_mov_S20.mat','idx_rsq_Mov');


% hold on;
scatter(ROI_temp2(idx_rsq_Mov,1),ROI_temp2(idx_rsq_Mov,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),1),ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),2));



%%% To get the ROIs of of the movements 

    
    idx_rsq_Mov_cleaned=intersect(idx_brain2,idx_rsq_Mov);

    
    
    
    %%% and for the multisense
    
    
    load('multisense_S20.mat','idx_multisense','High_corr_multisense_Nb_s20');

% hold on;
scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_multisense),1),ROI_temp2(intersect(idx_brain2,idx_multisense),2));
    
    
    idx_multisense_cleaned=intersect(idx_brain2,idx_multisense);
    
    
    %%% To get the ROIs of multisense subtype.
    multisense_s20_clusters={};
for i=1:length(unique(High_corr_multisense_Nb_s20))

    idx_temp1=find(High_corr_multisense_Nb_s20==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_multisense(idx_temp1);
    
    filename=strcat('multisense_s20_cluster_',num2str(i),'_cleaned');
    
    multisense_s20_clusters.(filename)=intersect(idx_brain2,idx_temp2);

 

end
    
save('s20_cleaned_idxs.mat','idx_rsq_test_s20short_cleaned','clust_s20_CL7_cleaned','clust_s20_CL4_cleaned','idx_rsq_Mov_cleaned','idx_multisense_cleaned','multisense_s20_clusters');

clear all



%%

%%% first for s60




load('BrainReg_S60.mat','ROI_temp2')
load('final_S60_step1.mat','gooodmaps','rawregressS60');
load('final_S60_step1.mat','idx_Fish_s60','idx_Plane_s60','idx_rsq_test_s60short','High_corr_Nb_s60','High_corr_Nb_s60_short');




% figure;
% hold on;
scatter(ROI_temp2(idx_rsq_test_s60short,1),ROI_temp2(idx_rsq_test_s60short,2));

%Removing the eyes

 

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_test_s60short),1),ROI_temp2(intersect(idx_brain2,idx_rsq_test_s60short),2));

idx_rsq_test_s60short_cleaned=intersect(idx_brain2,idx_rsq_test_s60short);




%%% To get the ROIs of each cluster. First with the raw 7.
clust_s60_CL7_cleaned={};
for i=1:length(unique(High_corr_Nb_s60))

    idx_temp1=find(High_corr_Nb_s60==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_s60short(idx_temp1);
    
    filename=strcat('clust_s60_CL7_',num2str(i),'_cleaned');
    
    clust_s60_CL7_cleaned.(filename)=intersect(idx_brain2,idx_temp2);
  

end

%%% to check 
%figure;plot(mean(ZS_s60(idx_temp2,ZS_short_S60)));


%%% To get the ROIs of each cluster. now iwht the 4 main ones.
clust_s60_CL4_cleaned={};
for i=gooodmaps

    idx_temp1=find(High_corr_Nb_s60_short==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_s60short(idx_temp1);
    
    filename=strcat('clust_s60_CL4_',num2str(i),'_cleaned');
    
    clust_s60_CL4_cleaned.(filename)=intersect(idx_brain2,idx_temp2);


end

%%% to check 
% idx_temp1=find(High_corr_Nb_s60_short==2);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)
% idx_temp2=idx_rsq_test_s60short(idx_temp1);
% figure;plot(mean(ZS_s60(idx_temp2,ZS_short_S60)));
% figure;scatter(ROI_temp2(idx_temp2,1),ROI_temp2(idx_temp2,2));


%%% now for the tail movements


load('Tail_mov_S60.mat','idx_rsq_Mov');


% hold on;
scatter(ROI_temp2(idx_rsq_Mov,1),ROI_temp2(idx_rsq_Mov,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),1),ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),2));



%%% To get the ROIs of of the movements 

    
    idx_rsq_Mov_cleaned=intersect(idx_brain2,idx_rsq_Mov);

    
    
    
    %%% and for the multisense
    
    
    load('multisense_S60.mat','idx_multisense','High_corr_multisense_Nb_s60');

% hold on;
scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_multisense),1),ROI_temp2(intersect(idx_brain2,idx_multisense),2));
    
    
    idx_multisense_cleaned=intersect(idx_brain2,idx_multisense);
    
    
    %%% To get the ROIs of multisense subtype.
    multisense_s60_clusters={};
for i=1:length(unique(High_corr_multisense_Nb_s60))

    idx_temp1=find(High_corr_multisense_Nb_s60==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_multisense(idx_temp1);
    
    filename=strcat('multisense_s60_cluster_',num2str(i),'_cleaned');
    
    multisense_s60_clusters.(filename)=intersect(idx_brain2,idx_temp2);

 

end
    
save('s60_cleaned_idxs.mat','idx_rsq_test_s60short_cleaned','clust_s60_CL7_cleaned','clust_s60_CL4_cleaned','idx_rsq_Mov_cleaned','idx_multisense_cleaned','multisense_s60_clusters');

clear all



%% for unity

%%%%% this second section is the same as above but saving the coordiantes for Unity 


%%% i first make a mask with all the brain regions
load('Zbrain_Masks.mat');


Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
 
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

figure;scatter(Zbrain_AllMask(:,1),Zbrain_AllMask(:,2),'.');

%scatter3(Zbrain_AllMask(:,1),Zbrain_AllMask(:,2),Zbrain_AllMask(:,3),'.');



%%% i need to load all the idx and ROI coordinates

%%

%%% first for f20


%%%% first you need to load what you need

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab

load('BrainReg_F20.mat','ROI_temp2')
load('final_F20_step1.mat','gooodmaps','rawregressF20');
load('final_F20_step1.mat','idx_Fish_f20','idx_Plane_f20','idx_rsq_test_f20short','High_corr_Nb_f20','High_corr_Nb_f20_short');




% figure;
% hold on;
scatter(ROI_temp2(idx_rsq_test_f20short,1),ROI_temp2(idx_rsq_test_f20short,2));

%Removing the eyes

 

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

%idxKmeans_final(~idx_brain)=0;  %%% this is how gilles had it but his only
                                 %%% work for indexes with all the ROIs. i need to adjust it.
 
idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_test_f20short),1),ROI_temp2(intersect(idx_brain2,idx_rsq_test_f20short),2));



%%% To get the ROIs of each cluster. First with the raw 7.
for i=1:length(unique(High_corr_Nb_f20))

    idx_temp1=find(High_corr_Nb_f20==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_f20short(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_f20_CL7_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end

%%% to check 
%figure;plot(mean(ZS_f20(idx_temp2,ZS_short_F20)));


%%% To get the ROIs of each cluster. now iwht the 4 main ones.
for i=gooodmaps

    idx_temp1=find(High_corr_Nb_f20_short==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_f20short(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_f20_CL4_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end

%%% to check 
% idx_temp1=find(High_corr_Nb_f20_short==2);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)
% idx_temp2=idx_rsq_test_f20short(idx_temp1);
% figure;plot(mean(ZS_f20(idx_temp2,ZS_short_F20)));
% figure;scatter(ROI_temp2(idx_temp2,1),ROI_temp2(idx_temp2,2));


%%% now for the tail movements


load('Tail_mov_F20.mat','idx_rsq_Mov');


% hold on;
scatter(ROI_temp2(idx_rsq_Mov,1),ROI_temp2(idx_rsq_Mov,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),1),ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),2));



%%% To get the ROIs of of the movements 

    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_Tail_mov_f20_','cleaned','.csv');

    csvwrite(filename,CSV_temp);
    
    
    
    %%% and for the multisense
    
    
    load('multisense_F20.mat','idx_multisense','High_corr_multisense_Nb_f20');

% hold on;
scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_multisense),1),ROI_temp2(intersect(idx_brain2,idx_multisense),2));
    
    
    
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_multisense),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_multisense_f20_','cleaned','.csv');

    csvwrite(filename,CSV_temp);
    
    
    %%% To get the ROIs of multisense subtype.
for i=1:length(unique(High_corr_multisense_Nb_f20))

    idx_temp1=find(High_corr_multisense_Nb_f20==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_multisense(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_multisense_f20_cluster_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end
    
clear all


%%

%%% for f60



load('final_F60_step1_2.mat','idx_Fish_f60','idx_Plane_f60','ZS_short_F60');
load('final_F60_step1_3_short_correction.mat');
load('final_F20_step1.mat','gooodmaps','rawregressF20');
load('BrainReg_F60.mat','ROI_temp2');




% figure;
% hold on;
scatter(ROI_temp2(idx_rsq_test_f60short3,1),ROI_temp2(idx_rsq_test_f60short3,2));

%Removing the eyes

 

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

%idxKmeans_final(~idx_brain)=0;  %%% this is how gilles had it but his only
                                 %%% work for indexes with all the ROIs. i need to adjust it.
 
idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_test_f60short3),1),ROI_temp2(intersect(idx_brain2,idx_rsq_test_f60short3),2));



%%% To get the ROIs of each cluster. First with the raw 7.
for i=1:length(unique(High_corr_Nb_f60_3))

    idx_temp1=find(High_corr_Nb_f60_3==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_f60short3(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_f60_CL7_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end

%%% to check 
%figure;plot(mean(ZS_f60(idx_temp2,ZS_short_F60)));


%%% To get the ROIs of each cluster. now iwht the 4 main ones.
for i=gooodmaps

    idx_temp1=find(High_corr_Nb_f60_short==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_f60short3(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_f60_CL4_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end

%%% to check 
% idx_temp1=find(High_corr_Nb_f60_short==2);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)
% idx_temp2=idx_rsq_test_f60short3(idx_temp1);
% figure;plot(mean(ZS_f60(idx_temp2,ZS_short_F60)));
% figure;scatter(ROI_temp2(idx_temp2,1),ROI_temp2(idx_temp2,2));


%%% now for the tail movements


load('Tail_mov_F60.mat','idx_rsq_Mov');


% hold on;
scatter(ROI_temp2(idx_rsq_Mov,1),ROI_temp2(idx_rsq_Mov,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),1),ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),2));



%%% To get the ROIs of of the movements 

    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_Tail_mov_f60_','cleaned','.csv');

    csvwrite(filename,CSV_temp);
    
    
    
    %%% and for the multisense
    
    
    load('multisense_F60.mat','idx_multisense','High_corr_multisense_Nb_f60');

    
    % hold on;
scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_multisense),1),ROI_temp2(intersect(idx_brain2,idx_multisense),2));
    
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_multisense),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_multisense_f60_','cleaned','.csv');

    csvwrite(filename,CSV_temp);
    
    
    %%% To get the ROIs of multisense subtype.
for i=1:length(unique(High_corr_multisense_Nb_f60))

    idx_temp1=find(High_corr_multisense_Nb_f60==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_multisense(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_multisense_f60_cluster_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end
    
clear all



%%

%%% first for s20




load('BrainReg_S20.mat','ROI_temp2')
load('final_S20_step1.mat','gooodmaps','rawregressS20');
load('final_S20_step1.mat','idx_Fish_s20','idx_Plane_s20','idx_rsq_test_s20short','High_corr_Nb_s20','High_corr_Nb_s20_short');




% figure;
% hold on;
scatter(ROI_temp2(idx_rsq_test_s20short,1),ROI_temp2(idx_rsq_test_s20short,2));

%Removing the eyes

 

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

%idxKmeans_final(~idx_brain)=0;  %%% this is how gilles had it but his only
                                 %%% work for indexes with all the ROIs. i need to adjust it.
 
idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_test_s20short),1),ROI_temp2(intersect(idx_brain2,idx_rsq_test_s20short),2));



%%% To get the ROIs of each cluster. First with the raw 7.
for i=1:length(unique(High_corr_Nb_s20))

    idx_temp1=find(High_corr_Nb_s20==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_s20short(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_s20_CL7_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end

%%% to check 
%figure;plot(mean(ZS_s20(idx_temp2,ZS_short_S20)));


%%% To get the ROIs of each cluster. now iwht the 4 main ones.
for i=gooodmaps

    idx_temp1=find(High_corr_Nb_s20_short==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_s20short(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_s20_CL4_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end

%%% to check 
% idx_temp1=find(High_corr_Nb_s20_short==2);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)
% idx_temp2=idx_rsq_test_s20short(idx_temp1);
% figure;plot(mean(ZS_s20(idx_temp2,ZS_short_S20)));
% figure;scatter(ROI_temp2(idx_temp2,1),ROI_temp2(idx_temp2,2));


%%% now for the tail movements


load('Tail_mov_S20.mat','idx_rsq_Mov');


% hold on;
scatter(ROI_temp2(idx_rsq_Mov,1),ROI_temp2(idx_rsq_Mov,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),1),ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),2));



%%% To get the ROIs of of the movements 

    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_Tail_mov_s20_','cleaned','.csv');

    csvwrite(filename,CSV_temp);
    
    
    
    %%% and for the multisense
    
    
    load('multisense_S20.mat','idx_multisense','High_corr_multisense_Nb_s20');

% hold on;
scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_multisense),1),ROI_temp2(intersect(idx_brain2,idx_multisense),2));
    
    
    
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_multisense),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_multisense_s20_','cleaned','.csv');

    csvwrite(filename,CSV_temp);
    
    
    %%% To get the ROIs of multisense subtype.
for i=1:length(unique(High_corr_multisense_Nb_s20))

    idx_temp1=find(High_corr_multisense_Nb_s20==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_multisense(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_multisense_s20_cluster_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end
    
clear all

%%

%%% first for s60




load('BrainReg_S60.mat','ROI_temp2')
load('final_S60_step1.mat','gooodmaps','rawregressS60');
load('final_S60_step1.mat','idx_Fish_s60','idx_Plane_s60','idx_rsq_test_s60short','High_corr_Nb_s60','High_corr_Nb_s60_short');




% figure;
% hold on;
scatter(ROI_temp2(idx_rsq_test_s60short,1),ROI_temp2(idx_rsq_test_s60short,2));

%Removing the eyes

 

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

%idxKmeans_final(~idx_brain)=0;  %%% this is how gilles had it but his only
                                 %%% work for indexes with all the ROIs. i need to adjust it.
 
idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_test_s60short),1),ROI_temp2(intersect(idx_brain2,idx_rsq_test_s60short),2));



%%% To get the ROIs of each cluster. First with the raw 7.
for i=1:length(unique(High_corr_Nb_s60))

    idx_temp1=find(High_corr_Nb_s60==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_s60short(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_s60_CL7_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end

%%% to check 
%figure;plot(mean(ZS_s60(idx_temp2,ZS_short_S60)));


%%% To get the ROIs of each cluster. now iwht the 4 main ones.
for i=gooodmaps

    idx_temp1=find(High_corr_Nb_s60_short==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_rsq_test_s60short(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_clust_s60_CL4_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end

%%% to check 
% idx_temp1=find(High_corr_Nb_s60_short==2);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)
% idx_temp2=idx_rsq_test_s60short(idx_temp1);
% figure;plot(mean(ZS_s60(idx_temp2,ZS_short_S60)));
% figure;scatter(ROI_temp2(idx_temp2,1),ROI_temp2(idx_temp2,2));


%%% now for the tail movements


load('Tail_mov_S60.mat','idx_rsq_Mov');


% hold on;
scatter(ROI_temp2(idx_rsq_Mov,1),ROI_temp2(idx_rsq_Mov,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),1),ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),2));



%%% To get the ROIs of of the movements 

    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_rsq_Mov),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_Tail_mov_s60_','cleaned','.csv');

    csvwrite(filename,CSV_temp);
    
    
    
    %%% and for the multisense
    
    
    load('multisense_S60.mat','idx_multisense','High_corr_multisense_Nb_s60');

% hold on;
scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2));


% hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_multisense),1),ROI_temp2(intersect(idx_brain2,idx_multisense),2));
    
    
    
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_multisense),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_multisense_s60_','cleaned','.csv');

    csvwrite(filename,CSV_temp);
    
    
    %%% To get the ROIs of multisense subtype.
for i=1:length(unique(High_corr_multisense_Nb_s60))

    idx_temp1=find(High_corr_multisense_Nb_s60==i);  %Or Whatever idx you?re using to identify the relevant clusters (idxKmeans_final for example)

    idx_temp2=idx_multisense(idx_temp1);
    
    CSV_temp=ROI_temp2(intersect(idx_brain2,idx_temp2),:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_multisense_s60_cluster_',num2str(i),'_cleaned','.csv');

    csvwrite(filename,CSV_temp);

end
    
clear all


