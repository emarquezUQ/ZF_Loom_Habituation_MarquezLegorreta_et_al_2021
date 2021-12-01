
%%%% this script is for the tsne visualization for Marquez-Legorreta et al
%%%% 2021 used in Supplementary figure 3e and Supplementary figure 4



f20_cleaned_idxs=load('f20_cleaned_idxs.mat');

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');

load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');


%%% for all fish together
fishgoodrsq=f20_cleaned_idxs.idx_rsq_test_f20short_cleaned;


%%

clustersF=fieldnames(mean_CL4_f20);

clust_f20_CL4_cleaned=f20_cleaned_idxs.clust_f20_CL4_cleaned;

clust_f20_CL4_cleaned_cell={};
clust=fieldnames(clust_f20_CL4_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL4_cleaned_cell.(clustersF{j,1})=clust_f20_CL4_cleaned.(clust{j});   
end    

clusterlabels=zeros(size(fishgoodrsq));
counter=1;

    for clust=1:length(clustersF)
    
    temp_idx_clust=intersect(fishgoodrsq,clust_f20_CL4_cleaned_cell.(clustersF{clust}));
   
   idx_temp2=find(ismember(fishgoodrsq,temp_idx_clust));
   
   clusterlabels(idx_temp2)=clust;
    end

%%

clustersF=fieldnames(mean_CL7_f20);

clust_f20_CL7_cleaned=f20_cleaned_idxs.clust_f20_CL7_cleaned;

clust_f20_CL7_cleaned_cell={};
clust=fieldnames(clust_f20_CL7_cleaned);
for j=1:size(clustersF,1)
 clust_f20_CL7_cleaned_cell.(clustersF{j,1})=clust_f20_CL7_cleaned.(clust{j});   
end    

clusterlabels2=zeros(size(fishgoodrsq));
counter=1;
  
    for clust=1:length(clustersF)
    
    temp_idx_clust=intersect(fishgoodrsq,clust_f20_CL7_cleaned_cell.(clustersF{clust}));
   
   idx_temp2=find(ismember(fishgoodrsq,temp_idx_clust));
   
   clusterlabels2(idx_temp2)=clust;
    end

%%


%%%% After optimizing the tsne
%%%% it is recommended that perplexity~sqrt(N), so for f20 N=33903, then perplexity ~184
%%%% I doubled the exaggeration and increased the MaxIter

options = statset('MaxIter',3000);


Y_all6 = tsne(ZS_f20(fishgoodrsq,:),'Perplexity',184,'Exaggeration',40,'Distance','correlation','Options',options);


figure; gscatter(Y_all6(:,1),Y_all6(:,2),clusterlabels);
figure; gscatter(Y_all6(:,1),Y_all6(:,2),clusterlabels2);


%% checking density


%%%%% checking density of t-sne tests
%%%%% I am using the tool from https://au.mathworks.com/matlabcentral/fileexchange/8430-flow-cytometry-data-reader-and-visualization?focused=6779476&tab=function
%%%%% using the function as "dscatter_EML"


%%% saved results from original run. 
load('tsne_f20_extra.mat');

figure; gscatter(Y_all6(:,1),Y_all6(:,2),clusterlabels,'krgbm','.',1);
figure; gscatter(Y_all6(:,1),Y_all6(:,2),clusterlabels2,'rgkcybm','.',1);



figure; dscatter_EML(Y_all6(:,1),Y_all6(:,2));colorbar;

%%%% for this figure it is needed to have in the path the folder with the inferno
%%%% colormap
figure; dscatter_EML(Y_all6(:,1),Y_all6(:,2));colorbar;colormap(inferno);
figure; dscatter_EML(Y_all6(:,1),Y_all6(:,2),'MARKER','o');colorbar;colormap(inferno);
%figure; 
%hold on;
figure;
dscatter_EML(Y_all6(:,1),Y_all6(:,2),'PLOTTYPE','contour');colorbar;colormap(inferno); 

figure;dscatter_EML(Y_all6(:,1),Y_all6(:,2),'PLOTTYPE','image');colorbar;colormap(inferno); 
 
%%%% NOTE: to also have a figure with coloured areas I change the dscatter
%%%% function so it used "contourf" instead of "contour". 

