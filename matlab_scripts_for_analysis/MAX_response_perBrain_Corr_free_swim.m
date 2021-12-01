
%%%% this script is to get calculate correlations between the brain regions Max
%%%% responses and the behavioural responses from the free-swiming
%%%% animals. Doing it per fish so I can calculate confidence
%%%% intervals  

%%% it's following up with the dataset made in
%%% MAX_response_perBrain_all_CL4_corrected.m

load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx', 'S_trim');

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);

load('Zbrain_Masks.mat');
load('All_More_BrainReg.mat','PerBrainRegions');



load('Max_response_perBrain_all_CL4_corrected.mat');


%%% I will also try to do this per brain region. 

load('All_means_maxResp_CL4_normalized.mat');

%%%% reordering the RegionList

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Hindbrain','Cerebellum'};

RegionList2=RegionList([1 2 3 4 5 6 7 9 8]);
RegionList=RegionList2;


%%%% get the data from the behaviour (I took it from prism and made a new
%%%% variable)
Behaviour=[];

%%%% then put in the variable the values from the escape probability for each loom, per columns, f20,f60,s20,s60
%f20
Behavior(:,1)=[0.972200000000000;0.916700000000000;0.833300000000000;0.944400000000000;0.750000000000000;0.777800000000000;0.777800000000000;0.722200000000000;0.777800000000000;0.722200000000000;0.944400000000000;0.666700000000000;0.722200000000000;0.638900000000000;0.583300000000000;0.666700000000000;0.611100000000000;0.500000000000000;0.472200000000000;0.444400000000000;0.750000000000000;0.750000000000000;0.694400000000000;0.638900000000000;0.527800000000000;0.611100000000000;0.583300000000000;0.555600000000000;0.472200000000000;0.416700000000000];
%f60
Behavior(:,2)=[1;0.861100000000000;0.916700000000000;0.805600000000000;0.888900000000000;0.777800000000000;0.861100000000000;0.833300000000000;0.833300000000000;0.750000000000000;0.944400000000000;0.500000000000000;0.527800000000000;0.722200000000000;0.611100000000000;0.527800000000000;0.333300000000000;0.472200000000000;0.527800000000000;0.416700000000000;0.611100000000000;0.583300000000000;0.500000000000000;0.500000000000000;0.444400000000000;0.416700000000000;0.416700000000000;0.444400000000000;0.444400000000000;0.416700000000000];
%s20
Behavior(:,3)=[0.888900000000000;0.833300000000000;0.805600000000000;0.750000000000000;0.611100000000000;0.638900000000000;0.500000000000000;0.388900000000000;0.277800000000000;0.388900000000000;0.861100000000000;0.611100000000000;0.777800000000000;0.527800000000000;0.388900000000000;0.388900000000000;0.444400000000000;0.388900000000000;0.222200000000000;0.305600000000000;0.777800000000000;0.805600000000000;0.555600000000000;0.416700000000000;0.333300000000000;0.361100000000000;0.277800000000000;0.333300000000000;0.250000000000000;0.222200000000000];
%s60
Behavior(:,4)=[0.944400000000000;0.888900000000000;0.777800000000000;0.666700000000000;0.638900000000000;0.694400000000000;0.638900000000000;0.638900000000000;0.583300000000000;0.555600000000000;0.694400000000000;0.583300000000000;0.500000000000000;0.416700000000000;0.555600000000000;0.527800000000000;0.444400000000000;0.444400000000000;0.388900000000000;0.333300000000000;0.472200000000000;0.555600000000000;0.472200000000000;0.472200000000000;0.388900000000000;0.333300000000000;0.333300000000000;0.250000000000000;0.250000000000000;0.333300000000000];



Correlations_BrainVSbehav=struct;

for brain=1:length(RegionList)

counter=1;


for clust=1:3

if ismember(clustersF{clust},fieldnames(Max_resp_f20_perfishNbrain2.(RegionList{brain})))
 for f=1:size(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust}),1)     
[R,P,RL,RU] =corrcoef(Max_resp_f20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(f,2:31),Behaviour(:,1));
temp(1,1)=R(1,2);temp(1,2)=P(1,2);temp(1,3)=RL(1,2);temp(1,4)=RU(1,2);
Correlations_BrainVSbehav.f20.(RegionList{brain}).(clustersF{clust})(f,:)=temp;
 clear R P RL RU
 end
else
    0 %%%% the zeros are just to flag that there are some cases where I wont have data.   
end 

if ismember(clustersF{clust},fieldnames(Max_resp_f60_perfishNbrain2.(RegionList{brain})))
for f=1:size(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust}),1)     
[R,P,RL,RU] =corrcoef(Max_resp_f60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(f,2:31),Behaviour(:,1));
temp(1,1)=R(1,2);temp(1,2)=P(1,2);temp(1,3)=RL(1,2);temp(1,4)=RU(1,2);
Correlations_BrainVSbehav.f60.(RegionList{brain}).(clustersF{clust})(f,:)=temp;
 clear R P RL RU
end
else
    0
end

if ismember(clustersF{clust},fieldnames(Max_resp_s20_perfishNbrain2.(RegionList{brain})))
for f=1:size(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust}),1)     
[R,P,RL,RU] =corrcoef(Max_resp_s20_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(f,2:31),Behaviour(:,1));
temp(1,1)=R(1,2);temp(1,2)=P(1,2);temp(1,3)=RL(1,2);temp(1,4)=RU(1,2);
Correlations_BrainVSbehav.s20.(RegionList{brain}).(clustersF{clust})(f,:)=temp;
 clear R P RL RU
end
else
    0
end

if ismember(clustersF{clust},fieldnames(Max_resp_s60_perfishNbrain2.(RegionList{brain})))
for f=1:size(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust}),1)     
[R,P,RL,RU] =corrcoef(Max_resp_s60_perfishNbrain2.(RegionList{brain}).(clustersF{clust})(f,2:31),Behaviour(:,1));
temp(1,1)=R(1,2);temp(1,2)=P(1,2);temp(1,3)=RL(1,2);temp(1,4)=RU(1,2);
Correlations_BrainVSbehav.s60.(RegionList{brain}).(clustersF{clust})(f,:)=temp;
 clear R P RL RU
end
else
    0
end

counter=counter+1;

end

end
clear temp

%%%% visualizing the results
%%% I am asking that there are at least 3 fish
Corr_BrainVSbehav_mat=struct;
Corr_mat_mean=[];
for data=1:length(datasets)
    
    for brain=1:length(RegionList)
        counter=1;
        for clust=[2 3 1] %%% to put the clusters in order (fasthab, slopehab, nonhab).
            if ismember(clustersF{clust},fieldnames(Correlations_BrainVSbehav.(datasets(data,:)).(RegionList{brain})))
            if 3<sum(~isnan(Correlations_BrainVSbehav.(datasets(data,:)).(RegionList{brain}).(clustersF{clust})(:,1)));
            temp=nanmean(Correlations_BrainVSbehav.(datasets(data,:)).(RegionList{brain}).(clustersF{clust})(:,1));
            else
            temp=NaN;
            end
            else
            temp=NaN;
            end
            Corr_BrainVSbehav_mat.(datasets(data,:))(counter,brain)=temp;
            counter=counter+1;
        end
    end
    
    Corr_mat_mean=cat(3,Corr_mat_mean,Corr_BrainVSbehav_mat.(datasets(data,:)));
    Corr_mat_mean=mean(Corr_mat_mean,3); %%%% note, as I am doing a normal mean I am leaving out the combinations that have NaNs. so this means that for it to be included there most be at least 3 fish in the 4 datasets. 
    
    figure;imagesc(Corr_BrainVSbehav_mat.(datasets(data,:))),caxis([0 1]);xtickangle(45);xticklabels(RegionList);yticks([1 2 3]);yticklabels(clustersF([2 3 1]));
end
clear temp
figure;imagesc(Corr_mat_mean),caxis([0 1]);xtickangle(45);xticklabels(RegionList);yticks([1 2 3]);yticklabels(clustersF([2 3 1]));

%%%% now with inferno colormap. need to add to path 
figure;imagesc(Corr_mat_mean),caxis([0 1]);xtickangle(45);xticklabels(RegionList);yticks([1 2 3]);yticklabels(clustersF([2 3 1]));colorbar;colormap(inferno);
saveas(gcf,'Max_resp_Brain_corr_freeswim.svg');

