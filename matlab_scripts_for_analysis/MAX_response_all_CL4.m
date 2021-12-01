


%%%% this script is to get the differences in loom responses across
%%%% datasets and cluster types. It will take the peak of each loom and compare it to the
%%%% first loom response in each ROI for all 4 datasets. Then it will get
%%%% the mean of each fish for each dataset. 


%%% first with f20


%%%% load what you need

load('final_F20_step1.mat','ZS_f20','idx_Fish_f20');
load('f20_cleaned_idxs.mat')
load('allfish_f20looms_tailmov.mat');

%%
 %%% this is to get the onset and the loom moments

%%% to take the loom onsets based on the movie

%%% for f20
Loomf20_onset=zeros(6720,1);
startpoints=[30,254,478]; %%in seconds
loom_times=[0,22,40,60,78,100,120,140,158,180]; %%in seconds
for p=1:3
for k=1:10
    Loomf20_onset(startpoints(p)*10+loom_times(k)*10)=2;
end
end

figure; plot(Loomf20_onset);
%%

Loomf20_onset_idx=find(Loomf20_onset==2);
Loomf20_onset_idx=round(Loomf20_onset_idx/5);
Loomf20_onset_idx(Loomf20_onset_idx<1)=[];
Loomf20_onset_idx(Loomf20_onset_idx>(length(allfish_f20looms_tailmov)-1))=[];        
  
%%
%%% this is to try to get the idx to normalize with each loom. 


loom_moments={};
for i=1:30
    if i==1||i==11%||i==21
    loom_moments{i}=(Loomf20_onset_idx(i)-59):Loomf20_onset_idx(i+1);
    %elseif i==21
    % loom_moments{i}=(Loomf20_onset_idx(i)-29):Loomf20_onset_idx(i+1); %%% i added this to not take into acount the sound response   
    elseif i==10||i==20||i==30
      loom_moments{i}=Loomf20_onset_idx(i)+1:Loomf20_onset_idx(i)+28;  
    else
     loom_moments{i}=Loomf20_onset_idx(i)+1:Loomf20_onset_idx(i+1);   
    end
end

%%% to check
for i=1:30
temp(1,i)=length(loom_moments{i});

end
sum(temp(1,:)) %%% i took away the 30s before the 21st loom to not pick up sound responses. now it should add 1344-60=1284


%%

%%% this is to get the means of each cluster per fish. 


clust_f20_CL4_cleaned2=struct2cell(clust_f20_CL4_cleaned);

Means_CL4_per_fish_f20={};
fish=unique(idx_Fish_f20);
for clust=1:length(clust_f20_CL4_cleaned2)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f20==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_f20_CL4_cleaned2{clust,1});
temp_mean=mean(ZS_f20(temp_idx,:));

Means_CL4_per_fish_f20{clust}(tempfish,:)=temp_mean;
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_f20_CL4_cleaned2)
figure;
for i=1:length(fish)
   plot(Means_CL4_per_fish_f20{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Means_CL4_per_fish_f20{1,clust}))

end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f20_perfish={};

for clust=1:length(clust_f20_CL4_cleaned2)

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL4_per_fish_f20{1,clust}(i,loom_moments{1,k}))-min(Means_CL4_per_fish_f20{1,clust}(i,Loomf20_onset_idx(k)+1));
end


Max_resp_f20_perfish{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_f20_CL4_cleaned2)
figure;
for i=1:length(fish)
   plot(Max_resp_f20_perfish{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Max_resp_f20_perfish{1,clust}))

end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Means_CL4_per_fish_f20 Max_resp_f20_perfish

%%

%%% now for f60



load('final_F60_step1_2.mat','ZS_f60','idx_Fish_f60','ZS_short_F60');
load('f60_cleaned_idxs.mat')


%%

%%% this is to get the means of each cluster per fish. 


clust_f60_CL4_cleaned2=struct2cell(clust_f60_CL4_cleaned);

Means_CL4_per_fish_f60={};
fish=unique(idx_Fish_f60);

fish(find(fish==47))=[]; %%% cause I also took out fish 47...

for clust=1:length(clust_f60_CL4_cleaned2)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_f60==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_f60_CL4_cleaned2{clust,1});
temp_mean=mean(ZS_f60(temp_idx,ZS_short_F60));

Means_CL4_per_fish_f60{clust}(tempfish,:)=temp_mean;
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_f60_CL4_cleaned2)
figure;
for i=1:length(fish)
   plot(Means_CL4_per_fish_f60{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Means_CL4_per_fish_f60{1,clust}))

end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_f60_perfish={};

for clust=1:length(clust_f60_CL4_cleaned2)

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL4_per_fish_f60{1,clust}(i,loom_moments{1,k}))-min(Means_CL4_per_fish_f60{1,clust}(i,Loomf20_onset_idx(k)+1));
end


Max_resp_f60_perfish{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_f60_CL4_cleaned2)
figure;
for i=1:length(fish)
   plot(Max_resp_f60_perfish{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Max_resp_f60_perfish{1,clust}))

end


clearvars -except Loomf20_onset Loomf20_onset_idx loom_moments Means_CL4_per_fish_f20 Max_resp_f20_perfish Means_CL4_per_fish_f60 Max_resp_f60_perfish




%%

%%% to trim the slow loom movies

S_trim=[1:448 453:901 906:1352];

%%


%%% now for s20

load('final_S20_step1.mat','ZS_s20','idx_Fish_s20','ZS_short_S60');
load('s20_cleaned_idxs.mat')


%%

%%% this is to get the means of each cluster per fish. 


clust_s20_CL4_cleaned2=struct2cell(clust_s20_CL4_cleaned);

Means_CL4_per_fish_s20={};
fish=unique(idx_Fish_s20);


for clust=1:length(clust_s20_CL4_cleaned2)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s20==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_s20_CL4_cleaned2{clust,1});
temp_mean=mean(ZS_s20(temp_idx,S_trim));

Means_CL4_per_fish_s20{clust}(tempfish,:)=temp_mean;
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_s20_CL4_cleaned2)
figure;
for i=1:length(fish)
   plot(Means_CL4_per_fish_s20{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Means_CL4_per_fish_s20{1,clust}))

end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_s20_perfish={};

for clust=1:length(clust_s20_CL4_cleaned2)

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL4_per_fish_s20{1,clust}(i,loom_moments{1,k}))-min(Means_CL4_per_fish_s20{1,clust}(i,Loomf20_onset_idx(k)+1));
end


Max_resp_s20_perfish{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_s20_CL4_cleaned2)
figure;
for i=1:length(fish)
   plot(Max_resp_s20_perfish{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Max_resp_s20_perfish{1,clust}))

end

clearvars -except Loomf20_onset loom_moments Loomf20_onset_idx S_trim Means_CL4_per_fish_f20 Max_resp_f20_perfish Means_CL4_per_fish_f60 Max_resp_f60_perfish Means_CL4_per_fish_s20 Max_resp_s20_perfish

%%


%%% now for s60

load('final_S60_step1.mat','ZS_s60','idx_Fish_s60','ZS_short_S60');
load('s60_cleaned_idxs.mat')


%%

%%% this is to get the means of each cluster per fish. 


clust_s60_CL4_cleaned2=struct2cell(clust_s60_CL4_cleaned);

Means_CL4_per_fish_s60={};
fish=unique(idx_Fish_s60);


for clust=1:length(clust_s60_CL4_cleaned2)
    
for tempfish=1:length(fish)
   
temp_idx_fish=find(idx_Fish_s60==fish(tempfish));
temp_idx=intersect(temp_idx_fish,clust_s60_CL4_cleaned2{clust,1});
temp_mean=mean(ZS_s60(temp_idx,ZS_short_S60(S_trim)));

Means_CL4_per_fish_s60{clust}(tempfish,:)=temp_mean;
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_s60_CL4_cleaned2)
figure;
for i=1:length(fish)
   plot(Means_CL4_per_fish_s60{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Means_CL4_per_fish_s60{1,clust}))

end

%%
%%% now lets look at the difference in max resposne. first make a matrix of
%%% max responses per loom, for each cluster per fish. 

Max_resp_s60_perfish={};

for clust=1:length(clust_s60_CL4_cleaned2)

for  i=1:length(fish)  
for k=1:30

   temp_Max_resp(1,k)= max(Means_CL4_per_fish_s60{1,clust}(i,loom_moments{1,k}))-min(Means_CL4_per_fish_s60{1,clust}(i,Loomf20_onset_idx(k)+1));
end


Max_resp_s60_perfish{clust}(i,:)=temp_Max_resp(1,:)/temp_Max_resp(1,1);
end
end

%%% to check it worked
%%% the graphs are also interesting to see the individual variability. 

for clust=1:length(clust_s60_CL4_cleaned2)
figure;
for i=1:length(fish)
   plot(Max_resp_s60_perfish{1,clust}(i,:)); 
    hold on;
end

figure;plot(mean(Max_resp_s60_perfish{1,clust}))

end


clearvars -except Loomf20_onset loom_moments Loomf20_onset_idx S_trim Means_CL4_per_fish_f20 Max_resp_f20_perfish Means_CL4_per_fish_f60 Max_resp_f60_perfish Means_CL4_per_fish_s20 Max_resp_s20_perfish Means_CL4_per_fish_s60 Max_resp_s60_perfish
%%


%%% this is to compare the datasets


datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20');
clustersF=fieldnames(mean_CL4_f20);
clustersS=fieldnames(mean_CL4_s20);

%%

%%% organizing the data
All_means_maxResp_CL4_normalized={};


counter=1;
figure;
for data=1:size(datasets,1)
  datasetname=strcat('Max_resp_',datasets(data,:),'_perfish');
  tempVar=eval(datasetname);

if ~isempty(regexp(datasets(data,:),'f'))
      
    for j=1:size(clustersF,1)


    All_means_maxResp_CL4_normalized.(datasets(data,:)).(clustersF{j,1})=nanmean(tempVar{1,j});

    subplot(4,7,counter);plot(All_means_maxResp_CL4_normalized.(datasets(data,:)).(clustersF{j,1}));

    counter=counter+1
    end

else

    for j=1:size(clustersS,1)


    All_means_maxResp_CL4_normalized.(datasets(data,:)).(clustersS{j,1})=nanmean(tempVar{1,j});

    subplot(4,7,counter);plot(All_means_maxResp_CL4_normalized.(datasets(data,:)).(clustersS{j,1}));

    counter=counter+1

 end

    
end

end


clear tempVar  datasetname

All_means_maxResp_CL4_normalized=struct2cell(All_means_maxResp_CL4_normalized);


for j=1:4
figure;
for i=1:4
   plot(All_means_maxResp_CL4_normalized{i,1}.(clustersF{j,1}))

hold on;
end 



end
        
        
 clearvars -except Loomf20_onset loom_moments Loomf20_onset_idx S_trim Means_CL4_per_fish_f20 Max_resp_f20_perfish Means_CL4_per_fish_f60 Max_resp_f60_perfish Means_CL4_per_fish_s20 Max_resp_s20_perfish Means_CL4_per_fish_s60 Max_resp_s60_perfish All_means_maxResp_CL4_normalized
       
save('All_means_maxResp_CL4_normalized.mat');


