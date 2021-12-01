


%%%%% this script is to do an analysis of the
%%%%% fmr1loomhab experiments. 

%%%%% Similarly to the wildtype loom habituation analysis, the results of the kmeans will be
%%%%% used as regressors to filter the ROIs and then they will be
%%%%% classified to the clusters using a correlation. 

%%%%% Note: as I am merging the ROIs fo all groups (ctrls, hets, mutants) I
%%%%% will also do a correlation to the regressors first to reduce the
%%%%% number of ROIs to test. 

%%%%% I will take some of the data alreay processed in the quick1stanalysis_fmr1loomhab.m script

%%% to load the previous data (e.g. ZS, Matfiles...)
load('s20_fmr1_loomhab_CN.mat','MatFiles','ZS_CN'); 

%clear GoodNoise %%% i dont need it anymore


cd C:\Emmanuel_temp\fmr1_loomhab\matlab_fmr1_loomhab

%%% to some more things
load('s20_fmr1_loomhab_CN_post_rsq01.mat','idx_rsq1');

load('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');

%%% to get the clusters from the Kmeans done to ALL the ROIs
load('s20_postKmeans_CN.mat','idxKmeans_ZS_CN','Cmap_ZS_CN','Model_ZS_CN1','GoodBetas_ZS_CN1');

load('s20_good_idx_Fish.mat','idx_Fish');
 
%%
%%% getting the regressors from the kmeans

GoodBetas_ZS_CN_selected=[13 21 29 32 35 36 38 40 43 44];

figure;
for i=1:length(GoodBetas_ZS_CN_selected)
subplot(4,3,i);plot(Cmap_ZS_CN(GoodBetas_ZS_CN_selected(i),:));
       
end

for i=1:length(GoodBetas_ZS_CN_selected)
GoodBetas_regress(i,:)=Cmap_ZS_CN(GoodBetas_ZS_CN_selected(i),:);

end

%%
%%% i will use this clusters to filter the data with a correlation of 0.5
%%% to reduce the number of ROIs to test. 

corrfilter=[];Threshold=0.5;
for i=1:length(GoodBetas_ZS_CN_selected)    
    
    corr_temp=zeros(1,length(ZS_CN));
    parfor jj=1:size(ZS_CN,1)
        temp_corr=corrcoef(GoodBetas_regress(i,:), ZS_CN(jj,:));
        corr_temp(jj)=temp_corr(1,2);
    end    
    %corrfilter(i).ZS=ZS_CN(find(corr_temp>Threshold),:);
    corrfilter(i).idx=find(corr_temp>Threshold);
    %corrfilter(i).mean=mean(corrfilter(i).ZS,1);
    %corrfilter(i).STD=std(corrfilter(i).ZS,1,1);       
end

idx_corr=[];
for i=1:length(GoodBetas_ZS_CN_selected)    

idx_corr=horzcat(idx_corr,corrfilter(i).idx);

end    
   
idx_corr=unique(idx_corr);

%%%% to check if it worked
figure;plot(mean(ZS_CN(corrfilter(1).idx,:)));

figure;
for i=1:length(GoodBetas_ZS_CN_selected)
subplot(4,3,i);plot(mean(ZS_CN(corrfilter(i).idx,:)));
       
end


%%

%%% now to do an individual linear regression with each regressor to the
%%% ROIs filtered after the correlation. 

ModelResults_shortS20_all={};
for j=1:length(GoodBetas_ZS_CN_selected)
    
    
ModelResults_shortS20=[];
parfor i=1:size(ZS_CN(idx_corr,:),1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(GoodBetas_regress(j,:)',ZS_CN(idx_corr(i),:));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS20(i).coef=mdl.Coefficients;
    ModelResults_shortS20(i).MSE=mdl.MSE;
    ModelResults_shortS20(i).Fitted=mdl.Fitted;
    ModelResults_shortS20(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortS20_all{j}=ModelResults_shortS20;

%rsquare_loom_shortS20=[ModelResults_shortS20.rsquared];


% figure;histogram(rsquare_loom_shortS20,'Normalization','probability');


%idx_rsq_test_s20short=find(rsquare_loom_shortS20>0.3 & rsquare_loom_shortS20<1); %%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
%proportion_s20=length(idx_rsq_test_s20short)/length(rsquare_loom_shortS20);

% figure;imagesc(ZS_s20(idx_rsq_test_s20short,:), [-0.5 4]);colormap hot

end

%% 

%%% this is to check the distribution of the rsq. but I have to keep in
%%% mind that i did a first filtering with a correlation so the
%%% disbribution will probably be different than in the original loom hab experiments

rsq_test_all=[];
for j=1:10
    temp_rsq=[ModelResults_shortS20_all{1,j}.rsquared];
      
rsq_test_all(:,j)=temp_rsq;

end


rsquare_loom=[];
for i=1:length(rsq_test_all)
    [M,I]=max(rsq_test_all(i,:));
    rsquare_loom(i,1)=M;
    
end

figure;histogram(rsquare_loom);
2*std(rsquare_loom) %% the result is 0.2867, actually very similar than with the other datasets.

save('rsquare_fmr1loomhab.mat','rsq_test_all','rsquare_loom');


%% filtering based on the rsq value

%%% using 0.3 as in the wildtype dataset
idx_rsq_test_s20short_all={};
for j=1:length(GoodBetas_ZS_CN_selected)
    temp_rsq=[ModelResults_shortS20_all{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
        
    
idx_rsq_test_s20short_all{j}=idx_rsq_temp;

end

idx_rsq_test_s20short=horzcat(idx_rsq_test_s20short_all{:});
idx_rsq_test_s20short=unique(idx_rsq_test_s20short);


%figure;imagesc(ZS_CN(idx_corr(idx_rsq_test_s20short),:), [-0.5 4]);colormap hot


%%% to put them in the idx format for ZS_CN

idx_rsq=idx_corr(idx_rsq_test_s20short);
%figure;imagesc(ZS_CN(idx_rsq,:), [-0.5 4]);colormap hot



%%

%%% to clasify them with a correlation

Correlation_group={};
counter=1;
for i=1:length(GoodBetas_ZS_CN_selected)
    Correlation_temp=[];
    for idx=1:size(ZS_CN(idx_rsq),2)
        temp_corr=corrcoef(GoodBetas_regress(i,:),ZS_CN(idx_rsq(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat=[];
for n=1:size(GoodBetas_ZS_CN_selected,2)
Correlation_group_mat(n,:)=Correlation_group{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb=zeros(length(Correlation_group_mat),1);
for i=1:length(Correlation_group_mat)
    [~,I]=max(Correlation_group_mat(:,i));
    High_corr_Nb(i,1)=I;
    
end

%%%%% NOTE: there was an error when making the Correlation_group_mat that I
%%%%% only noticed after saving the High_corr_Nb variable. a corrected
%%%%% version of it will be saved in another .mat file

save('s20_fmr1_loomhab_CN_part2_High_corr_Nb.mat','High_corr_Nb');
 

%%

save('s20_fmr1_loomhab_CN_part2_Model.mat','ModelResults_shortS20_all','idx_rsq_test_s20short_all','-v7.3');
 

idx_Fish_cat=categorical(idx_Fish);

save('s20_fmr1_loomhab_CN_part2.mat','GoodBetas_ZS_CN_selected','GoodBetas_regress','corrfilter','idx_corr','idx_rsq','Correlation_group','High_corr_Nb','idx_Fish_cat','-v7.3');
 

%%

%%% to visualize it


%%%% still need to work on it

% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1300, 900]);
% rows=size(GoodBetas_ZS_CN_selected,2)/2;counter=1;
% for i=1:size(GoodBetas_ZS_CN_selected,2)/2
%     
%     idx_temp=find(High_corr_Nb==i);
%     subplot(rows,4,counter);plot(mean(ZS_CN(idx_rsq(idx_temp),:),1)); %%%to plot the mean
%     subplot(rows,4,counter+1);imagesc(ZS_CN(idx_rsq(idx_temp),:),[0 3]);%%%for the raster plot
%     subplot(rows,4,counter+2);histogram(idx_Plane(idx_rsq(idx_temp))); %%%for the plane location   
%     subplot(rows,4,counter+3);histogram(idx_Fish_cat(idx_rsq(idx_temp)));%%% for the fish location
%     counter=counter+4;
% end



