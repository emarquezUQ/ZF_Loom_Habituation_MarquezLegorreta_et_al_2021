

%%%% s60 final_step1

%%%%% this script is to get the zscored calcium traces of the s60 group, 
%%%%% perform the linear regression to each of the regressors and clasify
%%%%% them based on highest correlation. 

%%% for s60

%%% now the ZS with the noise, matfiles,regressors, planes and fish. 
load('inhib_s20_regress_200CL.mat','inhib_s20_regress');


S60data_CN=load('s60_postKmeans_CN.mat','ZS_CN','MatFiles');
ZS_s60=S60data_CN.('ZS_CN');
MatFiles_s60=S60data_CN.('MatFiles');

S60data=load('s60_r2050_CL3.mat','idx_Plane','idx_Fish');

idx_Plane_s60=S60data.('idx_Plane');
idx_Fish_s60=S60data.('idx_Fish');


rawregressS20=load('rawregressS20.mat','rawregress');
rawregressS20 = rawregressS20.('rawregress');

rawregressS60=load('rawregressS60.mat','rawregress');
rawregressS60 = rawregressS60.('rawregress');


%%

%%% to take the timepoints of s60 to fit it in s20
 ZS_short_S60=zeros(size(rawregressS60(1,:)));
 startpoints=[0,586,1172]; %%in seconds

 loom_times=[0,96,150,210,264,330,390,450,504,570]; %%in seconds
 loom_length=[52,18,20,18,22,20,20,18,22,16]; %%% in seconds

for p=1:3
    for k=1:10
    ZS_short_S60(startpoints(p)*2+loom_times(k)*2+1:startpoints(p)*2+loom_times(k)*2+loom_length(k)*2)=1;
    end
end


ZS_short_S60=find(ZS_short_S60==1);

%figure;plot(rawregressS60(1,ZS_short_S60)); %%% to check



%%


%%% this is for a linear regression using the s20 denoised data regressors
%%% from the 50CL + the inhibition cluster

% counter=1;
% figure;
% for i=1:length(rawregressS20)
% subplot(4,2,counter);plot(rawregressS20(i,:));
% counter=counter+1;
% end


rawregressS20(7,:)=inhib_s20_regress;

ModelResults_shortS60_all={};
for j=1:length(rawregressS20)
    
    
ModelResults_shortS60=[];
parfor i=1:size(ZS_s60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressS20(j,:)',ZS_s60(i,ZS_short_S60));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS60(i).coef=mdl.Coefficients;
    ModelResults_shortS60(i).MSE=mdl.MSE;
    ModelResults_shortS60(i).Fitted=mdl.Fitted;
    ModelResults_shortS60(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortS60_all{j}=ModelResults_shortS60;

rsquare_loom_shortS60=[ModelResults_shortS60.rsquared];


% figure;histogram(rsquare_loom_shortS60,'Normalization','probability');


idx_rsq_test_s60short=find(rsquare_loom_shortS60>0.3 & rsquare_loom_shortS60<1); %%% 
proportion_s60=length(idx_rsq_test_s60short)/length(rsquare_loom_shortS60);

% figure;imagesc(ZS_s60(idx_rsq_test_s60short,ZS_short_S60), [-0.5 4]);colormap hot

end




%%

idx_rsq_test_s60short_all={};
for j=1:size(rawregressS20)
    temp_rsq=[ModelResults_shortS60_all{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
       
    
idx_rsq_test_s60short_all{j}=idx_rsq_temp;

end

idx_rsq_test_s60short=horzcat(idx_rsq_test_s60short_all{:});
idx_rsq_test_s60short=unique(idx_rsq_test_s60short);


figure;imagesc(ZS_s60(idx_rsq_test_s60short,ZS_short_S60), [-0.5 4]);colormap hot


%%

%%% to clasify them with a correlation

Correlation_group_s60={};
counter=1;
for i=1:size(rawregressS20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_s60(idx_rsq_test_s60short),2)
        temp_corr=corrcoef(rawregressS20(i,:),ZS_s60(idx_rsq_test_s60short(idx),ZS_short_S60));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_s60{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_s60=[];
for n=1:size(rawregressS20,1)
Correlation_group_mat_s60(n,:)=Correlation_group_s60{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_s60=zeros(length(Correlation_group_mat_s60),1);
for i=1:length(Correlation_group_mat_s60)
    [~,I]=max(Correlation_group_mat_s60(:,i));
    High_corr_Nb_s60(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    
    idx_temp=find(High_corr_Nb_s60==i);
    subplot(rows,4,counter);plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s60(idx_rsq_test_s60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s60(idx_rsq_test_s60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_S60_CN_each','-dpdf','-bestfit');

%%

High_corr_Nb_s60_short=High_corr_Nb_s60;
High_corr_Nb_s60_short(find(High_corr_Nb_s60==2))=1;
High_corr_Nb_s60_short(find(High_corr_Nb_s60==4))=1;

gooodmaps=[1 3 5 7];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(gooodmaps,2);counter=1;
for i=gooodmaps
    
    idx_temp=find(High_corr_Nb_s60_short==i);
    subplot(rows,4,counter);plot(mean(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s60(idx_rsq_test_s60short(idx_temp),ZS_short_S60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s60(idx_rsq_test_s60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s60(idx_rsq_test_s60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL4_S60_CN_each','-dpdf','-bestfit');

%%

save('final_S60_step1.mat','-v7.3');

