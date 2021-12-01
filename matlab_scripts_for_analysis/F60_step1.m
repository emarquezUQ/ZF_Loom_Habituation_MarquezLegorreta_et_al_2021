


%%%% f60 final_step1

%%%%% this script is to get the zscored calcium traces of the f60 group, 
%%%%% perform the linear regression to each of the regressors and clasify
%%%%% them based on highest correlation. 


%%% for f60



F60data_CN=load('f60_postKmeans_CN_2.mat','ZS_CN','MatFiles');
ZS_f60=F60data_CN.('ZS_CN');
MatFiles_f60=F60data_CN.('MatFiles');


F60data=load('f60_r2050_CL4.mat','idx_Plane','idx_Fish');

idx_Plane_f60=F60data.('idx_Plane');
idx_Fish_f60=F60data.('idx_Fish');


%%% i am taking fish 37 away cause it itself has half the neurons (>10K) of the slope cluster.
ZS_f60(idx_Fish_f60==37,:)=[];

idx_Fish_allf60=idx_Fish_f60;

idx_Fish_f60(idx_Fish_allf60==37,:)=[];

idx_Plane_f60(idx_Fish_allf60==37,:)=[];



%%

rawregressF20=load('rawregressF20.mat','rawregress');
rawregressF20 = rawregressF20.('rawregress');

rawregressF20_CN=load('f20_CN_r2050_CL5_extra.mat','rawregress_CN');
rawregressF20_CN = rawregressF20_CN.('rawregress_CN');

rawregressF60=load('f60_r2050_CL4.mat','rawregress');
rawregressF60 = rawregressF60.('rawregress');



Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(rawregressF20);

for i=1:length(rawregressF20)%% +2 %%%for the multisensregress
   
    subplot(4,2,counter);plot(rawregressF20(i,:));
     counter=counter+1;
end

rawregressF20(7,:)=rawregressF20_CN(1,:);


%%% to take the timepoints of f60 to fit it in f20
 ZS_short_F60=zeros(size(rawregressF60(1,:)));
 startpoints=[0,584,1168]; %%in seconds

 loom_times=[0,96,150,210,264,330,390,450,504,570]; %%in seconds
 loom_length=[52,18,20,18,22,20,20,18,22,14]; %%% in seconds

for p=1:3
    for k=1:10
    ZS_short_F60(startpoints(p)*2+loom_times(k)*2+1:startpoints(p)*2+loom_times(k)*2+loom_length(k)*2)=1;
    end
end


ZS_short_F60=find(ZS_short_F60==1);

%figure; plot(rawregressF60(1,ZS_short_F60)); %%% to check

%%

ModelResults_shortF60_all={};
for j=1:length(rawregressF20)
    
    
ModelResults_shortF60=[];
parfor i=1:size(ZS_f60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressF20(j,:)',ZS_f60(i,ZS_short_F60));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortF60(i).coef=mdl.Coefficients;
    ModelResults_shortF60(i).MSE=mdl.MSE;
    ModelResults_shortF60(i).Fitted=mdl.Fitted;
    ModelResults_shortF60(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortF60_all{j}=ModelResults_shortF60;


end




%%

idx_rsq_test_f60short_all={};
for j=1:size(rawregressF20)
    temp_rsq=[ModelResults_shortF60_all{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
        
    
idx_rsq_test_f60short_all{j}=idx_rsq_temp;

end

idx_rsq_test_f60short=horzcat(idx_rsq_test_f60short_all{:});
idx_rsq_test_f60short=unique(idx_rsq_test_f60short);


figure;imagesc(ZS_f60(idx_rsq_test_f60short,ZS_short_F60), [-0.5 4]);colormap hot


%%

%%% to clasify them with a correlation

Correlation_group_f60={};
counter=1;
for i=1:size(rawregressF20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_f60(idx_rsq_test_f60short),2)
        temp_corr=corrcoef(rawregressF20(i,:),ZS_f60(idx_rsq_test_f60short(idx),ZS_short_F60));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_f60{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_f60=[];
for n=1:size(rawregressF20,1)
Correlation_group_mat_f60(n,:)=Correlation_group_f60{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_f60=zeros(length(Correlation_group_mat_f60),1);
for i=1:length(Correlation_group_mat_f60)
    [~,I]=max(Correlation_group_mat_f60(:,i));
    High_corr_Nb_f60(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressF20,1);counter=1;
for i=1:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_F60_CN_each_2','-dpdf','-bestfit');


save('final_F60_step1_2.mat','-v7.3');


%%% some ROIs in the inhibitory cluster that
%%% are from fish 3 have some movement artifacts
%%% it is better to clean them

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=1;counter=1;
for i=7%:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short(idx_temp)));%%% for the fish location
    counter=counter+4;
end


idx_temp=idx_rsq_test_f60short(find(High_corr_Nb_f60==7));
idx_fish3_inh=intersect(find(idx_Fish_f60==3),idx_temp);
%figure;imagesc(ZS_f60(idx_fish3_inhb,ZS_short_F60),[0 3]);
idx_mov_fish3_plane1=intersect(idx_fish3_inh,find(idx_Plane==1));
%figure;imagesc(ZS_f60(idx_mov_fish3_plane1,ZS_short_F60),[0 3]);
idx_todelete=find(ismember(idx_rsq_test_f60short,idx_mov_fish3_plane1));



idx_rsq_test_f60short2=idx_rsq_test_f60short;
High_corr_Nb_f60_2=High_corr_Nb_f60;


idx_rsq_test_f60short2(idx_todelete)=[];
High_corr_Nb_f60_2(idx_todelete)=[];


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressF20,1);counter=1;
for i=1:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60_2==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short2(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short2(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_F60_CN_each_2_corrected','-dpdf','-bestfit');

%%

%%%% for CL4

High_corr_Nb_f60_short=High_corr_Nb_f60_2;
High_corr_Nb_f60_short(find(High_corr_Nb_f60_2==4))=2;
High_corr_Nb_f60_short(find(High_corr_Nb_f60_2==5))=2;

gooodmaps=[1 2 6 7];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(gooodmaps,2);counter=1;
for i=gooodmaps
    
    idx_temp=find(High_corr_Nb_f60_short==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short2(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short2(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short2(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL4_F60_CN_each_2','-dpdf','-bestfit');

%%

 save('final_F60_step1_2_short_correction.mat','idx_rsq_test_f60short2','High_corr_Nb_f60_2','idx_todelete','High_corr_Nb_f60_short','-v7.3');


%%% I also need to take away fish 47 as it has more than half of the
%%% responses for the fast hab broad cluster.
%%% i will do it through indexing to not repeat the whole analysis again.


load('final_F60_step1_2_short_correction.mat')

load('final_F60_step1_2.mat','MatFiles','ZS_f60','idx_Fish_f60','idx_Plane_f60','rawregressF20','ZS_short_F60');



%%% maybe not necesary
% ZS_f60(idx_Fish_f60==47,:)=[];
% 
% idx_Fish_allf60=idx_Fish_f60;
% 
% idx_Fish_f60(idx_Fish_allf60==47,:)=[];
% 
% idx_Plane_f60(idx_Fish_allf60==47,:)=[];



idx_todelete=find((idx_Fish_allf60(idx_rsq_test_f60short2)==47)>0);
idx_rsq_test_f60short3=idx_rsq_test_f60short2;
idx_rsq_test_f60short3(idx_todelete)=[];


High_corr_Nb_f60_3=High_corr_Nb_f60_2;

%idx_rsq_test_f60short3(idx_todelete)=[];
High_corr_Nb_f60_3(idx_todelete)=[];




Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressF20,1);counter=1;
for i=1:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f60_3==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short3(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short3(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short3(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short3(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_F60_CN_each_3_corrected','-dpdf','-bestfit');

High_corr_Nb_f60_short=High_corr_Nb_f60_3;
High_corr_Nb_f60_short(find(High_corr_Nb_f60_3==4))=2;
High_corr_Nb_f60_short(find(High_corr_Nb_f60_3==5))=2;

gooodmaps=[1 2 6 7];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(gooodmaps,2);counter=1;
for i=gooodmaps
    
    idx_temp=find(High_corr_Nb_f60_short==i);
    subplot(rows,4,counter);plot(mean(ZS_f60(idx_rsq_test_f60short3(idx_temp),ZS_short_F60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f60(idx_rsq_test_f60short3(idx_temp),ZS_short_F60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f60(idx_rsq_test_f60short3(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f60(idx_rsq_test_f60short3(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL4_F60_CN_each_3','-dpdf','-bestfit');

 save('final_F60_step1_3_short_correction.mat','idx_rsq_test_f60short3','High_corr_Nb_f60_3','idx_todelete','High_corr_Nb_f60_short','-v7.3');

