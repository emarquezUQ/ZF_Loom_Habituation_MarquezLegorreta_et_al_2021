
%%%%% this script is to get the zscored calcium traces of the f20 group, 
%%%%% perform the linear regression to each of the regressors and clasify
%%%%% them based on highest correlation. 


%%%% f20 step1

%%% Getting the ZS with the noise, matfiles,regressors, planes and fish. 


F20data_CN=load('f20_CN_r2050_CL5_extra.mat','ZS_CN','MatFiles');
ZS_f20=F20data_CN.('ZS_CN');
MatFiles_f20=F20data_CN.('MatFiles');


F20data=load('f20_r2050_CL5.mat','idx_Plane','idx_Fish');

idx_Plane_f20=F20data.('idx_Plane');
idx_Fish_f20=F20data.('idx_Fish');


rawregressF20=load('rawregressF20.mat','rawregress');
rawregressF20 = rawregressF20.('rawregress');

rawregressF20_CN=load('f20_CN_r2050_CL5_extra.mat','rawregress_CN');
rawregressF20_CN = rawregressF20_CN.('rawregress_CN');


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;rows=length(rawregressF20);

for i=1:length(rawregressF20)%% 
   
    subplot(4,2,counter);plot(rawregressF20(i,:));
     counter=counter+1;
end

rawregressF20(7,:)=rawregressF20_CN(1,:);  %% the inhibition cluster

%%


%%% this is for a linear regression using the f20 denoised data regressors
%%% from the 50CL + the inhibition cluster

ModelResults_shortF20_all={};
for j=1:length(rawregressF20)
    
    
ModelResults_shortF20=[];
parfor i=1:size(ZS_f20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressF20(j,:)',ZS_f20(i,:));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortF20(i).coef=mdl.Coefficients;
    ModelResults_shortF20(i).MSE=mdl.MSE;
    ModelResults_shortF20(i).Fitted=mdl.Fitted;
    ModelResults_shortF20(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortF20_all{j}=ModelResults_shortF20;

end


%%

idx_rsq_test_f20short_all={};
for j=1:size(rawregressF20)
    temp_rsq=[ModelResults_shortF20_all{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
        
    
idx_rsq_test_f20short_all{j}=idx_rsq_temp;

end

idx_rsq_test_f20short=horzcat(idx_rsq_test_f20short_all{:});
idx_rsq_test_f20short=unique(idx_rsq_test_f20short);


figure;imagesc(ZS_f20(idx_rsq_test_f20short,:), [-0.5 4]);colormap hot


%%

%%% to clasify them with a correlation

Correlation_group_f20={};
counter=1;
for i=1:size(rawregressF20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_f20(idx_rsq_test_f20short),2)
        temp_corr=corrcoef(rawregressF20(i,:),ZS_f20(idx_rsq_test_f20short(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_f20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_f20=[];
for n=1:size(rawregressF20,1)
Correlation_group_mat_f20(n,:)=Correlation_group_f20{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_f20=zeros(length(Correlation_group_mat_f20),1);
for i=1:length(Correlation_group_mat_f20)
    [~,I]=max(Correlation_group_mat_f20(:,i));
    High_corr_Nb_f20(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressF20,1);counter=1;
for i=1:size(rawregressF20,1)
    
    idx_temp=find(High_corr_Nb_f20==i);
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_rsq_test_f20short(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_rsq_test_f20short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_rsq_test_f20short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_F20_CN_each','-dpdf','-bestfit');

%%

High_corr_Nb_f20_short=High_corr_Nb_f20;
High_corr_Nb_f20_short(find(High_corr_Nb_f20==4))=2;  %%% merging strong hab
High_corr_Nb_f20_short(find(High_corr_Nb_f20==5))=2;  %%% merging strong hab

gooodmaps=[1 2 6 7];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(gooodmaps,2);counter=1;
for i=gooodmaps
    
    idx_temp=find(High_corr_Nb_f20_short==i);
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_rsq_test_f20short(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_rsq_test_f20short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_rsq_test_f20short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL4_F20_CN_each','-dpdf','-bestfit');

%%

 save('final_F20_step1.mat','-v7.3');






