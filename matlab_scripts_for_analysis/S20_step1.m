

%%%% s20 final_step1

%%%%% this script is to get the zscored calcium traces of the s20 group, 
%%%%% perform the linear regression to each of the regressors and clasify
%%%%% them based on highest correlation. 


%%% for s20

%%% loading the ZS with the noise, matfiles,regressors, planes and fish. 
load('inhib_s20_regress_200CL.mat','inhib_s20_regress');

S20data_CN=load('s20_postKmeans_CN.mat','ZS_CN','MatFiles');
ZS_s20=S20data_CN.('ZS_CN');
MatFiles_s20=S20data_CN.('MatFiles');

S20data=load('s20_r2050_CL6.mat','idx_Plane','idx_Fish');

idx_Plane_s20=S20data.('idx_Plane');
idx_Fish_s20=S20data.('idx_Fish');


rawregressS20=load('rawregressS20.mat','rawregress');
rawregressS20 = rawregressS20.('rawregress');


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


ModelResults_shortS20_all={};
for j=1:length(rawregressS20)
    
    
ModelResults_shortS20=[];
parfor i=1:size(ZS_s20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(rawregressS20(j,:)',ZS_s20(i,:));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS20(i).coef=mdl.Coefficients;
    ModelResults_shortS20(i).MSE=mdl.MSE;
    ModelResults_shortS20(i).Fitted=mdl.Fitted;
    ModelResults_shortS20(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortS20_all{j}=ModelResults_shortS20;


end




%%

idx_rsq_test_s20short_all={};
for j=1:size(rawregressS20)
    temp_rsq=[ModelResults_shortS20_all{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
        
    
idx_rsq_test_s20short_all{j}=idx_rsq_temp;

end

idx_rsq_test_s20short=horzcat(idx_rsq_test_s20short_all{:});
idx_rsq_test_s20short=unique(idx_rsq_test_s20short);


figure;imagesc(ZS_s20(idx_rsq_test_s20short,:), [-0.5 4]);colormap hot


%%

%%% to clasify them with a correlation

Correlation_group_s20={};
counter=1;
for i=1:size(rawregressS20,1)
    Correlation_temp=[];
    for idx=1:size(ZS_s20(idx_rsq_test_s20short),2)
        temp_corr=corrcoef(rawregressS20(i,:),ZS_s20(idx_rsq_test_s20short(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_s20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_mat_s20=[];
for n=1:size(rawregressS20,1)
Correlation_group_mat_s20(n,:)=Correlation_group_s20{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_Nb_s20=zeros(length(Correlation_group_mat_s20),1);
for i=1:length(Correlation_group_mat_s20)
    [~,I]=max(Correlation_group_mat_s20(:,i));
    High_corr_Nb_s20(i,1)=I;
    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(rawregressS20,1);counter=1;
for i=1:size(rawregressS20,1)
    
    idx_temp=find(High_corr_Nb_s20==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_rsq_test_s20short(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_rsq_test_s20short(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_rsq_test_s20short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_rsq_test_s20short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL7_S20_CN_each','-dpdf','-bestfit');

%%

High_corr_Nb_s20_short=High_corr_Nb_s20;
High_corr_Nb_s20_short(find(High_corr_Nb_s20==2))=1;
High_corr_Nb_s20_short(find(High_corr_Nb_s20==4))=1;

gooodmaps=[1 3 5 7];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(gooodmaps,2);counter=1;
for i=gooodmaps
    
    idx_temp=find(High_corr_Nb_s20_short==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_rsq_test_s20short(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_rsq_test_s20short(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_rsq_test_s20short(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_rsq_test_s20short(idx_temp)));%%% for the fish location
    counter=counter+4;
end

%print(gcf,'multigraph_rsq030_CL4_S20_CN_each','-dpdf','-bestfit');

%%

 save('final_S20_step1.mat','-v7.3');


