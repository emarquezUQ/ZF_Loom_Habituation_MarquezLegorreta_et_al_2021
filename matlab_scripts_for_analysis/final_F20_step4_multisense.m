%%% this script is to look for the multisensory respondent neurons. 


%%% for fish f20

%%%% first you need to load what you need

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab

load('BrainReg_F20.mat')
load('final_F20_step1.mat','MatFiles_f20','ZS_f20','idx_Fish_f20','idx_Plane_f20','rawregressF20','idx_rsq_test_f20short','High_corr_Nb_f20','High_corr_Nb_f20_short');
load('Zbrain_Masks.mat');
load('final_F20_step1.mat','gooodmaps');

%%%% this is to look for the multisense regressor

%%% building the regressor

idx_temp1=find(High_corr_Nb_f20_short==2);

multisens_regress=zeros(3,length(rawregressF20));

gooodmaps

for i=1:3
   idx_temp1=find(High_corr_Nb_f20_short==gooodmaps(i));
multisens_regress(i,:)=mean(ZS_f20(idx_rsq_test_f20short(idx_temp1),:));
multisens_regress(i,897:956)=rawregressF20(3,897:956);

end


figure;
for i=1:3
subplot(1,3,i)
    plot(multisens_regress(i,:));

end

%%

ModelResults_shortF20_all_multisense={};
for j=1:size(multisens_regress,1)
    
    
ModelResults_shortF20=[];
parfor i=1:size(ZS_f20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(multisens_regress(j,:)',ZS_f20(i,:));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortF20(i).coef=mdl.Coefficients;
    ModelResults_shortF20(i).MSE=mdl.MSE;
    ModelResults_shortF20(i).Fitted=mdl.Fitted;
    ModelResults_shortF20(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortF20_all_multisense{j}=ModelResults_shortF20;

%rsquare_loom_shortF20=[ModelResults_shortF20.rsquared];


% figure;histogram(rsquare_loom_shortF20,'Normalization','probability');


%idx_rsq_test_f20short=find(rsquare_loom_shortF20>0.3 & rsquare_loom_shortF20<1); %%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
%proportion_f20=length(idx_rsq_test_f20short)/length(rsquare_loom_shortF20);

% figure;imagesc(ZS_f20(idx_rsq_test_f20short,:), [-0.5 4]);colormap hot

end


%%
idx_rsq_test_f20short_all_multisense={};
for j=1:size(multisens_regress)
    temp_rsq=[ModelResults_shortF20_all_multisense{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
        
    
idx_rsq_test_f20short_all_multisense{j}=idx_rsq_temp;

end

idx_rsq_test_f20short_multisense=horzcat(idx_rsq_test_f20short_all_multisense{:});
idx_rsq_test_f20short_multisense=unique(idx_rsq_test_f20short_multisense);


figure;imagesc(ZS_f20(idx_rsq_test_f20short_multisense,:), [-0.5 4]);colormap hot


%%

%%% i still get some that are not responding to both stimuli... so i will
%%% finish up the cleaning with a correlation to the sound response.


Correlation_multisense=[];

temp_ZS=ZS_f20(idx_rsq_test_f20short_multisense,:);
    for idx=1:size(idx_rsq_test_f20short_multisense,2)
        
            temp=corrcoef(rawregressF20(3,897:956),temp_ZS(idx,897:956));
            Correlation_multisense(1,idx)=temp(1,2);            
        
    end
    
clearvars temp_ZS 


idx_multisense=[];

Threshold=0.5;

idx_multisense=idx_rsq_test_f20short_multisense(find(Correlation_multisense>Threshold));
    

figure;imagesc(ZS_f20(idx_multisense,:), [-0.5 4]);colormap hot

figure;plot(mean(ZS_f20(idx_multisense,:)));

scatter3(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2),ROI_temp2(idx_multisense,3),'.');


figure;
scatter3(ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),1),ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),2),ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),3),'.'); %%% this is my rotation



%%



%%% to clasify them with a correlation

Correlation_group_multisense_f20={};
counter=1;
for i=1:size(multisens_regress,1)
    Correlation_temp=[];
    for idx=1:size(ZS_f20(idx_multisense),2)
        temp_corr=corrcoef(multisens_regress(i,:),ZS_f20(idx_multisense(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_multisense_f20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_multisense_mat_f20=[];
for n=1:size(multisens_regress,1)
Correlation_group_multisense_mat_f20(n,:)=Correlation_group_multisense_f20{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_multisense_Nb_f20=zeros(length(Correlation_group_multisense_mat_f20),1);
for i=1:length(Correlation_group_multisense_mat_f20)
    [~,I]=max(Correlation_group_multisense_mat_f20(:,i));
    High_corr_multisense_Nb_f20(i,1)=I;
    
end

%%

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(multisens_regress,1);counter=1;
for i=1:size(multisens_regress,1)
    
    idx_temp=find(High_corr_multisense_Nb_f20==i);
    subplot(rows,4,counter);plot(mean(ZS_f20(idx_multisense(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_f20(idx_multisense(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_f20(idx_multisense(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_f20(idx_multisense(idx_temp)));%%% for the fish location
    counter=counter+4;
end


save('multisense_F20.mat','idx_multisense','idx_rsq_test_f20short_multisense','Correlation_multisense','ModelResults_shortF20_all_multisense','multisens_regress','Correlation_group_multisense_mat_f20','High_corr_multisense_Nb_f20','-v7.3');

