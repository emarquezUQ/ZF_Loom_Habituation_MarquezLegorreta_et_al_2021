

%%% this script is to look for the multisensory respondent neurons. 


%%% for fish s20

%%%% first you need to load what you need

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab

load('BrainReg_S20.mat')
load('final_S20_step1.mat','MatFiles_s20','ZS_s20','idx_Fish_s20','idx_Plane_s20','rawregressS20','idx_rsq_test_s20short','High_corr_Nb_s20','High_corr_Nb_s20_short');
load('Zbrain_Masks.mat');
load('final_S20_step1.mat','gooodmaps');

%%%% this is to look for the multisense regressor

%%% building the regressor

%idx_temp1=find(High_corr_Nb_s20_short==2);

multisens_regress=zeros(3,length(rawregressS20));

gooodmaps

for i=1:3
   idx_temp1=find(High_corr_Nb_s20_short==gooodmaps(i));
multisens_regress(i,:)=mean(ZS_s20(idx_rsq_test_s20short(idx_temp1),:));
multisens_regress(i,897:956)=rawregressS20(6,897:956);

end


figure;
for i=1:3
subplot(1,3,i)
    plot(multisens_regress(i,:));

end

%%

ModelResults_shortS20_all_multisense={};
for j=1:size(multisens_regress,1)
    
    
ModelResults_shortS20=[];
parfor i=1:size(ZS_s20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(multisens_regress(j,:)',ZS_s20(i,:));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS20(i).coef=mdl.Coefficients;
    ModelResults_shortS20(i).MSE=mdl.MSE;
    ModelResults_shortS20(i).Fitted=mdl.Fitted;
    ModelResults_shortS20(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortS20_all_multisense{j}=ModelResults_shortS20;

%rsquare_loom_shortS20=[ModelResults_shortS20.rsquared];


% figure;histogram(rsquare_loom_shortS20,'Normalization','probability');


%idx_rsq_test_s20short=find(rsquare_loom_shortS20>0.3 & rsquare_loom_shortS20<1); %%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
%proportion_s20=length(idx_rsq_test_s20short)/length(rsquare_loom_shortS20);

% figure;imagesc(ZS_s20(idx_rsq_test_s20short,:), [-0.5 4]);colormap hot

end


%%
idx_rsq_test_s20short_all_multisense={};
for j=1:size(multisens_regress)
    temp_rsq=[ModelResults_shortS20_all_multisense{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
        
    
idx_rsq_test_s20short_all_multisense{j}=idx_rsq_temp;

end

idx_rsq_test_s20short_multisense=horzcat(idx_rsq_test_s20short_all_multisense{:});
idx_rsq_test_s20short_multisense=unique(idx_rsq_test_s20short_multisense);


figure;imagesc(ZS_s20(idx_rsq_test_s20short_multisense,:), [-0.5 4]);colormap hot


%%

%%% i still get some that are not responding to both stimuli... so i will
%%% finish up the cleaning with a correlation to the sound response.


Correlation_multisense=[];

temp_ZS=ZS_s20(idx_rsq_test_s20short_multisense,:);
    for idx=1:size(idx_rsq_test_s20short_multisense,2)
        
            temp=corrcoef(rawregressS20(6,897:956),temp_ZS(idx,897:956));
            Correlation_multisense(1,idx)=temp(1,2);            
        
    end
    
clearvars temp_ZS 


idx_multisense=[];

Threshold=0.5;

idx_multisense=idx_rsq_test_s20short_multisense(find(Correlation_multisense>Threshold));
    

figure;imagesc(ZS_s20(idx_multisense,:), [-0.5 4]);colormap hot

figure;plot(mean(ZS_s20(idx_multisense,:)));

figure;scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2),'.');


figure;scatter3(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2),ROI_temp2(idx_multisense,3),'.');


% figure;
% scatter3(ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),1),ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),2),ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),3),'.'); %%% this is my rotation



%%



%%% to clasify them with a correlation

Correlation_group_multisense_s20={};
counter=1;
for i=1:size(multisens_regress,1)
    Correlation_temp=[];
    for idx=1:size(ZS_s20(idx_multisense),2)
        temp_corr=corrcoef(multisens_regress(i,:),ZS_s20(idx_multisense(idx),:));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_multisense_s20{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_multisense_mat_s20=[];
for n=1:size(multisens_regress,1)
Correlation_group_multisense_mat_s20(n,:)=Correlation_group_multisense_s20{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_multisense_Nb_s20=zeros(length(Correlation_group_multisense_mat_s20),1);
for i=1:length(Correlation_group_multisense_mat_s20)
    [~,I]=max(Correlation_group_multisense_mat_s20(:,i));
    High_corr_multisense_Nb_s20(i,1)=I;
    
end

%%

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(multisens_regress,1);counter=1;
for i=1:size(multisens_regress,1)
    
    idx_temp=find(High_corr_multisense_Nb_s20==i);
    subplot(rows,4,counter);plot(mean(ZS_s20(idx_multisense(idx_temp),:),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s20(idx_multisense(idx_temp),:),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s20(idx_multisense(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s20(idx_multisense(idx_temp)));%%% for the fish location
    counter=counter+4;
end


save('multisense_S20.mat','idx_multisense','idx_rsq_test_s20short_multisense','Correlation_multisense','ModelResults_shortS20_all_multisense','multisens_regress','Correlation_group_multisense_mat_s20','High_corr_multisense_Nb_s20','-v7.3');

