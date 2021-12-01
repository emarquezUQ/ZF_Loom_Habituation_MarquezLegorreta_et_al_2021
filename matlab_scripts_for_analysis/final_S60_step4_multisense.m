


%%% this script is to look for the multisensory respondent neurons. 


%%% for fish s60

%%%% first you need to load what you need

cd /QRISdata/Q0291/Emmanuel_MeDiCi/FvsS_20vs60_CNMF/matlab

load('BrainReg_S60.mat')
load('final_S60_step1.mat','MatFiles_s60','ZS_s60','idx_Fish_s60','idx_Plane_s60','rawregressS60','idx_rsq_test_s60short','High_corr_Nb_s60','High_corr_Nb_s60_short');
load('Zbrain_Masks.mat');
load('final_S60_step1.mat','gooodmaps','ZS_short_S60');
load('final_S20_step1.mat','rawregressS20');


%%%% this is to look for the multisense regressor

%%% building the regressor

%idx_temp1=find(High_corr_Nb_s60_short==2);

multisens_regress=zeros(3,length(rawregressS20));

gooodmaps

for i=1:3
   idx_temp1=find(High_corr_Nb_s60_short==gooodmaps(i));
multisens_regress(i,:)=mean(ZS_s60(idx_rsq_test_s60short(idx_temp1),ZS_short_S60));
multisens_regress(i,897:956)=rawregressS20(6,897:956);

end


figure;
for i=1:3
subplot(1,3,i)
    plot(multisens_regress(i,:));

end

%%

ModelResults_shortS60_all_multisense={};
for j=1:size(multisens_regress,1)
    
    
ModelResults_shortS60=[];
parfor i=1:size(ZS_s60,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(multisens_regress(j,:)',ZS_s60(i,ZS_short_S60));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults_shortS60(i).coef=mdl.Coefficients;
    ModelResults_shortS60(i).MSE=mdl.MSE;
    ModelResults_shortS60(i).Fitted=mdl.Fitted;
    ModelResults_shortS60(i).rsquared=mdl.Rsquared.Adjusted;
     
end

ModelResults_shortS60_all_multisense{j}=ModelResults_shortS60;

%rsquare_loom_shortS60=[ModelResults_shortS60.rsquared];


% figure;histogram(rsquare_loom_shortS60,'Normalization','probability');


%idx_rsq_test_s60short=find(rsquare_loom_shortS60>0.3 & rsquare_loom_shortS60<1); %%% I am using 0.3 with the noised data as it gave similar results than 0.5 of the denoised data.
%proportion_s60=length(idx_rsq_test_s60short)/length(rsquare_loom_shortS60);

% figure;imagesc(ZS_s60(idx_rsq_test_s60short,:), [-0.5 4]);colormap hot

end


%%
idx_rsq_test_s60short_all_multisense={};
for j=1:size(multisens_regress)
    temp_rsq=[ModelResults_shortS60_all_multisense{1,j}.rsquared];
    idx_rsq_temp=find(temp_rsq>0.3 & temp_rsq<1);
        
    
idx_rsq_test_s60short_all_multisense{j}=idx_rsq_temp;

end

idx_rsq_test_s60short_multisense=horzcat(idx_rsq_test_s60short_all_multisense{:});
idx_rsq_test_s60short_multisense=unique(idx_rsq_test_s60short_multisense);


figure;imagesc(ZS_s60(idx_rsq_test_s60short_multisense,ZS_short_S60), [-0.5 4]);colormap hot


%%

%%% i still get some that are not responding to both stimuli... so i will
%%% finish up the cleaning with a correlation to the sound response.


Correlation_multisense=[];

temp_ZS=ZS_s60(idx_rsq_test_s60short_multisense,ZS_short_S60);
    for idx=1:size(idx_rsq_test_s60short_multisense,2)
        
            temp=corrcoef(rawregressS20(6,897:956),temp_ZS(idx,897:956));
            Correlation_multisense(1,idx)=temp(1,2);            
        
    end
    
clearvars temp_ZS 


idx_multisense=[];

Threshold=0.5;

idx_multisense=idx_rsq_test_s60short_multisense(find(Correlation_multisense>Threshold));
    

figure;imagesc(ZS_s60(idx_multisense,ZS_short_S60), [-0.5 4]);colormap hot

figure;plot(mean(ZS_s60(idx_multisense,ZS_short_S60)));

figure;scatter(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2),'.');


figure;scatter3(ROI_temp2(idx_multisense,1),ROI_temp2(idx_multisense,2),ROI_temp2(idx_multisense,3),'.');


% figure;
% scatter3(ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),1),ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),2),ROI_temp2(intersect(PerBrainRegions.Pallium.idx,idx_multisense),3),'.'); %%% this is my rotation



%%



%%% to clasify them with a correlation

Correlation_group_multisense_s60={};
counter=1;
for i=1:size(multisens_regress,1)
    Correlation_temp=[];
    for idx=1:size(ZS_s60(idx_multisense),2)
        temp_corr=corrcoef(multisens_regress(i,:),ZS_s60(idx_multisense(idx),ZS_short_S60));
            Correlation_temp(i,idx)=temp_corr(1,2);
                        
        
    end
    Correlation_group_multisense_s60{i}=Correlation_temp;
    counter=counter+1;
end

Correlation_group_multisense_mat_s60=[];
for n=1:size(multisens_regress,1)
Correlation_group_multisense_mat_s60(n,:)=Correlation_group_multisense_s60{n}(n,:);

end

%%% now if I filter looking for the max correlation
%%% 

High_corr_multisense_Nb_s60=zeros(length(Correlation_group_multisense_mat_s60),1);
for i=1:length(Correlation_group_multisense_mat_s60)
    [~,I]=max(Correlation_group_multisense_mat_s60(:,i));
    High_corr_multisense_Nb_s60(i,1)=I;
    
end

%%

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
rows=size(multisens_regress,1);counter=1;
for i=1:size(multisens_regress,1)
    
    idx_temp=find(High_corr_multisense_Nb_s60==i);
    subplot(rows,4,counter);plot(mean(ZS_s60(idx_multisense(idx_temp),ZS_short_S60),1)); %%%to plot the mean
    subplot(rows,4,counter+1);imagesc(ZS_s60(idx_multisense(idx_temp),ZS_short_S60),[0 3]);%%%for the raster plot
    subplot(rows,4,counter+2);histogram(idx_Plane_s60(idx_multisense(idx_temp))); %%%for the plane location   
    subplot(rows,4,counter+3);histogram(idx_Fish_s60(idx_multisense(idx_temp)));%%% for the fish location
    counter=counter+4;
end


save('multisense_S60.mat','idx_multisense','idx_rsq_test_s60short_multisense','Correlation_multisense','ModelResults_shortS60_all_multisense','multisens_regress','Correlation_group_multisense_mat_s60','High_corr_multisense_Nb_s60','-v7.3');


