

%%% this script is to look at the difference in response to the first loom,
%%% second and to the 10th loom to make a heat map of this change. I will select 
%%% ROIs that responded to the first loom and then look how they changed.
%%% I will select this ROIs with a linear regression to the first loom. 

%%% I will do it in f20s 


load('final_F20_step1.mat','ZS_f20');


load('All_More_BrainReg2.mat');


%%

%%% this is if I want to get first the ROIs that reacted when the first
%%% loom presentation happen based on the strenght of the response
%%% filtering by 2 units of the z-score. then I get the delta of the first
%%% response and substract it to the delta of the 10th response
idx_resp=[];
for i=1:size(ZS_f20,1)

%if  max(ZS_f20(i,60:80))>2   
idx_resp(i)=max(ZS_f20(i,60:80))>2;
%elseif min(ZS_f20(i,60:80))<-2 
%idx_resp(i)=min(ZS_f20(i,60:80))<-2;
%end

delta_1(i)=max(ZS_f20(i,60:80))-mean(ZS_f20(i,10:59));
delta_2(i)=max(ZS_f20(i,105:140))-mean(ZS_f20(i,10:59));
delta_10(i)=max(ZS_f20(i,420:440))-mean(ZS_f20(i,10:59));

deltahab2(i)=delta_1(i)-delta_2(i);
deltahab10(i)=delta_1(i)-delta_10(i);
 
ratiohab2(i)=delta_2(i)/delta_1(i);
ratiohab10(i)=delta_10(i)/delta_1(i);

end
% 
% idx_resp=find(idx_resp);
% 
% figure;imagesc(ZS_f20(idx_resp,:)); colormap('hot');
% figure;histogram(deltahab2);
% figure;histogram(deltahab10);

%%%%%  

%%

%%%% generating a regressor

Stimuli=zeros(1,100);
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=60;
for i=1%:10
    Stimuli(i,(idxStart+(i-1)*60):(idxStart+(i-1)*60)+size(GCaMP6,1)-1)=GCaMP6;
end

figure;plot(Stimuli);

%%
%%% finally, I will try with a linear regression to the first loom
  
ModelResults=[];
parfor i=1:size(ZS_f20,1)  %%%parfor is to do it in parallel
    
    %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
    %%vector Y as a response variable and the columns of the matrix X as
    %%%predictor variables, performs stepwise regression, and returns the
    %%final result as the linear model LM.
    mdl=fitlm(Stimuli',ZS_f20(i,1:100));
   
    %%%this is to put the results in the ModelResulsts variable
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
     
end

rsquare_loom=[ModelResults.rsquared];%%%to take te rsquared field from ModelResults and make a variable with them

figure; histogram(rsquare_loom);
std_rsq=std(rsquare_loom);

idx_rsq=find(rsquare_loom>2*std_rsq & rsquare_loom<1); %%%then select the rsquare threshold
figure; 
imagesc(ZS_f20(idx_rsq,:), [-0.5 4]);colormap hot %%%to plot them in a raster plot

%%% with 0.5
idx_rsq=find(rsquare_loom>0.5 & rsquare_loom<1); %%%
figure; 
imagesc(ZS_f20(idx_rsq,:), [-0.5 4]);colormap hot


rsq2d1_2=deltahab2(idx_rsq);
Qs_rsq_2 = quantile(rsq2d1_2,[0.025 0.25 0.50 0.75 0.975]);
figure;histogram(rsq2d1_2);
temp_std=std(rsq2d1_2)
figure;histogram(rsq2d1_2(find(rsq2d1_2<-2*temp_std)))
figure;plot(mean(ZS_f20(idx_rsq(find(rsq2d1_2<-2*temp_std)),:))); %%% in this case it looks like noise
figure;imagesc(ZS_f20(idx_rsq(find(rsq2d1_2<-2*temp_std)),:)); 


rsq2d1_10=deltahab10(idx_rsq);
Qs_rsq_10 = quantile(rsq2d1_10,[0.025 0.25 0.50 0.75 0.975]);
temp_std=std(rsq2d1_10)
figure;histogram(rsq2d1_10(find(rsq2d1_10<-2*temp_std)))
figure;plot(mean(ZS_f20(idx_rsq(find(rsq2d1_10<-2*temp_std)),:))); %%% this is weird. I find some negative deltas (more response in the tenth loom) but when I plot them it doenst seem so 


rsq2r1_2=ratiohab2(idx_rsq);
Qs_rsq_ratio_2 = quantile(rsq2r1_2,[0.025 0.25 0.50 0.75 0.975]);
figure;histogram(rsq2r1_2); %%% i am geting some few but very weird results very far from 0-1. like -280.38 or 808.0434. they are most probably artifacts
max(rsq2r1_2)
high_idx=find(rsq2r1_2>1);
high=rsq2r1_2(high_idx);
figure;scatter(ROI_temp2.f20(idx_rsq(high_idx),1),ROI_temp2.f20(idx_rsq(high_idx),2),10,high,'filled');colormap('jet');colorbar; %caxis([Qs_rsq2(1) Qs_rsq2(5)]);



rsq2r1_10=ratiohab10(idx_rsq);
Qs_rsq_ratio_10 = quantile(rsq2r1_10,[0.025 0.25 0.50 0.75 0.975]);


del_1_Nrsq=delta_1(idx_rsq);
Qs_d1 = quantile(del_1_Nrsq,[0.025 0.25 0.50 0.75 0.975]);

del_2_Nrsq=delta_2(idx_rsq);
Qs_d2 = quantile(del_2_Nrsq,[0.025 0.25 0.50 0.75 0.975]);

del_10_Nrsq=delta_10(idx_rsq);
Qs_d10 = quantile(del_10_Nrsq,[0.025 0.25 0.50 0.75 0.975]);



figure;
subplot(1,5,1)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,del_1_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,2)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,del_2_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,3)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,del_10_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,4)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,rsq2r1_2,'filled');colormap('jet');colorbar; caxis([Qs_rsq_ratio_2(1) Qs_rsq_ratio_10(5)]);
subplot(1,5,5)
scatter(ROI_temp2.f20(idx_rsq,2),ROI_temp2.f20(idx_rsq,1),10,rsq2r1_10,'filled');colormap('jet');colorbar; caxis([Qs_rsq_ratio_2(1) Qs_rsq_ratio_10(5)]);

%% to clean the ROIs outside the brain. of the linear regression filtering

%%% i first make a mask with all the brain regions
load('Zbrain_Masks.mat');


Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
 
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

figure;scatter(Zbrain_AllMask(:,1),Zbrain_AllMask(:,2),'.');


idx_brain=ismember(ROI_temp2.f20,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks
 
idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain

idx_rsq_1stLoom_cleaned=intersect(idx_brain2,idx_rsq);

rsq2d1_2=deltahab2(idx_rsq_1stLoom_cleaned);
Qs_rsq_2 = quantile(rsq2d1_2,[0.025 0.25 0.50 0.75 0.975]);

rsq2d1_10=deltahab10(idx_rsq_1stLoom_cleaned);
Qs_rsq_10 = quantile(rsq2d1_10,[0.025 0.25 0.50 0.75 0.975]);

rsq2r1_2=ratiohab2(idx_rsq_1stLoom_cleaned);
Qs_rsq_ratio_2 = quantile(rsq2r1_2,[0.025 0.25 0.50 0.75 0.975]);

rsq2r1_10=ratiohab10(idx_rsq_1stLoom_cleaned);
Qs_rsq_ratio_10 = quantile(rsq2r1_10,[0.025 0.25 0.50 0.75 0.975]);

del_1_Nrsq=delta_1(idx_rsq_1stLoom_cleaned);
Qs_d1 = quantile(del_1_Nrsq,[0.025 0.25 0.50 0.75 0.975]);

del_2_Nrsq=delta_2(idx_rsq_1stLoom_cleaned);
Qs_d2 = quantile(del_2_Nrsq,[0.025 0.25 0.50 0.75 0.975]);

del_10_Nrsq=delta_10(idx_rsq_1stLoom_cleaned);
Qs_d10 = quantile(del_10_Nrsq,[0.025 0.25 0.50 0.75 0.975]);



figure;
subplot(1,5,1)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,del_1_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,2)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,del_2_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,3)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,del_10_Nrsq,'filled');colormap('hot');colorbar; caxis([Qs_d10(1) Qs_d1(5)]);
subplot(1,5,4)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,rsq2r1_2,'filled');colormap('jet');colorbar; caxis([Qs_rsq_ratio_10(1) Qs_rsq_ratio_10(5)]);
subplot(1,5,5)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,rsq2r1_10,'filled');colormap('jet');colorbar; caxis([Qs_rsq_ratio_10(1) Qs_rsq_ratio_10(5)]);


%%% for the colorbar of the substractions. I made a a colormap of 18 values
%%% (0:17) and upload it to unity. so I need that range for my colorbar. 

%c = jet(18);
%c = plasma(18);
%c = inferno(18);
c = viridis(18);


% filename=strcat('jet_colormap2.csv');
% filename=strcat('plasma_colormap.csv');
% filename=strcat('inferno_colormap.csv');
 filename=strcat('viridis_colormap.csv');
 
 csvwrite(filename,c);

figure;
 subplot(1,2,1)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,rsq2r1_2,'filled');colormap(plasma);colorbar; caxis([0 1]);
subplot(1,2,2)
scatter(ROI_temp2.f20(idx_rsq_1stLoom_cleaned,2),ROI_temp2.f20(idx_rsq_1stLoom_cleaned,1),10,rsq2r1_10,'filled');colormap(plasma);colorbar; caxis([0 1]);


%% generate coordenates for Unity

%%%% I will put the coords of x,y and z. and then the rsq value and then
%%%% the deltas or ratios. 


%%% for the linear regression results

rsquare_loom_short=rsquare_loom(idx_rsq_1stLoom_cleaned);

    idx_temp2=idx_rsq_1stLoom_cleaned;
    
    CSV_temp=ROI_temp2.f20(idx_temp2,:);

    CSV_temp(:,3)=CSV_temp(:,3); 

    CSV_temp(:,5)=rsquare_loom_short;
    
%     CSV_temp(:,4)=del_1_Nrsq;
%     CSV_temp(:,4)=del_2_Nrsq;
%     CSV_temp(:,4)=del_10_Nrsq;
%     CSV_temp(:,4)=rsq2d1_2;
%     CSV_temp(:,4)=rsq2d1_10;
%    CSV_temp(:,4)=rsq2r1_2;
     CSV_temp(:,4)=rsq2r1_10;
    
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d1_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d2_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d10_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_df2_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_df10_new.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_rf2_new.csv');
    filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_rf10_new.csv');
    
    csvwrite(filename,CSV_temp);

    save('f20_1stloom_deltaNreg.mat');

%%

%%%  with the deltas 

idx_temp2=idx_rsq_1stLoom_cleaned;
    
    CSV_temp=ROI_temp2.f20(idx_temp2,:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% 

    CSV_temp(:,4)=rsquare_loom_short;
    
    %CSV_temp(:,5)=del_1_Nrsq;
    %CSV_temp(:,5)=del_2_Nrsq;
    CSV_temp(:,5)=rsq2d1_10;
    
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d1.csv');
    %filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_d2.csv');
    filename=strcat('__Coords_f20_idx_rsq_1stLoom_cleaned_df.csv');
    
    csvwrite(filename,CSV_temp);
    
    

