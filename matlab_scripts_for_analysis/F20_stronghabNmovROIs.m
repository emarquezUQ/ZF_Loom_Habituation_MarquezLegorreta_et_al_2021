



%%% this script is to compare strong habituating neurons and the tail movements. 
%%% using the f20 dataset

load('All_More_BrainReg2.mat')
load('Zbrain_brainMask2D.mat')

f20_idx=load('f20_cleaned_idxs.mat');


load('ZS_N_idx_Fish_all.mat','ZS_f20','idx_Fish_f20');


%%%% to check the movement responses per fish. 
counter=1;
figure;
for i=unique(idx_Fish_f20)'
    subplot(3,4,counter);
    temp_idx=find(idx_Fish_f20==i);
    temp_idx=intersect(temp_idx,f20_idx.idx_rsq_Mov_cleaned);
    imagesc(ZS_f20(temp_idx,:));
counter=counter+1;
    
end

%% correlation between strongly habituating neurons (green) neurons and the movement neurons. 
%%%% a correlation of each fish fasthab ROIs to its movements to find neurons related

%%%% first loom is excluded because it could have too much weight on the correltion results
%%%% using spearman correlation 

fastNmovCorr=zeros(length(idx_Fish_f20(fasthab_idx)),2);
for i=1:length(idx_Fish_f20(fasthab_idx))
    
    temp_f=idx_Fish_f20(fasthab_idx(i));
    temp_idx=find(idx_Fish_f20==temp_f);
    temp_idx=intersect(temp_idx,f20_idx.idx_rsq_Mov_cleaned);
    temp_mov=mean(ZS_f20(temp_idx,:));
    
    [temp_corr temp_p]=corr(ZS_f20(fasthab_idx(i),120:end)',temp_mov(1,120:end)','Type','Spearman');
    temp_corr=temp_corr(1,1);
    temp_p=temp_p(1,1);
    
   fastNmovCorr(i,1)=temp_corr; 
    fastNmovCorr(i,2)=temp_p;  
end

figure;histogram(fastNmovCorr(:,1));

h = kstest(fastNmovCorr(:,1))
h = adtest(fastNmovCorr(:,1))
h = jbtest(fastNmovCorr(:,1))

length(find(isnan(fastNmovCorr(:,1))));%%% there are some NANs because there are fish without movements. 
nanmean(fastNmovCorr(:,1));
nanstd(fastNmovCorr(:,1));

threshold=nanmean(fastNmovCorr(:,1))+nanstd(fastNmovCorr(:,1))*1; 

%%% to selected based on corr coef value
idx_fastNmovCorr=find(fastNmovCorr(:,1)>threshold); 
idx_fastNmovCorr=fasthab_idx(idx_fastNmovCorr);

%%% where are these neurons? 
figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(ROI_temp2.f20(idx_fastNmovCorr,1),ROI_temp2.f20(idx_fastNmovCorr,2),10,'filled')
view(-90,90)

%%%% comparing them with the total strong hab
figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
hold on;
scatter(ROI_temp2.f20(fasthab_idx,1),ROI_temp2.f20(fasthab_idx,2),10,'filled')
hold on;
scatter(ROI_temp2.f20(idx_fastNmovCorr,1),ROI_temp2.f20(idx_fastNmovCorr,2),10,'filled')
view(-90,90)



%%% how many and wich fish are included? 

unique(idx_Fish_f20(idx_fastNmovCorr))


%%%% to check the movement responses per fish. 
counter=1;
figure;
for i=unique(idx_Fish_f20)'
    subplot(3,4,counter);
    temp_idx=find(idx_Fish_f20==i);
    temp_idx=intersect(temp_idx,idx_fastNmovCorr);
    imagesc(ZS_f20(temp_idx,:));
counter=counter+1;
    
end

%%%% to check the movement responses per fish in the brain. 
counter=1;
figure;
for i=unique(idx_Fish_f20)'
    subplot(3,4,counter);
    temp_idx=find(idx_Fish_f20==i);
    temp_idx=intersect(temp_idx,idx_fastNmovCorr);
    plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');
    hold on;
    scatter(ROI_temp2.f20(temp_idx,1),ROI_temp2.f20(temp_idx,2),10,'filled')
    view(-90,90)
       
counter=counter+1;
    
end


%%%% to see how many ROIs are in each brain region. 

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

fasthab_mov_brain=[];
 for brain=1:length(RegionList)   
     
      idx_brain_temp=PerBrainRegions.f20.(RegionList{brain}).idx;
      temp_idx=intersect(idx_fastNmovCorr,idx_brain_temp);
      
      if ~isempty(temp_idx)
      fasthab_mov_brain(brain)=size(temp_idx,1);
      
      else       
        fasthab_mov_brain(brain)=0;  
      end
 end
 
 figure;bar(fasthab_mov_brain);
 
 
 %%%  the proportions compared to the strong hab with in each brain region.
 fasthab_mov_brain_prop=[];
 for brain=1:length(RegionList)   
     
      idx_brain_temp=PerBrainRegions.f20.(RegionList{brain}).idx;
      temp_idx1=intersect(idx_fastNmovCorr,idx_brain_temp);
      temp_idx2=intersect(fasthab_idx,idx_brain_temp);
      
      if ~isempty(temp_idx1)
      fasthab_mov_brain_prop(brain)=size(temp_idx1,1)/size(temp_idx2,1);
      
      else       
        fasthab_mov_brain_prop(brain)=0;  
      end
 end
 
 figure;bar(fasthab_mov_brain_prop);
 
%%%% getting the ROIs locations for Unity
CSV_temp=ROI_temp2.f20(idx_fastNmovCorr,:);

    CSV_temp(:,3)=CSV_temp(:,3); %%% why?

    CSV_temp(:,4)=1;

    filename=strcat('__Coords_fastNmovCorr_f20_','.csv');

    csvwrite(filename,CSV_temp);


