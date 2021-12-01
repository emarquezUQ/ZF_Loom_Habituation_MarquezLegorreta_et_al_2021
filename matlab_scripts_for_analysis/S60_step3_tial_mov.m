

%%% this script is to get the ROIs that correlated to the movements. I do
%%% it with a linear regression using rsq=0.2 as threshold.

%%% NOTE: In the long movies (F60 and S60) detecting the tail movements in
%%% the calcium traces was not as effective with this method. 

%%% for fish s60

%%%% first you need to load what you need

load('BrainReg_S60.mat')
load('final_S60_step1.mat','ZS_s60','idx_Fish_s60','idx_Plane_s60','rawregressS60','idx_rsq_test_s60short','High_corr_Nb_s60','High_corr_Nb_s60_short','ZS_short_S60');
load('Zbrain_Masks.mat');

load('allfish_s60looms_tailmov.mat','-mat');

%%

allfish_s60looms_tailmov=cell2mat(allfish_s60looms_tailmov(2:end,:));


%%%first to generate the regressors for the movements. 
Fish_movements={}; %%%first we make a structure to put the data in
for fish=1:size(allfish_s60looms_tailmov,2) %%%we make a variable 'fish' as long as the columns from the Behavioural_responses variable which has the data copied from excel
    for i=1:2 %%%% Note: here I had a loop to work with the different strenght responses but we decided to only take the strong responses (1 and 2).
        %%% cause I have resposnes from 1-4 depending on the strenght of the moment.
        
        temp=find(allfish_s60looms_tailmov(:,fish)==i); %%%to find the type of responses in turn. 
        temp=round(temp/5);%%%% i this is to round it to the frame rate that we use for the SPIM videos. the tail movies have a 10 Hz and the SPIM is 2Hz. 
        temp(temp<1)=[];
        temp(temp>(length(allfish_s60looms_tailmov)-1))=[];        
        Fish_movements{fish,i}=temp;
        
        
    end
    
    %%% this is loop is to locate the moments the loom is ON.
    %%% but I need to make a vector with the momennts the loom is on. 
    if fish==13 %%at the moment is 13 (for s60)
        temp=find(allfish_s60looms_tailmov(:,fish)==2);
        temp=round(temp/5);
        temp(temp<1)=[];
        temp(temp>(length(allfish_s60looms_tailmov)-1))=[];        
        Fish_movements{fish,1}=temp;
    end
end


%%%this is to put Gcamp spikes to the fish movements. 
Movements=zeros(size(allfish_s60looms_tailmov,2)-1,1,size(ZS_s60,2)); %%%to mae a 3D array the size of the number of fish (18 at the moment), the type of movement (1) and the lenght of the SPIM movies (depending on the variable). 
%GCaMP6=[5.13796058542217,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502]';
GCaMP6=[0,0.5,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
for fish=1:size(allfish_s60looms_tailmov,2)-1
    for movement=1  %movement=1:4
        idx=Fish_movements{fish,movement};%%%%to take the timepoints of the movements.  
        if idx
            for i=idx'
                Movements(fish,movement,i:i+size(GCaMP6,1)-1)=GCaMP6/6;%%% to place the Gcamp spikes in the places where the movements occur
            end
        end
    end
end
Movements=Movements(:,:,1:size(ZS_s60,2));


%%

%%% this is the fish list. 

% Fish_with_behav_f20=[1 7 13 17 21 25 29 34 38 40 44 48];
% 
% Fish_with_behav_f60=[3 4 8 9 15 19 26 32 37 41 47];
% 
% Fish_with_behav_s20=[5 10 14 18 22 24 28 33 36 43 45];
% 
Fish_with_behav_s60=[6 11 12 16 20 23 27 30 31 35 39 42 46];

Fish_with_behav_s60=intersect(unique(idx_Fish_s60),Fish_with_behav_s60);

%%% to show where when all the looms are presented


loomtimes=zeros(1,size(ZS_s60,2)); %%%to mae a 3D array the size of the number of fish (18 at the moment), the type of movement (1) and the lenght of the SPIM movies (depending on the variable). 

idx=Fish_movements{13,1};%%%%to take the timepoints of the movements.  
       
  for i=idx'
   loomtimes(1,i)=0.5;%%% to place the Gcamp spikes in the places where the movements occur
  end
        
 


figure;plot(loomtimes(1,ZS_short_S60));


%%

%%% linear regression with the all the movements to the
%%%%Gcamps of the actual brain imaging of each fish.I am taking away the
%%%%first 3 looms. for s60 it would be from 180 onwards. This is because
%%%%otherwise the strong first responses dominate the results and it's hard
%%%%to pick up the later movements. 

ModelResults_allmovPerFish={};
counter=1;
rsquare_loom_allmov={};
for fish=1:length(Fish_with_behav_s60)
   if ismember(Fish_with_behav_s60(fish),unique(idx_Fish_s60))  
%for fish=Fish_with_behav_s60
  % if ismember(fish,unique(idx_Fish_s60))  
    temp_ZS=ZS_s60(find(idx_Fish_s60==Fish_with_behav_s60(fish)),ZS_short_S60);ModelResults_temp=[];
    %temp_ZS=ZS_s60(find(idx_Fish_s60==fish),:);ModelResults_temp=[];
        parfor i=1:size(temp_ZS,1)  %%%parfor is to do it in parallel
    
        %%% this is to to the linear regression: LM = stepwiselm(X,Y) fits a linear regression model using the column
        %%vector Y as a response variable and the columns of the matrix X as
        %%%predictor variables, performs stepwise regression, and returns the
        %%final result as the linear model LM.
        
        mdl=fitlm(squeeze(squeeze(Movements(find(Fish_with_behav_s60==Fish_with_behav_s60(fish)),ZS_short_S60)))',temp_ZS(i,:));
         
        
        %%%this is to put the results in the ModelResulsts variable
        ModelResults_temp(i).coef=mdl.Coefficients;
        ModelResults_temp(i).MSE=mdl.MSE;
        ModelResults_temp(i).Fitted=mdl.Fitted;
        ModelResults_temp(i).rsquared=mdl.Rsquared.Adjusted;
        end
        
        ModelResults_allmovPerFish{counter}=ModelResults_temp;
       
        rsquare_loom_allmov{counter}=[ModelResults_temp.rsquared];
        
         counter=counter+1;
        
   else
   end
end


%%

%%%this is to get the idx of the ROIs that passed the rsq2 thereshold that we choose to the movements
%%%(regardless of the loom). it solves the problem of the indexing. 
%%% 0.1-0.2 is a good threshold 


idx_rsq_Mov={};
rsquare_loom_allmov2={};
counter=1;
for fish=1:length(Fish_with_behav_s60)
   if ismember(Fish_with_behav_s60(fish),unique(idx_Fish_s60))  
    temp_idx=find(idx_Fish_s60==Fish_with_behav_s60(fish)); idx_Mov_temp=[];   
    rsq_temp=[ModelResults_allmovPerFish{1,counter}.rsquared];
    idx_Mov_temp=temp_idx(find(rsq_temp>0.2 & rsq_temp<1),:);
    idx_rsq_Mov{counter}=idx_Mov_temp;
    rsquare_loom_allmov2{counter}=rsq_temp;
    counter=counter+1;
   else
   end    
end

clear idx_Mov_temp temp_idx rsq_temp counter

idx_rsq_Mov=vertcat(idx_rsq_Mov{:});
figure; imagesc(ZS_s60(idx_rsq_Mov,ZS_short_S60), [-0.5 4]);colormap hot


rsquare_loom_allmov2=horzcat(rsquare_loom_allmov2{:});

figure;histogram(rsquare_loom_allmov2);


save('Tail_mov_S60.mat','idx_rsq_Mov','rsquare_loom_allmov2','ModelResults_allmovPerFish','loomtimes','Fish_with_behav_s60','Movements','-v7.3');



%%
%%%%and this is to plot the means of the traces of the ROIs that had an r2
%%%%value above the threshold set above.
%%%% It makes a figure with a subplot of each fish. 
 


   
    counter=1;
    figure; 
    xplot=floor(sqrt(length(Fish_with_behav_s60)));yplot=ceil(length(Fish_with_behav_s60)/xplot);
    for fish=1:length(Fish_with_behav_s60)
        if ismember(Fish_with_behav_s60(fish),unique(idx_Fish_s60)) %%isempty(Correlation_movement_post3{fish})~=1 % idx_Fish_s60==Fish_with_behav_s60(fish)
        temp_ZS=ZS_s60(find(idx_Fish_s60==Fish_with_behav_s60(fish)),ZS_short_S60);        
        rsq_temp=[ModelResults_allmovPerFish{1,counter}.rsquared];
        subplot(xplot,yplot,counter);plot(mean(temp_ZS(find(rsq_temp>0.1 & rsq_temp<1),:),1));
        hold on
%       Movements2=squeeze(Movements(:,i,:)); %%% this was already done before
%       Movements2=squeeze(Movements2);
        plot(Movements2(find(Fish_with_behav_s60==Fish_with_behav_s60(fish)),ZS_short_S60))
        plot(loomtimes(1,ZS_short_S60));
        hold off
        
        %mdl=stepwiselm(Stimuli',mean(temp_ZS(find(Correlation(i,:)>Threshold),:),1),'linear','Criterion','adjrsquared','Intercept',true,'Upper','interactions','Verbose',0);
        %title(num2str([mdl.Rsquared.Adjusted]));
        title(strcat(' nb of ROIs : ', num2str(length(find(rsq_temp>0.1 & rsq_temp<1)))));
        
         counter=counter+1;
        else
        end
    end

    clear temp_ZS rsq_temp counter



