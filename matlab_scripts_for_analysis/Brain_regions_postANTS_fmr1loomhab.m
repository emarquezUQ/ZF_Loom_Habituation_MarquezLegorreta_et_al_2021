


%%%% this script is to get the ANTs results and asign to each ROI its
%%%% proper indexing. I had a few problems with a previous verison but now
%%%% its working well. What i did was to first generate another list of
%%%% ROIs per fish where each ROI was tagged with its indexing (see script
%%%% ROIs_for_ANTS_Emm2_fmr1loomhab_tagged.m). Then I could use that
%%%% indexing to sort them properly. Then in this script I also use the
%%%% zbrain masks to locate the ROIs in their brain areas. 


load('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');

load('s20_good_idx_Fish.mat','idx_Fish');

%%% to get the clusters from the Kmeans done to ALL the ROIs to check
%%% distribution of ROIs
load('s20_postKmeans_CN.mat','idxKmeans_ZS_CN');

load('s20_fmr1_loomhab_CN_part2.mat','GoodBetas_ZS_CN_selected','idx_rsq');

GoodBetas_ZS_CN_selected=[13 21 29 32 35 36 38 40 43 44];

%%

%%% this part was to check for the new ROIs for ANTs

CSV_Files2=dir('_ROIsFish*.csv');
%CSV_Files2=dir('_ROIsFish1.csv'); %% to test with first one
%M = csvread('_ROIsFish1.csv');
ROIs_new=struct();truth_2=[];
for i=1:length(CSV_Files2)
    temp_fishROI=csvread(CSV_Files2(i).name,0); %%%% NOTE: to read the csv file from the begining I need to put a 0 here instead of 1...  
    Fishname=regexp(CSV_Files2(i).name,'_ROIsFish(\d+)_tagged.csv','tokens');Fishname=Fishname{1}{1};    
    ROIs_new(i).name=Fishname;    
    %ROIs_new(i).coord=temp_fishROI(:,1:3);
    ROIs_new(i).idx=temp_fishROI(:,5);
    ROIs_new(i).good_idx=temp_fishROI(:,6);
    truth_2(i)=size(temp_fishROI,1)==sum(idx_Fish==str2num(Fishname))
end
clearvars i temp CSV_Files Fishname

%%% this is to check for errors. it looks if the number of ROIs in the csv
%%% files is the same than the goodnumber ROIs of each fish. If there are
%%% errors they need to be fixed before to move on. 


for i=1:length(truth_2)
    if ~truth_2(i)
        strcat('Fish',num2str(ROIs_new(i).name),': error = ',num2str(length(ROIs_new(i).idx)-sum(idx_Fish==str2num(ROIs_new(i).name))))
    end    
end

clearvars i temp_fishROI CSV_Files2 Fishname 



%%



%%%% now, this is the part where i used the warped datasets to get the
%%%% brain locations of the ROIs.

Fish_list=unique(idx_Fish);

CSV_Files=dir('_2Warped*.csv');
%CSV_Files=dir('_ROIsFish*.csv'); %%% to check if it was working
ROIs=struct();truth=[]; counter=1;
for i=1:length(CSV_Files)
      
    Fishname=regexp(CSV_Files(i).name,'_2Warped(\d+)_fish(\d+).csv','tokens');Fishname=strcat('2018',Fishname{1}{1},Fishname{1}{2});    
    %Fishname=regexp(CSV_Files(i).name,'_ROIsFish(\d+)_tagged.csv','tokens');Fishname=Fishname{1}{1}; %%% this was to check if it
    %was working 
   
    
    if ismember(str2num(Fishname),Fish_list)
    
    temp_warp=csvread(CSV_Files(i).name,1); %%%% NOTE: in this case I put 1 cause the first row are the heathers.   
     
    ROIs(counter).name=Fishname;    
    ROIs(counter).coord=temp_warp(:,1:3);
    ROIs(counter).idx=temp_warp(:,5);
    %ROIs(counter).good_idx=temp_warp(:,6);  %%% this was to check if it
    %was working
    
    truth(counter)=size(temp_warp,1)==sum(idx_Fish==str2num(Fishname)) %%% important!! i need to have an idx_Fish with all the fish.
    
     counter=counter+1;
    else
    end
      
end


%%% this is to get the ROIs of all the fish (or the ones i used the warped files off) and put them all together.
i=1;ROI_pool=ROIs(i).coord;
for i=2:length(ROIs)
    ROI_pool=[ROI_pool; ROIs(i).coord];
end

%%% this is to get the good_idx to try to sort the ROIs properly
i=1;ROI_idx=ROIs_new(i).good_idx;
for i=2:length(ROIs_new)
    ROI_idx=[ROI_idx; ROIs_new(i).good_idx];
end



%% sorting the ROIs

ROI_fish(ROI_idx,:)=ROI_pool; %%% Sorting ROIs based on good_idx
ROI_temp=round(ROI_fish);           %%% rounding the values of the coordenates because in the zbrain masks they are integers 
ROI_temp(:,1)=round(ROI_fish(:,2)); %%% this is to make the x axis the new y axis in case i need to rotate it
ROI_temp(:,2)=round(ROI_fish(:,1));  %%% this is to make the y axis the new x axis in case i need to rotate it 

figure;scatter(ROI_temp(:,1),ROI_temp(:,2));

figure;
scatter(ROI_temp(intersect(idx_rsq,find(idxKmeans_ZS_CN==GoodBetas_ZS_CN_selected(8))),1),ROI_temp(intersect(idx_rsq,find(idxKmeans_ZS_CN==GoodBetas_ZS_CN_selected(8))),2),'filled');
hold on;
scatter(ROI_temp(intersect(idx_rsq,find(idxKmeans_ZS_CN==GoodBetas_ZS_CN_selected(6))),1),ROI_temp(intersect(idx_rsq,find(idxKmeans_ZS_CN==GoodBetas_ZS_CN_selected(6))),2),'filled');
hold on;
scatter(ROI_temp(intersect(idx_rsq,find(idxKmeans_ZS_CN==GoodBetas_ZS_CN_selected(9))),1),ROI_temp(intersect(idx_rsq,find(idxKmeans_ZS_CN==GoodBetas_ZS_CN_selected(9))),2),'filled');

%%

%%% I need to get the Zbrain_Masks

load('Zbrain_Masks.mat');


%%

ROI_temp2=ROI_temp;  %%%  making this variable  to be able to play with it.
 
figure;scatter(ROI_temp2(:,1),ROI_temp2(:,2));


%%% checking if I need to do a rotation

figure;scatter(ROI_temp(find(idx_Fish==Fish_list(1)),1),ROI_temp(find(idx_Fish==Fish_list(1)),2)); 

figure;scatter(ROI_temp2(find(idx_Fish==Fish_list(1)),1),ROI_temp2(find(idx_Fish==Fish_list(1)),2)); 


%%% it seems they are looking right already.


figure;scatter(Zbrain_Masks{294,3}(:,1),Zbrain_Masks{294,3}(:,2)); %%% this is telencephalon in Zbrain
%%% and I need them to look to the right. 

figure;scatter(ROI_temp2(find(idx_Fish==Fish_list(1)),1),ROI_temp2(find(idx_Fish==Fish_list(1)),2)); %%% this is my rotation
%%%as I have it at the moment is looking to the right. 
hold on;
scatter(Zbrain_Masks{294,3}(:,1),Zbrain_Masks{294,3}(:,2)); %%% this is telencephalon in Zbrain


%%% now checking on the z axis that the ROIs are not outside the range. 

%%% zbrain OT
figure;histogram(Zbrain_Masks{105,3}(:,3));
%%% so the OT reaches till a bit more than 130. the maximum in zbrain is
%%% 137


max(Zbrain_Masks{105,3}(:,3))


%%% one of the fish
figure;histogram(ROI_temp2(find(idx_Fish==Fish_list(1)),3));
%%% but I am getting much higher values... from -49 to 129185.
%%% but it seems to vary from fish to fish.
%%% so i will need to adjust it to each fish. 

%%% I will test a formula for rescaling. 

min(ROI_temp2(:,3))
max(ROI_temp2(find(idx_Fish==Fish_list(1)),3))
min(ROI_temp2(find(idx_Fish==Fish_list(1)),3))

old_Zfish1=ROI_temp2(find(idx_Fish==Fish_list(1)),3);
new_Zfish1=(((137-0)*(old_Zfish1-min(old_Zfish1)))/(max(old_Zfish1)-min(old_Zfish1)))+0; 

%%% or:

new_Zfish1=(old_Zfish1-min(old_Zfish1))/(max(old_Zfish1)-min(old_Zfish1))*(137-0)+0; 

figure;histogram(new_Zfish1);

min(new_Zfish1)
max(new_Zfish1)

%%% it worked, so now a loop for each fish

%%% NOTE: my fish dont go as deep as the zbrain stack... so i will try to
%%% fit it from 45-135 instead. it seems that it didnt work... maybe cause
%%% is already taken into account after the warping. I will try again with
%%% 0-137


ROI_fish2=ROI_fish;
ROI_fish2 = ROI_fish2(any(ROI_fish2,2),:);%%% 


for i=1:length(Fish_list)
    
    temp_fish=Fish_list;
    temp_fish=temp_fish(i);
    temp_Zvalues=ROI_fish(find(idx_Fish==temp_fish),3);
    new_Zvalues=(temp_Zvalues-min(temp_Zvalues))/(max(temp_Zvalues)-min(temp_Zvalues))*(137-0)+0;
    ROI_fish2(find(idx_Fish==temp_fish),3)=new_Zvalues;
end


clear temp_fish temp_Zvalues new_Zvalues


ROI_temp2(:,3)=round(ROI_fish2(:,3));
 ROI_temp2(:,1)=round(ROI_fish2(:,2)); 
 ROI_temp2(:,2)=round(ROI_fish2(:,1)); 

%%% checking rotation and Z values


figure;scatter(ROI_temp2(find(idx_Fish==Fish_list(1)),1),ROI_temp2(find(idx_Fish==Fish_list(1)),2)); %%% this is my rotation
%%%as I have it at the moment is looking to the right. 

figure;histogram(ROI_temp2(:,3));

min(ROI_temp2(:,3))
max(ROI_temp2(:,3))



figure;scatter3(ROI_temp2(find(idx_Fish==Fish_list(1)),1),ROI_temp2(find(idx_Fish==Fish_list(1)),2),ROI_temp2(find(idx_Fish==Fish_list(1)),3)); %%% this is my rotation
%%%as I have it at the moment is looking to the right. 
hold on;
scatter3(Zbrain_Masks{294,3}(:,1),Zbrain_Masks{294,3}(:,2),Zbrain_Masks{294,3}(:,3)); %%% this is telencephalon in Zbrain

%%%% it seems that it worked.


%%

%%%% Finding ROIs in brain regions

PerBrainRegions=struct();
RegionList={'Telencephalon','Pallium','Subpallium','Caudal Hypothalamus','Dorsal Thalamus','Ventral Thalamus','Eminentia Thalami','(M1)','(M2)','AF4','AF6','AF8','Thalamus','Cerebellum','Semicircularis','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Msystem'};

%RegionList={'Thalamus','Cerebellum','Semicircularis','Telencephalon','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain'};
progressbar;
for i=1:length(RegionList)
    progressbar(i/length(RegionList));
    regionName=RegionList{i}%;
    if strcmp(regionName,'Telencephalon')
        Mask=Zbrain_Masks{294,3};
    elseif strcmp(regionName,'Hindbrain')
        Hindbrain_Mask=Zbrain_Masks{259,3};
        Mask=Zbrain_Masks{131,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove cerebellum
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Zbrain_Masks{295,3};
        IsInEyes_temp=ismember(Hindbrain_Mask,Mask,'rows');IsInEyes_temp=find(IsInEyes_temp==1);%remove MON
        Hindbrain_Mask(IsInEyes_temp,:)=[];
        Mask=Hindbrain_Mask;
    elseif strcmp(regionName,'Msystem')
        Mask=[];
        Mask_temp=[];
        Msystem_masks=[184 186 187];
        for j=1:3
           %Mask_temp=Zbrain_Masks{Msystem_masks(j),3};
           Mask=vertcat(Mask,Zbrain_Masks{Msystem_masks(j),3});
        end
    
        clear Mask_temp
    else
        Mask=[];
        IndexC=strfind({Zbrain_Masks{:,2}}, regionName);
        IndexC=find(not(cellfun('isempty', IndexC)));
        for j=IndexC
            if isempty(Mask)
                Mask=Zbrain_Masks{j,3};
            else
                Mask=vertcat(Mask,Zbrain_Masks{j,3});
            end
        end
    end
    Mask=unique(Mask,'rows');
    IsInBrainRegion=ismember(ROI_temp2,Mask,'rows');
    if strcmp(regionName,'(M1)')
        regionName='M1';
    elseif strcmp(regionName,'(M2)')
        regionName='M2';
    end
    PerBrainRegions.(regionName(~isspace(regionName))).idx=find(IsInBrainRegion==1);    
end


%%

figure;scatter(ROI_temp2(:,1),ROI_temp2(:,2));
hold on;
scatter(ROI_temp2(PerBrainRegions.Pallium.idx,1),ROI_temp2(PerBrainRegions.Pallium.idx,2));  %%% the brain location seems to be working
hold on;
scatter(ROI_temp2(idx_rsq,1),ROI_temp2(idx_rsq,2));  

%% cleaning the ROIs from nonbrain ones

Zbrain_AllMask=vertcat(Zbrain_Masks{[1:1:77 79:1:294],3});
 
Zbrain_AllMask=unique(Zbrain_AllMask,'rows');

figure;scatter(Zbrain_AllMask(:,1),Zbrain_AllMask(:,2),'.');

%Removing the eyes

idx_brain=ismember(ROI_temp2,Zbrain_AllMask,'rows'); %% to find the ROIs inside the brain masks

idx_brain2=find(idx_brain);  %% now i am getting the indexes of the ROIs inside teh brain


figure;
scatter(ROI_temp2(idx_rsq,1),ROI_temp2(idx_rsq,2));
hold on;
scatter(ROI_temp2(intersect(idx_brain2,idx_rsq),1),ROI_temp2(intersect(idx_brain2,idx_rsq),2));

idx_rsq_cleaned=intersect(idx_brain2,idx_rsq);



%%
save('fmr1loomhab_BrainRegNclean.mat','PerBrainRegions','RegionList','ROI_temp2','idx_rsq_cleaned','-v7.3');



