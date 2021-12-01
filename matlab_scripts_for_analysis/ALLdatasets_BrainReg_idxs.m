%%%%% this script has two parts where we organized the location of ROIs in
%%%%% some brain regions of interest. The output datasets are then used in
%%%%% multiple other scripts for analysis. 


%%%% this the first script is to collect the idxs of some brain regions. 


datasets=['f20'; 'f60'; 's20'; 's60'];

load('Zbrain_Masks.mat');
BrainReg_F20=load('BrainReg_F20.mat');
BrainReg_F60=load('BrainReg_F60.mat');
BrainReg_S20=load('BrainReg_S20.mat');
BrainReg_S60=load('BrainReg_S60.mat');

ROI_temp2.f20=BrainReg_F20.ROI_temp2;
ROI_temp2.f60=BrainReg_F60.ROI_temp2;
ROI_temp2.s20=BrainReg_S20.ROI_temp2;
ROI_temp2.s60=BrainReg_S60.ROI_temp2;


PerBrainRegions=struct();
RegionList={'Telencephalon','Pallium','Subpallium','Thalamus','Cerebellum','Semicircularis','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Msystem'};


%%

%%% this is to look at brain regions

for data=1:4
%progressbar;
for i=1:length(RegionList)
    %progressbar(i/length(RegionList2));
    regionName=RegionList{i};
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
    IsInBrainRegion=ismember(ROI_temp2.(datasets(data,:)),Mask,'rows');
    PerBrainRegions.(datasets(data,:)).(regionName).idx=find(IsInBrainRegion==1);    
end

end


save('All_More_BrainReg.mat','PerBrainRegions','RegionList','ROI_temp2','-v7.3');



%% with more brain regions

%%%% to add the caudal hypothalamus, EmT, M1, M3 (this one might not be in the atlas), divide
%%%% the thalamus in dorsal and ventral and the M2 (future preglomerular), some of the AF (4, 6 and 8).


PerBrainRegions=struct();
RegionList={'Telencephalon','Pallium','Subpallium','Caudal Hypothalamus','Dorsal Thalamus','Ventral Thalamus','Eminentia Thalami','(M1)','(M2)','AF4','AF6','AF8','Thalamus','Cerebellum','Semicircularis','Tectum','Tegmentum','Habenula','Pretectum','MON','Hindbrain','Msystem'};


%%

%%% this is to look at brain regions

for data=1:4
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
    IsInBrainRegion=ismember(ROI_temp2.(datasets(data,:)),Mask,'rows');
    if strcmp(regionName,'(M1)')
        regionName='M1';
    elseif strcmp(regionName,'(M2)')
        regionName='M2';
    end
    PerBrainRegions.(datasets(data,:)).(regionName(~isspace(regionName))).idx=find(IsInBrainRegion==1);    
end

end


save('All_More_BrainReg2.mat','PerBrainRegions','RegionList','ROI_temp2','-v7.3');


