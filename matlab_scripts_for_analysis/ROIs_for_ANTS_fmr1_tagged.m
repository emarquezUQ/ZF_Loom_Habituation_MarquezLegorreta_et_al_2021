
%%% this scritp is to take the coordinates of ROIs of the fmr1 loomhab fish. 
%%% the purpose is to then warp them to the zbrain with ANTS.
%%% It's based on the script ROIs_for_ANTS_Emm2.m

%%% in this version I will tag the ROIs to make the sorting easier once
%%% they get them back from the ANTS



%%% to get the data (I need to have a path to the folder with the caiman CNMF
%%% results.
MatFiles=dir('*analysis_matlab.mat');


%%% NOTE: the *correlation.mat from slice 49 file of the 22 fish (201810254)
%%% came empty from the cnmf. In this script it is used for the size of the
%%% image so I used another slice.

%%
%%%to get the individual fish names or number


load('s20_good_NumbersNidx_Plane.mat','idx_Plane','Numbers');
 
%%%%% to get the corrected idx_Fish
 
 idx_Fish_name={};
for i=1:length(MatFiles) %%%%to take slices one by one	
    name=strcat(MatFiles(i).name);
    
    [name2,~]=regexp(name,'loomhab_(\d+)_','tokens','match'); %%%to get the number of the fish 
    [name3,~]=regexp(name,'fish(\d+)_','tokens','match'); %%%to get the number of the fish
    
    Fish=strcat(name2{1}{1},name3{1}{1}); %%%to get the number of the fish 
    idx_Fish_name{i}=Fish; %%%to put the number of the plane on the correspondent goodnumbers in the idx_Fish
end
clearvars i Fish Plane name counter %%%to get rid of vairables we will not use anymore


%Fish_list=unique(idx_Fish_name);
Fish_list=unique(idx_Fish);


% To Create ROI csv files

Errored_ROI={};

progressbar(0,0,0);

for fish_nb=1:length(Fish_list) %%%% it will crash at 22 due to one that fish having an error file.
    temp_day=char(Fish_list(fish_nb));
    fish_name=strcat('_loomhab_',temp_day(1:8),'_fish',temp_day(9:end),'_');
    IndexC=strfind({MatFiles.name},fish_name);
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    progressbar(fish_nb/length(Fish_list));
    Centroids=zeros(1,6);    

    for file_nb=1:length(MatFiles_fish)
        progressbar([],file_nb/length(MatFiles_fish));
        name=MatFiles(MatFiles_fish(file_nb)).name;
        Rs=load(name, 'ROIs');
        Rs=Rs.ROIs;
        F=load(name, 'idx_components');%Include all the ROIs
        F=F.idx_components+1;
        Rs=Rs(:,F);
        cor_name=strrep(name,'analysis_matlab','correlation');
        cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
        dims=size(cor_im);
        if file_nb==1
            temp_roi=0;
        else
            temp_roi=temp_roi+size(ROI,3);
        end
        ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
        [slice,~]=regexp(name,'Slice(\d+)_','tokens','match');slice=str2num(slice{1}{1});
        for roi_nb=1:size(ROI,3)
            progressbar([],[],roi_nb/size(ROI,3));
            temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
            Centroids(roi_nb+temp_roi,5)=roi_nb+temp_roi;
            %Centroids(roi_nb+temp_roi,1:2)=temp.Centroid; With the resize
            %need to multiply coords by 1.5
            temp=temp.Centroid;
            Centroids(roi_nb+temp_roi,1:2)=temp;
            Centroids(roi_nb+temp_roi,3)=slice;
            if MatFiles_fish(file_nb)==1
           
            Centroids(roi_nb+temp_roi,6)=roi_nb;
            else
            Centroids(roi_nb+temp_roi,6)=MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+roi_nb;   
            end
        end     
    end    
    if iscell(Fish_list)
        image_name=strcat('_ROIsFish',Fish_list{fish_nb},'_tagged.csv');
    else
        image_name=strcat('_ROIsFish',num2str(Fish_list(fish_nb)),'_tagged.csv');
    end
    dlmwrite(image_name,Centroids,'precision',8); 
end

clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb

%%
 
 %%% for fish 201810254,  the slice49 correlation file is empty so I
 %%% will use a previous slice. this is only used for the size of the image
 %%% so i will use the same one for all slices.
 
 progressbar(0,0,0);

for fish_nb=22%:length(Fish_list)
    temp_day=char(Fish_list(fish_nb));
    fish_name=strcat('_loomhab_',temp_day(1:8),'_fish',temp_day(9:end),'_');
    IndexC=strfind({MatFiles.name},fish_name);
    MatFiles_fish = find(not(cellfun('isempty', IndexC)));
    progressbar(fish_nb/length(Fish_list));
    Centroids=zeros(1,6);    

    name4=MatFiles(MatFiles_fish(1)).name; %%% to take the name of the first slice
    
    for file_nb=1:length(MatFiles_fish)
        progressbar([],file_nb/length(MatFiles_fish));
        name=MatFiles(MatFiles_fish(file_nb)).name;
        Rs=load(name, 'ROIs');
        Rs=Rs.ROIs;
        F=load(name, 'idx_components');%Include all the ROIs
        F=F.idx_components+1;
        Rs=Rs(:,F);
        cor_name=strrep(name4,'analysis_matlab','correlation'); %%% to take always the data of the first slice for the size of the image
        cor_im=load(cor_name);cor_im=cor_im.Correlation_image;
        dims=size(cor_im);
        if file_nb==1
            temp_roi=0;
        else
            temp_roi=temp_roi+size(ROI,3);
        end
        ROI=reshape(full(Rs),dims(1),dims(2),size(Rs,2));counter=1;
        [slice,~]=regexp(name,'Slice(\d+)_','tokens','match');slice=str2num(slice{1}{1});
        for roi_nb=1:size(ROI,3)
            progressbar([],[],roi_nb/size(ROI,3));
            temp=regionprops(uint16(squeeze(ROI(:,:,roi_nb)))==max(max(uint16(squeeze(ROI(:,:,roi_nb))))),'Centroid');
            Centroids(roi_nb+temp_roi,5)=roi_nb+temp_roi;
            %Centroids(roi_nb+temp_roi,1:2)=temp.Centroid; With the resize
            %need to multiply coords by 1.5
            temp=temp.Centroid;
            Centroids(roi_nb+temp_roi,1:2)=temp;
            Centroids(roi_nb+temp_roi,3)=slice;
            if MatFiles_fish(file_nb)==1
           
            Centroids(roi_nb+temp_roi,6)=roi_nb;
            else
            Centroids(roi_nb+temp_roi,6)=MatFiles(MatFiles_fish(file_nb)-1).GoodNumber+roi_nb;   
            end
            
        end     
    end    
    if iscell(Fish_list)
        image_name=strcat('_ROIsFish',Fish_list{fish_nb},'_tagged.csv');
    else
        image_name=strcat('_ROIsFish',num2str(Fish_list(fish_nb)),'_tagged.csv');
    end
    dlmwrite(image_name,Centroids,'precision',8);
end

clearvars slice fish_nb roi_nb ROI Centroids IndexC file_nb



