


%%%% this script is to look at the ROIs sizes to confirm if we are
%%%% getting single cell neurons. 

%%%% Then also taking an example slice and present it. This is for
%%%% supplementary figure 2

%%%% 1 pixel is equal to 0.3226um at binning of one. the images were
%%%% aquired at binning of 4 (4x4). so each pixel has 1.2904um. 
 
%%%% To plot the ROIs in an example plane it is needed to have the temporal
%%%% mean images of each plane (an ouput of Caiman)
%%%% they were added in a folder "means_f20"

%% for f20



load('final_F20_step1.mat','MatFiles_f20','idx_Fish_f20','idx_rsq_test_f20short');

load('f20_cleaned_idxs.mat')

%%%% distinguishable_colors funciton from:
%%%% Tim Holy (2018). Generate maximally perceptually-distinct colors (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors), MATLAB Central File Exchange. Retrieved November 4, 2021.
%%%% https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors

Numbers=[0 [MatFiles_f20.GoodNumber]]; %%%to take make a vector with the GoodNumber field
counter=1;


%% first checking all the ROIs (from CNMF output)


destdirectory = strcat('D:/Emmanuel/matlab_temp_2020/ROIs_plot_test/');
mkdir(destdirectory);   %create the directory
clearvars idx filename


colors = distinguishable_colors(10,[1 1 1; 0 0 0]); %%%here we use a script from Matlab (downloaded, and needs to be in the folder) to generate colors
colors = colors*256;

Area_ROIs=[];

for idx=1:length(MatFiles_f20)
    filename=MatFiles_f20(idx).name;%%%to put the respective name of the files (in this case the slices)
    ROIsNb=[];ClusterNb=[];%%% make the variables we are going to use
    
    
        imagename=regexp(filename,'_output_analysis','split');
        
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*64;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles_f20(idx).ROI;       
      
         %%%%% to try to put random colours (brainbow like)
         ClusterNb=randi(10,size(ROIs,2),1);
         
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;
            image2=(image2/max(max(image2)));%%%%% this is to get only the highest values of the ROI. or the "core"
            image2=uint8(image2);
            
            area_temp=length(find(image2));
            Area_ROIs=vertcat(Area_ROIs,area_temp);
            
            for j=1:3
                
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
      
        name=strcat(destdirectory,'_Test_',imagename(4:end));
   

end
clearvars idx i tempFileNb fileNb AVG_files filename image counter image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster


edges=[0:10:200];
 figure;histogram(Area_ROIs,edges);
 
 std(Area_ROIs)
 mean(Area_ROIs)
 
 %%%% to get the quantiles of the radious of ROIs in um
(sqrt(quantile(Area_ROIs,[0.025 0.25 0.50 0.75 0.975])))*1.2904
 
%%%%% to get the distribution bases on area of a circle
 (sqrt(Area_ROIs)/pi)*1.2904; %%% this would be the radious
 ((sqrt(Area_ROIs)/pi)*1.2904)*2; %%% this would be the diameter
 
 edges=[0:0.5:20];
 figure;histogram(ans,edges);
 
 
 %% to check the filtered brain ROIs (with linear regression)
 %%%%% as some of the CNMF non-filtered ROIs could be outside the brain. 
 
 
Area_ROIs_rsq=[];

for idx=1:length(MatFiles_f20)
    filename=MatFiles_f20(idx).name;%%%to put the respective name of the files (in this case the slices)
    ROIsNb=[];ClusterNb=[];%%% make the variables we are going to use
    
    
    %%%%%this is to locate in every plane 
    
        tempROIsNb=find(idx_rsq_test_f20short_cleaned<=Numbers(idx+1) & idx_rsq_test_f20short_cleaned>Numbers(idx)); 
        if tempROIsNb            
            ROIsNb=[ROIsNb; idx_rsq_test_f20short_cleaned(tempROIsNb)];%%%%            
            
            ClusterNb=[ClusterNb; randi(10,length(tempROIsNb),1)];
        end
    
    
    %%%% and this part is to make and image using the _mean.tif image as
    %%%% base, where we will add the ROIs located before and color them 
    
    
        imagename=regexp(filename,'_output_analysis','split');
        
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*64;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles_f20(idx).ROI;       
         
         ROIsNb=ROIsNb-(Numbers(idx));
         ROIs=ROIs(:,ROIsNb);
         
         %%%%% to try to put random colours (brainbow like)

         
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;
            image2=(image2/max(max(image2)));
            image2=uint8(image2);
            
            area_temp=length(find(image2));
            Area_ROIs_rsq=vertcat(Area_ROIs_rsq,area_temp);
            
            for j=1:3
                
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
      
        name=strcat(destdirectory,'_Test_',imagename(4:end));
   

end
clearvars idx i tempFileNb fileNb AVG_files filename image counter image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster


edges=[0:10:200];
 figure;histogram(Area_ROIs_rsq,edges);
 
 std(Area_ROIs_rsq)
 mean(Area_ROIs_rsq)
 
 %%%% to get the quantiles of the radious of ROIs in um
(sqrt(quantile(Area_ROIs_rsq,[0.025 0.25 0.50 0.75 0.975])))*1.2904

%%%%% to get the distribution bases on area of a circle
 (sqrt(Area_ROIs_rsq)/pi)*1.2904; %%%% for radious
 ((sqrt(Area_ROIs_rsq)/pi)*1.2904)*2; %%%% for diameter
 
 Q_rsq=quantile(ans,[0.025 0.25 0.50 0.75 0.975]);
 
 edges=[0:0.5:20];
 figure;histogram(ans,edges);
 
 
 %%
 
 save('ROIs_size_f20.mat','Area_ROIs','Area_ROIs_rsq');
 
 %% ploting some of the slices of an average fish
 
 figure;histogram(idx_Fish_f20);
 
 cat_fish=categorical(idx_Fish_f20);
 figure;histogram(cat_fish);
 
 counter=1;
 for i=unique(idx_Fish_f20)'
    ROIs_perfish(counter)=length(find(idx_Fish_f20==i));
    counter=counter+1;
 end
 
 mean(ROIs_perfish)
 %%%% 
 
 %%%%% fish 7 and 44 are very close
 
 
destdirectory2 = strcat('D:/Emmanuel/matlab_temp_2020/ROIs_plot_test2/');
mkdir(destdirectory2);   %create the directory
clearvars idx filename
 
 %%%% plotting the all ROIs and after rsq
for idx=1:length(MatFiles_f20)
    filename=MatFiles_f20(idx).name;%%%to put the respective name of the files (in this case the slices)
    
    if ~isempty(regexp(filename,'fish44_'))
    
    ROIsNb=[];ClusterNb=[];%%% make the variables we are going to use
      
    %%%% and this part is to make and image using the _mean.tif image 
        imagename=regexp(filename,'_output_analysis','split');
        
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*180;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles_f20(idx).ROI;       
      
         
         %%%%% to try to put random colours (brainbow like)
         ClusterNb=randi(10,size(ROIs,2),1);
         
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;
            image2=(image2/max(max(image2)));
            image2=uint8(image2);
                                 
            for j=1:3
                
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
      
        name=strcat(destdirectory,'_Test_',imagename(4:end));
    imwrite(image3,name,'tif');
    
   
    
    ROIsNb=[];ClusterNb=[];%%% make the variables we are going to use
      
    
        tempROIsNb=find(idx_rsq_test_f20short_cleaned<=Numbers(idx+1) & idx_rsq_test_f20short_cleaned>Numbers(idx)); %%%to put the data filtered of the selected clusters on a new variable. but i dont understand the sentece... why not just ask for Numbers(idx)?
        if tempROIsNb            
            ROIsNb=[ROIsNb; idx_rsq_test_f20short_cleaned(tempROIsNb)];%%%%            
            
            ClusterNb=[ClusterNb; randi(10,length(tempROIsNb),1)];
        end
    
    
    %%%% and this part is to make and image using the _mean.tif image 
        imagename=regexp(filename,'_output_analysis','split');
        
        imagename=strcat(imagename{1},'_mean.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*180;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=MatFiles_f20(idx).ROI;       
         
         ROIsNb=ROIsNb-(Numbers(idx));
         ROIs=ROIs(:,ROIsNb);
         
         %%%%% to try to put random colours (brainbow like)

         
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;
            image2=(image2/max(max(image2)));
            image2=uint8(image2);
                                 
            for j=1:3
                
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
      
        name=strcat(destdirectory2,'_Test_',imagename(4:end));
    imwrite(image3,name,'tif');
    
    
    
    else
    end

end
clearvars idx i tempFileNb fileNb AVG_files filename image counter image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

