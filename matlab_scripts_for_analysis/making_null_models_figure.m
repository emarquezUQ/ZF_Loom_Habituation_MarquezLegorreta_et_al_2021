


%%%% this script is to generate a null models of using the WT nodes from the
%%%% f20 experiment and to make figures for the paper. 


%%% i am using the aaft function from : https://au.mathworks.com/matlabcentral/fileexchange/16062-test-of-non-linearity
%%% which is used in this paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1369647
%%% and in turn bases its aaft in this one: D. Kugiumtzis, Surrogate data test for nonlinearity including monotonic transformations, Phys Rev E, vol. 62, 1, 2000


load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

datasets=['f20'; 'f60'; 's20'; 's60'];


cbrewer()

[RdBu]=cbrewer('div','RdBu',101);


%% doing the null graph
%%% I am first taking the timelines of the nodes to use them as sample to
%%% generate a surrogates with aaft
%%% as some nodes are empty I wil do an average of all of them 

%wtFish=fieldnames(NodesFmr1.wt.NodeMats); %%% to grab the WT from the fmr1
%dataset

f20Fish=fieldnames(Nodes.f20.NodeMats);

avgWTmat=[];
for i=1:length(f20Fish)
  
  tempMat=Nodes.f20.NodeMats.(f20Fish{i});  
  avgWTmat=cat(3,avgWTmat,tempMat);  
    
end

avgWTmat=nanmean(avgWTmat,3);

avgWTmat=avgWTmat(keep,:);

figure;imagesc(avgWTmat);

figure; imagesc(avgWTmat,[-2 11]);colormap hot; %colorbar;
%saveas(gcf,'rasterplot_avgf20nodes.tif');
saveas(gcf,'rasterplot_avgf20nodes.svg');
saveas(gcf,'rasterplot_avgf20nodes_wcolorbar.svg');

%%
%%% surrogate timeline for each node

very_good_surrMat=[];
for i=1:length(NodesFmr1.Nod_clustID)

    tempMean=avgWTmat(i,:); 
   %figure;plot(tempMean)
   temp_surrMat=aaft(tempMean,1);
   
   very_good_surrMat=vertcat(very_good_surrMat,temp_surrMat');
   
end


figure;imagesc(very_good_surrMat);

figure; imagesc(very_good_surrMat,[-2 11]);colormap hot; %colorbar;
%saveas(gcf,'rasterplot_aaft_f20nodes.tif');
saveas(gcf,'rasterplot_aaft_f20nodes.svg');


%%
%%%%% finally, I will try doing a surrogate for each loom

%%% I need to do it per loom. as previously
%%%from the wildtype loomhab dataset
load('All_means_maxResp_CL4_normalized.mat','Loomf20_onset', 'loom_moments', 'Loomf20_onset_idx','S_trim');

super_good_surrMat=avgWTmat;

 for node=1:size(avgWTmat,1)              
        
     temp_mean=avgWTmat(node,:);
     
        for k=1:31
        
        if k==1
            temp_mean(:,10:45)=aaft(temp_mean(:,10:45),1);
        elseif k==11||k==21 ||k==31
            temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28)= aaft(temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28),1);
        else
            
            temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36)= aaft(temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36),1); 
                           
        end
        end
        
        super_good_surrMat(node,:)=temp_mean;
 end

figure;imagesc(super_good_surrMat);

figure; imagesc(super_good_surrMat,[-2 11]);colormap hot; %colorbar;
%saveas(gcf,'rasterplot_aaft_looms_f20nodes.tif');
saveas(gcf,'rasterplot_aaft_looms_f20nodes.svg');


%% making the graph

%%% just testing
%R_temp=corrcoef(good_surrMat');


%%%% for the average model
Data_corrMat_superNull=struct;

temp_mean=avgWTmat;
                
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat_superNull.loomsR=temp_R;
    


%%%% for the aaft model of each node
Data_corrMat_Null=struct;

temp_mean=very_good_surrMat;
                
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat_Null.loomsR=temp_R;
    
    
    
%%%% for the aaft model of each loom
Data_corrMat_Null2=struct;

temp_mean=super_good_surrMat;
                
        temp_R={};
        
        for k=1:21
        
        if k==1
            temp_R{1,k}=corrcoef((temp_mean(:,10:45))');
        elseif k==11||k==21 
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+28))');
        else
            
            temp_R{1,k}= corrcoef((temp_mean(:,Loomf20_onset_idx(k-1):Loomf20_onset_idx(k-1)+36))'); 
                           
        end
        end
        
    Data_corrMat_Null2.loomsR=temp_R;
        
%% doing graph analysis

%%%% plotting the matrices of the node average model
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
counter=1;
subplot(4,8,counter);imagesc(Data_corrMat_superNull.loomsR{1,1});  pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat_superNull.loomsR{1,2}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat_superNull.loomsR{1,3}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+3);imagesc(Data_corrMat_superNull.loomsR{1,4}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+4);imagesc(Data_corrMat_superNull.loomsR{1,5}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+5);imagesc(Data_corrMat_superNull.loomsR{1,6}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+6);imagesc(Data_corrMat_superNull.loomsR{1,11}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat_superNull.loomsR{1,12}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title('node average model');


 %%%% plotting the matrices of the null model (aaft of each node)
 counter=counter+8;
%figure;
%counter=1;
subplot(4,8,counter);imagesc(Data_corrMat_Null.loomsR{1,1});  pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat_Null.loomsR{1,2}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat_Null.loomsR{1,3}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+3);imagesc(Data_corrMat_Null.loomsR{1,4}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+4);imagesc(Data_corrMat_Null.loomsR{1,5}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+5);imagesc(Data_corrMat_Null.loomsR{1,6}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+6);imagesc(Data_corrMat_Null.loomsR{1,11}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat_Null.loomsR{1,12}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title('null model aaft 1');


 %%%% plotting the matrices of the null model (aaft of each loom)
 counter=counter+8;
%figure;
%counter=1;
subplot(4,8,counter);imagesc(Data_corrMat_Null2.loomsR{1,1});  pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat_Null2.loomsR{1,2}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat_Null2.loomsR{1,3}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+3);imagesc(Data_corrMat_Null2.loomsR{1,4}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+4);imagesc(Data_corrMat_Null2.loomsR{1,5}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+5);imagesc(Data_corrMat_Null2.loomsR{1,6}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+6);imagesc(Data_corrMat_Null2.loomsR{1,11}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat_Null2.loomsR{1,12}); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title('null model aaft 2');

%%%% plotting the matrices i used in the paper for f20
counter=counter+8;
%counter=1;
%figure;
data=1; %%%% for f20
     datatemp=datasets(data,:);
subplot(4,8,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1}(keep,keep));  pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for pre loom
subplot(4,8,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 1st loom
subplot(4,8,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,8,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(4,8,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title(datatemp);

saveas(gcf,'nullMatrices_f20nodes.svg');


%% BCT analysis
%%%%% now further analyze the null models graph using
%%%%% the Brain Connectivity Toolbox
%%%%% I will look at density and recovery of each null model and compare it
%%%%% to the FnS datasets

%% for node average model
MatAll_corrected_avg=struct;

    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat_superNull.loomsR{1,moment(m)}),0.75);
     MatAll_corrected_avg.(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end

     
%%% density
loom=fieldnames(MatAll_corrected_avg);

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected_avg.(loom{i}).Mat);

MatAll_corrected_avg.(loom{i}).kden=temp_kden;
end

%%% participation and gateway coef

loom=fieldnames(MatAll_corrected_avg);

for i=1:length(loom)

temp_mat=MatAll_corrected_avg.(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Nod_clustID(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected_avg.(loom{i}).Mat,Nodes.Nod_clustID(keep),1);

MatAll_corrected_avg.(loom{i}).P=P;
MatAll_corrected_avg.(loom{i}).Gpos=Gpos;
end


%% for aaft node model

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_null=struct;

    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat_Null.loomsR{1,moment(m)}),0.75);
     MatAll_corrected_null.(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end

%%% density

loom=fieldnames(MatAll_corrected_null);

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected_null.(loom{i}).Mat);

MatAll_corrected_null.(loom{i}).kden=temp_kden;
end


%%% participation and gateway coef

loom=fieldnames(MatAll_corrected_null);

for i=1:length(loom)

temp_mat=MatAll_corrected_null.(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Nod_clustID(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected_null.(loom{i}).Mat,Nodes.Nod_clustID(keep),1);

MatAll_corrected_null.(loom{i}).P=P;
MatAll_corrected_null.(loom{i}).Gpos=Gpos;
end

%% having a look at the aaft null model 2
%%%% surrogate timelines just at the loom timepoints

MatAll_corrected_null2=struct;

    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat_Null2.loomsR{1,moment(m)}),0.75);
     MatAll_corrected_null2.(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end


%%% density
loom=fieldnames(MatAll_corrected_null2);

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected_null2.(loom{i}).Mat);

MatAll_corrected_null2.(loom{i}).kden=temp_kden;
end


%%%% participation and gateway coef

loom=fieldnames(MatAll_corrected_null2);

for i=1:length(loom)

temp_mat=MatAll_corrected_null2.(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Nod_clustID(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected_null2.(loom{i}).Mat,Nodes.Nod_clustID(keep),1);

MatAll_corrected_null2.(loom{i}).P=P;
MatAll_corrected_null2.(loom{i}).Gpos=Gpos;
end


%% checking the FnS datasets to compare
load('graphs_FnS_all.mat');

MatAll_corrected=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end


%%% density

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected.(datatemp).(loom{i}).Mat);

MatAll_corrected.(datatemp).(loom{i}).kden=temp_kden;
end
end


%%% participation and gateway coef


for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_mat=MatAll_corrected.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes.Nod_clustID(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected.(datatemp).(loom{i}).Mat,Nodes.Nod_clustID(keep),1);

MatAll_corrected.(datatemp).(loom{i}).P=P;
MatAll_corrected.(datatemp).(loom{i}).Gpos=Gpos;
end
end


%% making figures

%%%%% density plot 
figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end

%%% avg model
temp=[];
for i=1:length(loom)   
    temp(1,i)=MatAll_corrected_avg.(loom{i}).kden;
end
plot(temp);
hold on;

%%% null 1
temp=[];
for i=1:length(loom)   
    temp(1,i)=MatAll_corrected_null.(loom{i}).kden;
end
plot(temp);
hold on;

%%% null 3
temp=[];
for i=1:length(loom)   
    temp(1,i)=MatAll_corrected_null2.(loom{i}).kden;
end
plot(temp);
legend('f20','f60','s20','s60','node average model','null model aaft 1','null model aaft 2');
hold off;

saveas(gcf,'density_allDataNnulls.svg');

%%%%% average participation

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected.(datatemp).(loom{i}).P);
end
plot(temp);
hold on;

end

%%% avg model
temp=[];
for i=1:length(loom)   
    temp(1,i)=mean(MatAll_corrected_avg.(loom{i}).P);
end
plot(temp);
hold on;

%%% null 1
temp=[];
for i=1:length(loom)   
    temp(1,i)=mean(MatAll_corrected_null.(loom{i}).P);
end
plot(temp);
hold on;

%%% null 3
temp=[];
for i=1:length(loom)   
    temp(1,i)=mean(MatAll_corrected_null2.(loom{i}).P);
end
plot(temp);
legend('f20','f60','s20','s60','node average model','null model aaft 1','null model aaft 2');
hold off;

saveas(gcf,'participation_allDataNnulls.svg');

%%%%% participation brain plots

figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]);

low=0;high=0.8;
   
counter=1;
for i=1:length(loom)
    
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_avg.(loom{i}).P,'filled');colormap(inferno);caxis([low high]);view(-90,90);%colorbar; 
 title(strcat('avg model','/',(loom{i})));
counter=counter+1;
end
for i=1:length(loom)
         
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_null.(loom{i}).P,'filled');colormap(inferno);caxis([low high]);view(-90,90);%colorbar; 
 title(strcat('null model 1','/',(loom{i})));
counter=counter+1;
end
for i=1:length(loom)   
       
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected_null2.(loom{i}).P,'filled');colormap(inferno);caxis([low high]);view(-90,90);%colorbar; 
 title(strcat('null model 2','/',(loom{i})));
counter=counter+1;
end
for i=1:length(loom)
       
  subplot(4,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes.Nod_coor(keep,1),Nodes.Nod_coor(keep,2),25,MatAll_corrected.f20.(loom{i}).P,'filled');colormap(inferno);caxis([low high]);view(-90,90);%colorbar; 
  title(strcat('f20','/',(loom{i})));  
counter=counter+1;
end

saveas(gcf,'participations_allDataNnulls_brains.svg');

%% cross validation figures


for data=1%:4

names = fieldnames(Data_corrMat2.(datasets(data,:)));
CorrMatrices_mean2=zeros(31,length(names)-1,99,99);
for loom=1:31        
    for fish_rem_nb=1:length(names)-1
        temp=nan(length(names)-2,99,99);    
        counter=1;
        for fish_nb=1:length(names)-1
            if fish_nb ~= fish_rem_nb
                fish_name=names(fish_nb);
                temp(counter,:,:)=Data_corrMat2.(datasets(data,:)).(fish_name{1}).loomsR{1,loom}(keep,keep);
                counter=counter+1;            
            end
        end
        CorrMatrices_mean2(loom,fish_rem_nb,:,:)=squeeze(nanmean(temp,1));
    end       
end

end



%%%% 5 random group of matrices without a fish at random

figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
counter=1;
for fish=randperm(length(names)-1,5)
for loom=[1 2 3 4 5 6 11 12]

    subplot(5,8,counter);imagesc(squeeze(squeeze(CorrMatrices_mean2(loom,fish,:,:)))); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
    title(strcat('minus_',names(fish)));
    counter=counter+1;
end
end
saveas(gcf,'crossVal_5randomMat.svg');

%% to look at the changes in density

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end

%% density
%


%kden = density_und(CIJ);

for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)

temp_kden=density_und(MatAll_corrected.(datatemp).(loom{i}).Mat);

MatAll_corrected.(datatemp).(loom{i}).kden=temp_kden;
end
end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end



%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_crossval=[];

for fish=1:length(names)-1

    moment=[1 2 3 4 5 6 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    %loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(squeeze(squeeze(CorrMatrices_mean2(moment(m),fish,:,:)))),0.75);
     MatAll_corrected_crossval(m,fish,:,:)=Mat;
     end

end

%%%% calculating the density


crossval_density=[];
for fish=1:length(names)-1
for m=1:length(moment)

temp_kden=density_und(squeeze(squeeze(MatAll_corrected_crossval(m,fish,:,:))));

crossval_density(m,fish)=temp_kden;
end
end


figure;
temp2=[];
for fish=1:length(names)-1
for m=1:length(moment)-1    
%scatter(m,crossval_density(m+1,fish),'filled');

temp2(1,m)=crossval_density(m+1,fish);
plot(temp2,'r','LineWidth',1);ylim([0 0.8]);
hold on;
end
end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp,'k','LineWidth',2);ylim([0 0.8]);
hold on;

end
hold off;

saveas(gcf,'crossVal_f20_density.svg');

% for fish=1:length(names)-1
% for m=1:length(moment)-1    
% scatter(m,crossval_density(m+1,fish),'filled');
% hold on;
% end
% end
% hold off;

