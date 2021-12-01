
%%%%% This script is to get the means and rasterplots for figure 2

%%%% for f20, which is the one shown in figure 2.

load('final_f20_step1.mat','ZS_f20','rawregressf20','idx_rsq_test_f20short','High_corr_Nb_f20','High_corr_Nb_f20_short','gooodmaps');


%%
%%% to get the means of the 4 main clusters and the 3 kinds of fast hab
%%% clusters for figures. 


%%% first for CL4
for j=gooodmaps
idx_temp=find(High_corr_Nb_f20_short==j);

if j==2
mean_CL4_f20.fasthab=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));
elseif j==6 
mean_CL4_f20.slopehab=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));
elseif j==1
mean_CL4_f20.nonhab=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));   
else
mean_CL4_f20.inhib=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));
end

end

%figure;plot(mean_CL4_f20.slopehab);
%figure;plot(mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:)));

for j=1:size(rawregressf20,1)
idx_temp=find(High_corr_Nb_f20==j);

if j==2
mean_CL7_f20.fasthab_med=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));
mean_CL7_f20_short.fasthab_med=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),55:180));
elseif j==3
mean_CL7_f20.sound=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));  
mean_CL7_f20_short.sound=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),55:180));  
elseif j==4
mean_CL7_f20.fasthab_sharp=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));
mean_CL7_f20_short.fasthab_sharp=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),55:180));
elseif j==5
mean_CL7_f20.fasthab_broad=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));
mean_CL7_f20_short.fasthab_broad=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),55:180));
elseif j==6 
mean_CL7_f20.slopehab=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));
mean_CL7_f20_short.slopehab=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),55:180));
elseif j==1
mean_CL7_f20.nonhab=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:)); 
mean_CL7_f20_short.nonhab=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),55:180)); 
else
mean_CL7_f20.inhib=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),:));
mean_CL7_f20_short.inhib=mean(ZS_f20(idx_rsq_test_f20short(idx_temp),55:180));
end

end

%figure;plot(mean_CL7_f20.fasthab_sharp);
%figure;plot(mean_CL7_f20_short.slopehab);

save('means_f20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20','mean_CL7_f20_short','-v7.3');



%%% now getting the raster plots of the 4 main clusters

for j=gooodmaps
idx_temp=find(High_corr_Nb_f20_short==j);

if j==2
figure;imagesc(ZS_f20(idx_rsq_test_f20short(idx_temp),:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_CL4_fasthab.tif');

elseif j==6 
 figure;imagesc(ZS_f20(idx_rsq_test_f20short(idx_temp),:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_CL4_slopehab.tif');   

elseif j==1
  figure;imagesc(ZS_f20(idx_rsq_test_f20short(idx_temp),:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot
%saveas(gcf,'rasterplot_f20_CL4_nonhab.tif');  

else
  figure;imagesc(ZS_f20(idx_rsq_test_f20short(idx_temp),:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot; 
%saveas(gcf,'rasterplot_f20_CL4_inhib.tif');    
    
figure;imagesc(ZS_f20(idx_rsq_test_f20short(idx_temp),:),[min(min(mean_CL7_f20.inhib))-0.5 max(max(mean_CL7_f20.fasthab_sharp))]);colormap hot; colorbar
%saveas(gcf,'rasterplot_f20_CL4_inhib_wbar.tif');   

end

end


%%% to checkit worked
names=fieldnames(mean_CL4_f20);
counter=1;
figure;
for i=1:length(names)-1
subplot(4,2,counter);plot(mean_CL4_f20.(names{i}));

counter=counter+1;
end


%%% to checkit worked
names=fieldnames(mean_CL4_f20);

figure;
for i=1:length(names)-1
plot(mean_CL4_f20.(names{i}));
hold on;

end


%% example also with s20 dataset 

load('final_S20_step1.mat','ZS_S20','rawregressS20','idx_rsq_test_S20short','High_corr_Nb_S20','High_corr_Nb_S20_short','gooodmaps');


%%
%%% to get the means of the 4 main clusters and the 3 kinds of fast hab
%%% clusters for figures. 


%%% first for CL4
for j=gooodmaps
idx_temp=find(High_corr_Nb_S20_short==j);

if j==1
mean_CL4_s20.fasthab=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));
elseif j==5 
mean_CL4_s20.slopehab=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));
elseif j==3
mean_CL4_s20.nonhab=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));   
else
mean_CL4_s20.inhib=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));
end

end

%figure;plot(mean_CL4_S20.inhib);
% j=4
% idx_temp=find(High_corr_Nb_S20_short==j);
%figure;plot(mean(ZS_S20(idx_rsq_test_S20short(idx_temp),:)));

for j=1:size(rawregressS20,1)
idx_temp=find(High_corr_Nb_S20==j);

if j==4
mean_CL7_s20.fasthab_med=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));
mean_CL7_s20_short.fasthab_med=mean_CL7_s20.fasthab_med(1,55:180);
elseif j==6
mean_CL7_s20.sound=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));  
mean_CL7_s20_short.sound=mean_CL7_s20.sound(1,55:180);  
elseif j==2
mean_CL7_s20.fasthab_sharp=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));
mean_CL7_s20_short.fasthab_sharp=mean_CL7_s20.fasthab_sharp(1,55:180);
elseif j==1
mean_CL7_s20.fasthab_broad=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));
mean_CL7_s20_short.fasthab_broad=mean_CL7_s20.fasthab_broad(1,55:180);
elseif j==5
mean_CL7_s20.slopehab=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));
mean_CL7_s20_short.slopehab=mean_CL7_s20.slopehab(1,55:180);
elseif j==3
mean_CL7_s20.nonhab=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20)); 
mean_CL7_s20_short.nonhab=mean_CL7_s20.nonhab(1,55:180);
else
mean_CL7_s20.inhib=mean(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20));
mean_CL7_s20_short.inhib=mean_CL7_s20.inhib(1,55:180);
end

end

%figure;plot(mean_CL7_S20.nonhab);
%figure;plot(mean_CL7_S20_short.slopehab);

save('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20','mean_CL7_s20_short','-v7.3');



%%% now getting the raster plots of the 4 main clusters

for j=gooodmaps
idx_temp=find(High_corr_Nb_S20_short==j);

if j==1
figure;imagesc(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20),[min(min(mean_CL7_s20.inhib))-0.5 max(max(mean_CL7_s20.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_S20_CL4_fasthab.tif');

elseif j==5 
 figure;imagesc(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20),[min(min(mean_CL7_s20.inhib))-0.5 max(max(mean_CL7_s20.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_S20_CL4_slopehab.tif');   

elseif j==3
  figure;imagesc(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20),[min(min(mean_CL7_s20.inhib))-0.5 max(max(mean_CL7_s20.fasthab_sharp))]);colormap hot
saveas(gcf,'rasterplot_S20_CL4_nonhab.tif');  

else
  figure;imagesc(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20),[min(min(mean_CL7_s20.inhib))-0.5 max(max(mean_CL7_s20.fasthab_sharp))]);colormap hot; 
saveas(gcf,'rasterplot_S20_CL4_inhib.tif');    
    
figure;imagesc(ZS_S20(idx_rsq_test_S20short(idx_temp),ZS_short_S20),[min(min(mean_CL7_s20.inhib))-0.5 max(max(mean_CL7_s20.fasthab_sharp))]);colormap hot; colorbar
saveas(gcf,'rasterplot_S20_CL4_inhib_wbar.tif');   

end

end

