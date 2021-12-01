
%%%%% this script is to make figures of the sensitivity test to the nomber
%%%%% of nodes (nodes_sensitivity_analysis_f20_x2n4.m) and the per fish
%%%%% graph analysis (perfish_graph_analysis_f20.m).

%%%% i will need to load the data and maybe other things that I will need

cbrewer()

[RdBu]=cbrewer('div','RdBu',101);

%%%%% for the matrices matrices. 
%%%%% with the f20, ~200, ~400 and a single fish matrices

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
for data=1%:4
     datatemp=datasets(data,:);
subplot(4,6,counter);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,1}(OriginalData_corrMat.keep,OriginalData_corrMat.keep));  pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for pre loom
subplot(4,6,counter+1);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,2}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 1st loom
subplot(4,6,counter+2);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,3}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,6,counter+3);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,4}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(4,8,counter+4);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,5}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(4,8,counter+5);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,6}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,6,counter+4);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,11}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(4,6,counter+5);imagesc(OriginalData_corrMat.Data_corrMat2.(datatemp).Mean_corrMat{1,12}(OriginalData_corrMat.keep,OriginalData_corrMat.keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title(datatemp);

counter=counter+6;
end

%counter=1;
for data=1%:4
     datatemp=datasets(data,:);
subplot(4,6,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,1}(keep2,keep2));  pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for pre loom
subplot(4,6,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep2,keep2)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 1st loom
subplot(4,6,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,3}(keep2,keep2)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,6,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,4}(keep2,keep2)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(4,8,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep2,keep2)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(4,8,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,6}(keep2,keep2)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,6,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep2,keep2)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(4,6,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep2,keep2)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title(strcat(datatemp,'/',num2str(length(keep2))));

counter=counter+6;
end

%counter=1;
for data=1%:4
     datatemp=datasets(data,:);
subplot(4,6,counter);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,1}(keep4,keep4));  pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for pre loom
subplot(4,6,counter+1);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,2}(keep4,keep4)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 1st loom
subplot(4,6,counter+2);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,3}(keep4,keep4)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,6,counter+3);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,4}(keep4,keep4)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(4,8,counter+4);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,5}(keep4,keep4)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(4,8,counter+5);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,6}(keep4,keep4)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,6,counter+4);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,11}(keep4,keep4)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(4,6,counter+5);imagesc(Data_corrMat4.(datatemp).Mean_corrMat{1,12}(keep4,keep4)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title(strcat(datatemp,'/',num2str(length(keep4))));

counter=counter+6;
end

%counter=1;
 tempfish=2;
for data=1%:4
    datatemp=datasets(data,:);
       
subplot(4,6,counter);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,1});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for pre loom
subplot(4,6,counter+1);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,2});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 1st loom
subplot(4,6,counter+2);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,3});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,6,counter+3);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,4});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(4,8,counter+4);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,5});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(4,8,counter+5);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,6});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(4,6,counter+4);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,11});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(4,6,counter+5);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,12});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title(fish{tempfish});

counter=counter+6;
end

saveas(gcf,'SensNsinglefish_Matrices_f20.svg');


%% density plot for # of nodes comparison 


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp);ylim([0 0.8]);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected2.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected4.(datatemp).(loom{i}).kden;
end
plot(temp);
hold on;

end
title('density');
legend('99 nodes (f20)','197 nodes','368 nodes');
hold off;

saveas(gcf,'density_f20_NodeSensitivity.svg');

%% participation plot for # of nodes comparison



figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).P);
end
plot(temp);ylim([0 0.8]);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected2.(datatemp).(loom{i}).P);
end
plot(temp);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected4.(datatemp).(loom{i}).P);
end
plot(temp);
hold on;

end
title('Participation');
legend('99 nodes (f20)','197 nodes','368 nodes');
hold off;

saveas(gcf,'Participation_f20_NodeSensitivity.svg');


%% density for single fish


figure;
temp=[];
temp_mean=[];
for data=1%:4
    datatemp=datasets(data,:);
    
    for tempfish=1:length(fish)
        
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp(1,i)=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).kden;
        end
        
        
        temp_mean(tempfish,:)=temp;
        
        if tempfish==2
        plot(temp,'r','LineWidth',3);ylim([0 0.8]);
        else
        plot(temp);ylim([0 0.8]);
        end
        
        hold on;
        
    end
    
    
end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp,'k','LineWidth',3);
hold on;

end
title('Density per fish');
legend('fish_1','fish_7','fish_13','fish_17','fish_21','fish_25','fish_29','fish_34','fish_40','fish_44','fish_48','f20');
hold off;

saveas(gcf,'density_f20_perFish.svg');

%% participation per single fish


figure;
temp_mean=[];
for data=1%:4
    datatemp=datasets(data,:);
    
    
    for tempfish=1:length(fish)
        temp=[];
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp(:,i)=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).Ccoef;
        end
        
        temp_mean(tempfish,:)=nanmean(temp);
        
               
        if tempfish==2
        plot(nanmean(temp),'r','LineWidth',3);ylim([0 0.8]);
        else
        plot(nanmean(temp));ylim([0 0.8]);
        end
        
        hold on;
        
    end
    
    
end

temp=[];
for data=1%:4
    datatemp=datasets(data,:);
    loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));
    
    for i=1:length(loom)
        
        temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).Ccoef;
    end
    plot(nanmean(temp),'k','LineWidth',3);
    hold on;
        
end
title('Participation per fish');
legend('fish_1','fish_7','fish_13','fish_17','fish_21','fish_25','fish_29','fish_34','fish_40','fish_44','fish_48','f20');
hold off;

saveas(gcf,'participation_f20_perFish.svg');


%% merging node sensitivity and single fish plots

%%%%% density


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp,'LineWidth',3);ylim([0 0.8]);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected2.(datatemp).(loom{i}).kden;
end
plot(temp,'LineWidth',3);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=MatAll_corrected4.(datatemp).(loom{i}).kden;
end
plot(temp,'LineWidth',3);
hold on;

end

temp=[];
temp_mean=[];
for data=1%:4
    datatemp=datasets(data,:);
    
    for tempfish=1:length(fish)
        
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp(1,i)=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).kden;
        end
        
        
        temp_mean(tempfish,:)=temp;
        
        if tempfish==2
        plot(temp,'r','LineWidth',3);ylim([0 0.8]);
        else
        plot(temp);ylim([0 0.8]);
        end
        
        hold on;
        
    end
    
    
end

temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp,'k','LineWidth',3);
hold on;

end
title('density');
legend('99 nodes (f20)','197 nodes','368 nodes','fish_1','fish_7','fish_13','fish_17','fish_21','fish_25','fish_29','fish_34','fish_40','fish_44','fish_48','f20');
hold off;

saveas(gcf,'density_f20_nodesNperFish.svg');



%%%% participation


figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).P);
end
plot(temp,'LineWidth',3);ylim([0 0.8]);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected2.(datatemp).(loom{i}).P);
end
plot(temp,'LineWidth',3);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected4.(datatemp).(loom{i}).P);
end
plot(temp,'LineWidth',3);
hold on;

end

temp_mean=[];
for data=1%:4
    datatemp=datasets(data,:);
    
    
    for tempfish=1:length(fish)
        temp=[];
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp(:,i)=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).P;
        end
        
        temp_mean(tempfish,:)=nanmean(temp);
        
               
        if tempfish==2
        plot(nanmean(temp),'r','LineWidth',3);ylim([0 0.8]);
        else
        plot(nanmean(temp));ylim([0 0.8]);
        end
        
        hold on;
        
    end
    
    
end

temp=[];
for data=1%:4
    datatemp=datasets(data,:);
    loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));
    
    for i=1:length(loom)
        
        temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).P;
    end
    plot(nanmean(temp),'k','LineWidth',3);
    hold on;
        
end
title('Participation');
legend('99 nodes (f20)','197 nodes','368 nodes','fish_1','fish_7','fish_13','fish_17','fish_21','fish_25','fish_29','fish_34','fish_40','fish_44','fish_48','f20');
hold off;

saveas(gcf,'participation_f20_nodesNperFish.svg');


%%%% Cluster coeficient

figure;
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).Ccoef);
end
plot(temp,'LineWidth',3);ylim([0 0.8]);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected2.(datatemp).(loom{i}).Ccoef);
end
plot(temp,'LineWidth',3);
hold on;

end
temp=[];
for data=1%:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected4.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(MatAll_corrected4.(datatemp).(loom{i}).Ccoef);
end
plot(temp,'LineWidth',3);
hold on;

end

temp_mean=[];
for data=1%:4
    datatemp=datasets(data,:);
    
    
    for tempfish=1:length(fish)
        temp=[];
        loom=fieldnames(MatAll_corrected_singleFish.(datatemp).(fish{tempfish}));
        
        for i=1:length(loom)
            
            temp(:,i)=MatAll_corrected_singleFish.(datatemp).(fish{tempfish}).(loom{i}).Ccoef;
        end
        
        temp_mean(tempfish,:)=nanmean(temp);
        
               
        if tempfish==2
        plot(nanmean(temp),'r','LineWidth',3);ylim([0 0.8]);
        else
        plot(nanmean(temp));ylim([0 0.8]);
        end
        
        hold on;
        
    end
    
    
end

temp=[];
for data=1%:4
    datatemp=datasets(data,:);
    loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));
    
    for i=1:length(loom)
        
        temp(:,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).Ccoef;
    end
    plot(nanmean(temp),'k','LineWidth',3);
    hold on;
        
end
title('Cluster coeficient');
legend('99 nodes (f20)','197 nodes','368 nodes','fish_1','fish_7','fish_13','fish_17','fish_21','fish_25','fish_29','fish_34','fish_40','fish_44','fish_48','f20');
hold off;

saveas(gcf,'Ccoef_f20_nodesNperFish.svg');


%% just checking participation/density ratio


figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=mean(OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).P);
end
plot(temp,'LineWidth',3);%ylim([0 0.8]);
hold on;

end

figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden;
end
plot(temp,'LineWidth',3);%ylim([0 0.8]);
hold on;

end


figure;
temp=[];
for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(OrginalMatAll_corrected.MatAll_corrected.(datatemp));

for i=1:length(loom)
    
    temp(1,i)=OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).kden/mean(OrginalMatAll_corrected.MatAll_corrected.(datatemp).(loom{i}).P);
end
plot(temp,'LineWidth',3);%ylim([0 0.8]);
hold on;

end


%% single fish correlation matrix

%%% in case we need it


counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
for tempfish=1:length(fish)


for data=1%:4
    datatemp=datasets(data,:);
       
subplot(6,11,counter);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,1});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for pre loom
subplot(6,11,counter+11);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,2});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 1st loom
subplot(6,11,counter+22);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,3});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(6,11,counter+33);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,4});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(11,8,counter+4);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,5});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
%subplot(11,8,counter+5);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,6});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);
subplot(6,11,counter+44);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,11});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);%% for 10th loom
subplot(6,11,counter+55);imagesc(Data_corrMat_singleFish.(datatemp).(fish{tempfish}).loomsR{1,12});pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu); %% for 11th loom
title(fish{tempfish});

counter=counter+1;
end
end

%%%%saveas(gcf,'SingleFishMatrices.svg');%%%%% i get erros strying to save
%%%%it as an SVG file
saveas(gcf,'SingleFishMatrices.png');