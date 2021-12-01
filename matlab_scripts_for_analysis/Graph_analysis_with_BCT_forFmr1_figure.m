
%%%%% this script is to make the graphs that will go on the paper. 

%%% this is for the fmr1 dataset whle using the nodes and cluster from the wild type dataset
load('NodesNgraphFmr1Loomhab4.mat','Nodes4','goodorder_clust','Data_corrMat4','keep','discard');


%load('Nodes_N_means_alldatasets2.mat','Zbrain_brainMask2D');

%save('Zbrain_brainMask2D.mat','Zbrain_brainMask2D');
load('Zbrain_brainMask2D.mat','Zbrain_brainMask2D');


%%% making colormaps
cbrewer()

[RdYlBu]=cbrewer('div','RdYlBu',101);
[RdBu]=cbrewer('div','RdBu',101);
[RdBu2]=cbrewer('div','RdBu',1001);
[PRGn]=cbrewer('div','PRGn',101);
[PiYG]=cbrewer('div','PiYG',101);

ddd=[0 0 0;RdBu];%%% for black at the end


%%


groupnames=fieldnames(Nodes4.ROIs_idx);


%% for controls vs fmr1

%% first the matrices
%%%%
counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %% controls
     group=groupnames{g,1};
subplot(2,5,counter);imagesc(Data_corrMat4.(group).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for first loom
subplot(2,5,counter+1);imagesc(Data_corrMat4.(group).Mean_corrMat{1,3}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 2nd loom
subplot(2,5,counter+2);imagesc(Data_corrMat4.(group).Mean_corrMat{1,4}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 3rd loom
subplot(2,5,counter+3);imagesc(Data_corrMat4.(group).Mean_corrMat{1,11}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 10th loom
subplot(2,5,counter+4);imagesc(Data_corrMat4.(group).Mean_corrMat{1,12}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar; %% for 11th loom
%subplot(2,4,counter+3);imagesc(Data_corrMat4.(group).Mean_corrMat{1,3}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar; %% for the 2nd loom
title(group);
g=2; %% fmr1
     group=groupnames{g,1};
subplot(2,5,counter+5);imagesc(Data_corrMat4.(group).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for first loom
subplot(2,5,counter+6);imagesc(Data_corrMat4.(group).Mean_corrMat{1,3}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 2nd loom
subplot(2,5,counter+7);imagesc(Data_corrMat4.(group).Mean_corrMat{1,4}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 3rd loom
subplot(2,5,counter+8);imagesc(Data_corrMat4.(group).Mean_corrMat{1,11}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 10th loom
subplot(2,5,counter+9);imagesc(Data_corrMat4.(group).Mean_corrMat{1,12}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar; %% for 11th loom
%subplot(2,4,counter+7);imagesc(Data_corrMat4.(group).Mean_corrMat{1,3}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar; %% for the 2nd loom
title(group);
%counter=counter+8;

%%%% ploting the color of the nodes in the axes of the matrices
cluster_colormap=zeros(length(Nodes4.Mod_clust(keep)),3);
for i=1:length(Nodes4.Mod_clust(keep))

    if Nodes4.Mod_clust(keep(i))== 1 | Nodes4.Mod_clust(keep(i))== 2 | Nodes4.Mod_clust(keep(i))== 3
    cluster_colormap(i,:)=[0 1 0];
    elseif Nodes4.Mod_clust(keep(i))== 4
    cluster_colormap(i,:)=[0 0 1];
    elseif Nodes4.Mod_clust(keep(i))== 5
    cluster_colormap(i,:)=[1 0 0];
    elseif Nodes4.Mod_clust(keep(i))== 6
    cluster_colormap(i,:)=[1 0 1];
    end
end


counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=1;
    group=groupnames{g,1};
subplot(2,4,counter);imagesc(Data_corrMat4.(group).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(cluster_colormap);colorbar;%% for first loom
title(group);




%%% and now with the flitered matrix above 0.75

MatAll_corrected2=struct;
for g=1:3
    group=groupnames{g,1};
    moment=[2 3 4 5 6 7 8 9 10 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 6 7 8 9 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat4.(group).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected2.(group).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end


%% making some matrices substraction of substractions

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %% controls
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,3,counter);imagesc(MatAll_corrected2.(group).(loom{1}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{1}).Mat); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for pre loom
subplot(2,3,counter+1);imagesc(MatAll_corrected2.(group).(loom{10}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{10}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,3,counter+2);imagesc(MatAll_corrected2.(group).(loom{11}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{11}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
title(group);


%%%% to try to put the empty places with black boxes. 

%%% substractions of ctrl-fmr1
sub1=MatAll_corrected2.(group).(loom{1}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{1}).Mat;
sub2=MatAll_corrected2.(group).(loom{2}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{2}).Mat;
sub3=MatAll_corrected2.(group).(loom{3}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{3}).Mat;
sub4=MatAll_corrected2.(group).(loom{4}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{4}).Mat;
sub5=MatAll_corrected2.(group).(loom{5}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{5}).Mat;
sub10=MatAll_corrected2.(group).(loom{10}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{10}).Mat;
sub11=MatAll_corrected2.(group).(loom{11}).Mat-MatAll_corrected2.(groupnames{g-1,1}).(loom{11}).Mat;

max(max(sub1))
min(min(sub11))

sub1(find(sub1==0)) = -1;
sub2(find(sub2==0)) = -1;
sub3(find(sub3==0)) = -1;
sub4(find(sub4==0)) = -1;
sub5(find(sub5==0)) = -1;
sub10(find(sub10==0)) = -1;
sub11(find(sub11==0)) = -1;

ddd=[0 0 0;RdBu];

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %%controls
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,3,counter);imagesc(sub1); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar;%% for pre loom
subplot(2,3,counter+1);imagesc(sub10);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar;%% for 10th loom
subplot(2,3,counter+2);imagesc(sub11);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar; %% for 11th loom
title('control minus fmr1');




%%%% these are the matrices used to make the graphs in figure 5c
counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %%controls
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,5,counter);imagesc(sub1); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 1st loom
subplot(2,5,counter+1);imagesc(sub2); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 2nd loom
subplot(2,5,counter+2);imagesc(sub3); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 3rd loom
subplot(2,5,counter+3);imagesc(sub10);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 10th loom
subplot(2,5,counter+4);imagesc(sub11);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar; %% for 11th loom
title('control minus fmr1');


%%
%%% substractions of ctrl-hets

sub1_hets=MatAll_corrected2.(group).(loom{1}).Mat-MatAll_corrected2.(groupnames{g-2,1}).(loom{1}).Mat;
sub2_hets=MatAll_corrected2.(group).(loom{2}).Mat-MatAll_corrected2.(groupnames{g-2,1}).(loom{2}).Mat;
sub3_hets=MatAll_corrected2.(group).(loom{3}).Mat-MatAll_corrected2.(groupnames{g-2,1}).(loom{3}).Mat;
sub4_hets=MatAll_corrected2.(group).(loom{4}).Mat-MatAll_corrected2.(groupnames{g-2,1}).(loom{4}).Mat;
sub5_hets=MatAll_corrected2.(group).(loom{5}).Mat-MatAll_corrected2.(groupnames{g-2,1}).(loom{5}).Mat;
sub10_hets=MatAll_corrected2.(group).(loom{10}).Mat-MatAll_corrected2.(groupnames{g-2,1}).(loom{10}).Mat;
sub11_hets=MatAll_corrected2.(group).(loom{11}).Mat-MatAll_corrected2.(groupnames{g-2,1}).(loom{11}).Mat;

max(max(sub1_hets))
min(min(sub11_hets))

sub1_hets(find(sub1_hets==0)) = -1;
sub2_hets(find(sub2_hets==0)) = -1;
sub3_hets(find(sub3_hets==0)) = -1;
sub4_hets(find(sub4_hets==0)) = -1;
sub5_hets(find(sub5_hets==0)) = -1;
sub10_hets(find(sub10_hets==0)) = -1;
sub11_hets(find(sub11_hets==0)) = -1;


counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %%controls
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,5,counter);imagesc(sub1_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 1st loom
subplot(2,5,counter+1);imagesc(sub2_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 2nd loom
subplot(2,5,counter+2);imagesc(sub3_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 3rd loom
subplot(2,5,counter+3);imagesc(sub10_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 10th loom
subplot(2,5,counter+4);imagesc(sub11_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar; %% for 11th loom
title('control minus hets');

%%% ploting substraction of ctrls vs fmr1 and hets

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 g=3; %%controls
     group=groupnames{g,1};
     loom=fieldnames(MatAll_corrected2.(group));
subplot(2,5,counter);imagesc(sub1); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 1st loom
subplot(2,5,counter+1);imagesc(sub2); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 2nd loom
subplot(2,5,counter+2);imagesc(sub3); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 3rd loom
subplot(2,5,counter+3);imagesc(sub10);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 10th loom
subplot(2,5,counter+4);imagesc(sub11);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar; %% for 11th loom
title('control minus fmr1');

subplot(2,5,counter+5);imagesc(sub1_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 1st loom
subplot(2,5,counter+6);imagesc(sub2_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 2nd loom
subplot(2,5,counter+7);imagesc(sub3_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 3rd loom
subplot(2,5,counter+8);imagesc(sub10_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar;%% for 10th loom
subplot(2,5,counter+9);imagesc(sub11_hets); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);%colorbar; %% for 11th loom
title('control minus hets');

%% substraction of connections
%%%% I am trying control minus fmr1. 


%%
%%%%  getting the connections from the previous filtering
  
for k=[1:21]
  
    temp_mat1=Data_corrMat4.control.Mean_corrMat{1,k}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat1_idx=find(threshold_absolute(abs(temp_mat1),0.75));
    
    temp_mat2=Data_corrMat4.fmr1.Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat2_idx=find(threshold_absolute(abs(temp_mat2),0.75));
    
    temp_mat=temp_mat1-temp_mat2;
    
    temp_mat_idx=union(temp_mat1_idx,temp_mat2_idx);
    
    temp_mat_cleaned=zeros(size(temp_mat));
    temp_mat_cleaned(temp_mat_idx)=temp_mat(temp_mat_idx);
    
    ctrl_fmr1_subs_mat.Subs_Mean_corrMat_cleaned{1,k}=temp_mat_cleaned;
    
end

    
    %% now to plot some graphs. 
    
%     figure;
%      count=1;
%     for g=1:3
%  group=groupnames{g,1};
 %figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%sgtitle(group);
 count=1;
    for k=[2 3 4 11 12]
    figure;
    set(gcf, 'Position',  [200, 200, 700, 900]);
    %set(gcf, 'Position',  [200, 200, 1200, 900]);
    
    
    R=ctrl_fmr1_subs_mat.Subs_Mean_corrMat_cleaned{1,k}; 
    
    
n=length(Nodes4.Mod_loc(keep));

% set the source of the lines:
s = repelem(1:n-1,n-1:-1:1);

% set the target of the lines:
t = nonzeros(triu(repmat(2:n,n-1,1)).').';

%[~,~,weights] = find(tril(R,-1));
weights = nonzeros(tril(R,-1));

%quantile(abs(weights),[0.025 0.25 0.50 0.75 0.975])


% create the graph object:
%G = graph(s,t,weights,n);
G = graph(R);

% mark the lines to remove from the graph:
%threshold = 0; %  minimum correlation to plot
%threshold = 0; %  minimum correlation to plot
%line_to_remove = isnan(weights) | abs(weights)<threshold;
% remove the lines from the graph:
%G = G.rmedge(find(line_to_remove)); %#ok<FNDSB>


%subplot(1,3,count);
% plot it:
p = plot(G); % for labeling the lines uncomment add: 'EdgeLabel',G.Edges.Weight
p.NodeColor = 'k';
% color positive in blue and negative in red:
  colormap(RdBu) ;caxis([-1 1]);colorbar;
  %colormap jet;caxis([0 3]);%colorbar;
  p.EdgeCData=G.Edges.Weight;
%p.EdgeColor = [G.Edges.Weight>0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight<0.']; %% red high, blue low
%p.EdgeColor = [G.Edges.Weight<0.' zeros(numel(G.Edges.Weight),1) G.Edges.Weight>0.'];

% set the thickness of the lines:
p.LineWidth = abs(G.Edges.Weight)*2;
%p.LineWidth = 1;
%axis off

p.NodeLabel=[];
p.EdgeAlpha=0.75;
% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes4.Mod_loc(keep,1);
y = Nodes4.Mod_loc(keep,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
%gscatter(Nodes4.Mod_loc(keep,1),Nodes4.Mod_loc(keep,2),Nodes4.Mod_clust,'rgggbm','.',20,'off');
gscatter(Nodes4.Mod_loc(keep,1),Nodes4.Mod_loc(keep,2),Nodes4.Mod_clust(keep),'gggbrm','.',15,'off');


view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'w');xlim([300 1350])
set(gca,'Color','k');%%% to make the background black
hold off;

%saveas(gcf,strcat('subsf20vsf60_',num2str(k),'.svg'));

count=count+1;

    end

    %end

%% now degrees, strenth and participation

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll_corrected2.(group));

for i=1:length(loom)

deg=degrees_und(MatAll_corrected2.(group).(loom{i}).Mat);
str=strengths_und_sign(abs(MatAll_corrected2.(group).(loom{i}).Mat));

MatAll_corrected2.(group).(loom{i}).deg=deg;
MatAll_corrected2.(group).(loom{i}).str=str;
end
end


%%%% to plot with brains
figure;
counter=1;

for g=[3 1 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=[1:5 10 11]
    if counter==1|counter==8|counter==15|counter==22
    low=0;high=100;
    else
    low=0;high=50; 
    end
    
  subplot(3,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes4.Mod_loc(keep,1),Nodes4.Mod_loc(keep,2),25,MatAll_corrected2.(group).(loom{i}).deg,'filled');colormap(inferno);caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end


figure;
counter=1;

for g=3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=[1:5 10 11]
    
    low=-50;high=50;
    
subplot(1,7,counter);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]); 
hold on; scatter(Nodes4.Mod_loc(keep,1),Nodes4.Mod_loc(keep,2),25,(MatAll_corrected2.(group).(loom{i}).deg-MatAll_corrected2.(groupnames{g-1,1}).(loom{i}).deg),'filled');colormap(RdBu);caxis([low high]);view(-90,90);%colorbar; 
 
counter=counter+1;
end
end

%%%% degrees as rasterplot

figure;
counter=1;
for g=[3 1 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected2.(group).(loom{i}).str)';
    

end 

subplot(1,3,counter);
imagesc(temp);caxis([0 25]);colormap(inferno); colorbar;
counter=counter+1;
end

%%% raster substraction
figure;
counter=1;
for g=3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected2.(group).(loom{i}).deg-MatAll_corrected2.(groupnames{g-1,1}).(loom{i}).deg)';
    

end 

imagesc(temp);caxis([-25 25]);colormap(RdBu); colorbar;
counter=counter+1;
end


%%% participation merging the fasthab cluster subtypes

unique(Nodes4.Mod_clust)
Nodes4.Mod_clust_Fhab_pooled=Nodes4.Mod_clust;
Nodes4.Mod_clust_Fhab_pooled(find(Nodes4.Mod_clust == 2 | Nodes4.Mod_clust == 3))=1;

unique(Nodes4.Mod_clust_Fhab_pooled)

for g=1:3
group=groupnames{g,1};
loom=fieldnames(MatAll_corrected2.(group));

for i=1:length(loom)

temp_mat=MatAll_corrected2.(group).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes4.Mod_clust_Fhab_pooled(keep));
%[Gpos,~]=gateway_coef_sign(MatAll_corrected2.(group).(loom{i}).Mat,Nodes4.Mod_clust_Fhab_pooled(keep),1);

MatAll_corrected2.(group).(loom{i}).P2=P;
%MatAll_corrected2.(group).(loom{i}).Gpos2=Gpos;
end
end


%%%%% this is to make the brain heatmaps in figure 5d
figure;
counter=1;

for g=[3 1 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=[1:5 10 11]
    
    low=0;high=0.8;
     
  subplot(3,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes4.Mod_loc(keep,1),Nodes4.Mod_loc(keep,2),25,MatAll_corrected2.(group).(loom{i}).P2,'filled');colormap(inferno);view(-90,90);caxis([low high]);%colorbar; %
 title(group);
counter=counter+1;
end
end



%%% raster

figure;
counter=1;
for g=[3 1 2]
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected2.(group).(loom{i}).P2)';
    

end 

subplot(1,3,counter);
imagesc(temp);caxis([0 0.8]);colormap(inferno); colorbar;
counter=counter+1;
end

%%% raster substraction (for figure 5e)
figure;
counter=1;
subs_Par_perLoom_ctrl_fmr1=[];
for g=3
    group=groupnames{g,1};
    loom=fieldnames(MatAll_corrected2.(group));

for i=1:length(loom)

    subs_Par_perLoom_ctrl_fmr1(:,i)=(MatAll_corrected2.(group).(loom{i}).P2-MatAll_corrected2.(groupnames{g-1,1}).(loom{i}).P2)';
    

end 

imagesc(subs_Par_perLoom_ctrl_fmr1);caxis([-0.8 0.8]);colormap(RdBu); colorbar;
counter=counter+1;
end


%%


