
%%%%% this script is to make the graphs that will go on the paper. 

%%% first the comparison between f20 and f60


load('Nodes_N_means_alldatasets2.mat','Nodes2','ROI_temp2_all','Zbrain_brainMask2D');

load('graph_loomNdataset2.mat','Data_corrMat2');

datasets=['f20'; 'f60'; 's20'; 's60'];
load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);

load('graph_analysis_loomhab3.mat')


cbrewer()

[RdYlBu]=cbrewer('div','RdYlBu',101);
[RdBu]=cbrewer('div','RdBu',101);
[RdBu2]=cbrewer('div','RdBu',1001);
[PRGn]=cbrewer('div','PRGn',101);
[PiYG]=cbrewer('div','PiYG',101);




%%

%%% first, the matrices of the 1st, 10th, 11th and loom similar to 11th. 
%%% with a correlation  
similar_to_rec=struct;
for data=1:4
     datatemp=datasets(data,:);
     
     temp_mat1=Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;

    temp=[];
    temp2=[];
    
for k=1:11
    temp_mat2=Data_corrMat2.(datatemp).Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    temp_corr=corrcoef(temp_mat1,temp_mat2);
    temp(k)=temp_corr(1,2);
    temp_mean=mean(mean(temp_mat2));
    temp2(k)=temp_mean;
end
similar_to_rec.(datatemp).corr=temp;
similar_to_rec.(datatemp).Mean_corr=temp2;
end

 
for data=1:4
     datatemp=datasets(data,:);
    
    [M,I] = max(similar_to_rec.(datatemp).corr);
    %similar_to_rec.(datatemp).corr
    I-1 %%% i take one out cause the first one is the pre loom moment. so this result is the actual loom
end


%% for f20 and f60

%%%% to plot the 1st, 4th, 10th and 11th
counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 data=1;
     datatemp=datasets(data,:);
subplot(2,4,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for first loom
subplot(2,4,counter+1);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,4,counter+2);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
subplot(2,4,counter+3);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 4th loom
title(datatemp);
data=2;
     datatemp=datasets(data,:);
subplot(2,4,counter+4);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for first loom
subplot(2,4,counter+5);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,5}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,4,counter+6);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,11}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
subplot(2,4,counter+7);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,12}(keep,keep));pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for the 6th loom
title(datatemp);
%counter=counter+8;

%%%% ploting the color of the nodes in the axes of the matrices

cluster_colormap=zeros(length(Nodes2.Mod_clust(keep)),3);
for i=1:length(Nodes2.Mod_clust(keep))

    if Nodes2.Mod_clust(keep(i))== 2 | Nodes2.Mod_clust(keep(i))== 4 | Nodes2.Mod_clust(keep(i))== 5
    cluster_colormap(i,:)=[0 1 0];
    elseif Nodes2.Mod_clust(keep(i))== 6
    cluster_colormap(i,:)=[0 0 1];
    elseif Nodes2.Mod_clust(keep(i))== 1
    cluster_colormap(i,:)=[1 0 0];
    elseif Nodes2.Mod_clust(keep(i))== 7
    cluster_colormap(i,:)=[1 0 1];
    end
end


counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 data=1;
     datatemp=datasets(data,:);
subplot(2,4,counter);imagesc(Data_corrMat2.(datatemp).Mean_corrMat{1,2}(keep,keep)); pbaspect([1 1 1]);caxis([-1 1]); colormap(cluster_colormap);colorbar;%% for first loom
title(datatemp);


%%% and now with the flitered matrix above 0.75

MatAll_corrected2=struct;
for data=1:4
    datatemp=datasets(data,:);
    moment=[2 3 4 5 6 7 8 9 10 11 12]; %%% looms 1,2,3,10 and 11 (cause 1 is pre loom)
    loom=[1 2 3 4 5 6 7 8 9 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(Data_corrMat2.(datatemp).Mean_corrMat{1,moment(m)}(keep,keep)),0.75);
     MatAll_corrected2.(datatemp).(strcat('loom',num2str(loom(m)))).Mat=Mat;
     end
end


%%%%%
counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 data=1;
     datatemp=datasets(data,:);
     loom=fieldnames(MatAll_corrected2.(datatemp));
subplot(2,4,counter);imagesc(MatAll_corrected2.(datatemp).(loom{1}).Mat); pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar;%% for pre loom
subplot(2,4,counter+1);imagesc(MatAll_corrected2.(datatemp).(loom{10}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar;%% for 10th loom
subplot(2,4,counter+2);imagesc(MatAll_corrected2.(datatemp).(loom{11}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar; %% for 11th loom
subplot(2,4,counter+3);imagesc(MatAll_corrected2.(datatemp).(loom{4}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar; %% for the 4th loom
title(datatemp);
data=2;
     datatemp=datasets(data,:);
     loom=fieldnames(MatAll_corrected2.(datatemp));
subplot(2,4,counter+4);imagesc(MatAll_corrected2.(datatemp).(loom{1}).Mat); pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar;%% for pre loom
subplot(2,4,counter+5);imagesc(MatAll_corrected2.(datatemp).(loom{10}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar;%% for 10th loom
subplot(2,4,counter+6);imagesc(MatAll_corrected2.(datatemp).(loom{11}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar; %% for 11th loom
subplot(2,4,counter+7);imagesc(MatAll_corrected2.(datatemp).(loom{6}).Mat);pbaspect([1 1 1]);caxis([0.7 1]); colormap(inferno);colorbar; %% for the 6th loom
title(datatemp);


%% making some matrices substractions

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 data=1;
     datatemp=datasets(data,:);
     loom=fieldnames(MatAll_corrected2.(datatemp));
subplot(2,3,counter);imagesc(MatAll_corrected2.(datatemp).(loom{1}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{1}).Mat); pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for pre loom
subplot(2,3,counter+1);imagesc(MatAll_corrected2.(datatemp).(loom{10}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{10}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar;%% for 10th loom
subplot(2,3,counter+2);imagesc(MatAll_corrected2.(datatemp).(loom{11}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{11}).Mat);pbaspect([1 1 1]);caxis([-1 1]); colormap(RdBu);colorbar; %% for 11th loom
title(datatemp);


%%%% to put the empty places with black boxes. 

sub1=MatAll_corrected2.(datatemp).(loom{1}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{1}).Mat;
sub10=MatAll_corrected2.(datatemp).(loom{10}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{10}).Mat;
sub11=MatAll_corrected2.(datatemp).(loom{11}).Mat-MatAll_corrected2.(datasets(data+1,:)).(loom{11}).Mat;

max(max(sub1))
min(min(sub11))

sub1(find(sub1==0)) = -1;
sub10(find(sub10==0)) = -1;
sub11(find(sub11==0)) = -1;

ddd=[0 0 0;RdBu];

counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
 data=1;
     datatemp=datasets(data,:);
     loom=fieldnames(MatAll_corrected2.(datatemp));
subplot(2,3,counter);imagesc(sub1); pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar;%% for pre loom
subplot(2,3,counter+1);imagesc(sub10);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar;%% for 10th loom
subplot(2,3,counter+2);imagesc(sub11);pbaspect([1 1 1]);caxis([-1 1]); colormap(ddd);colorbar; %% for 11th loom
title(datatemp);


%% substraction of connections
%%%%  f20 minus f60 first. 


%%
%%%% here i am trying to get the connections from the previous filtering
  
for k=[1:31]
  
    temp_mat1=Data_corrMat2.f20.Mean_corrMat{1,k}(keep,keep);
    temp_mat1(isnan(temp_mat1))=0;
    
    temp_mat1_idx=find(threshold_absolute(abs(temp_mat1),0.75));
    
    temp_mat2=Data_corrMat2.f60.Mean_corrMat{1,k}(keep,keep);
    temp_mat2(isnan(temp_mat2))=0;
    
    temp_mat2_idx=find(threshold_absolute(abs(temp_mat2),0.75));
    
    temp_mat=temp_mat1-temp_mat2;
    
    temp_mat_idx=union(temp_mat1_idx,temp_mat2_idx);
    
    temp_mat_cleaned=zeros(size(temp_mat));
    temp_mat_cleaned(temp_mat_idx)=temp_mat(temp_mat_idx);
    
    f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,k}=temp_mat_cleaned;
    
end


    %% now to plot some graphs. 
    

 count=1;
    for k=[2 11 12]
    figure;
    set(gcf, 'Position',  [200, 200, 700, 900]);
    
    
    R=f20_f60_subs_mat.Subs_Mean_corrMat_cleaned{1,k}; 
   
    
n=length(Nodes2.Mod_loc(keep));

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
p.LineWidth = abs(G.Edges.Weight)*3;
%p.LineWidth = 1;
%axis off

p.NodeLabel=[];
p.EdgeAlpha=0.75;
% get the grid coordinates for all nodes
%[x,y] = ndgrid(1:ceil(sqrt(n)),1:ceil(sqrt(n)));
x = Nodes2.Mod_loc(keep,1);
y = Nodes2.Mod_loc(keep,2);
% set the nodes in a 'grid' structure
p.XData = x(1:n);
p.YData = y(1:n);
%axis ij % flip the plot so it will be orderd like in a matrix

hold on;
%%%% Ploting the nodes and the brain layout
%gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust,'rgggbm','.',20,'off');
gscatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),Nodes2.Mod_clust(keep),'rgggbm','.',15,'off');


view(-90,90)
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'w');xlim([300 1350])
set(gca,'Color','k');%%% to make the background black
hold off;

%saveas(gcf,strcat('subsf20vsf60_',num2str(k),'.svg'));

count=count+1;

    end


%% Participation f20 vs f60


%%% Participation with merged stronghab clusters

unique(Nodes2.Mod_clust)
Nodes2.Mod_clust_Fhab_pooled=Nodes2.Mod_clust;
Nodes2.Mod_clust_Fhab_pooled(find(Nodes2.Mod_clust == 4 | Nodes2.Mod_clust == 5))=2;

unique(Nodes2.Mod_clust_Fhab_pooled)


for data=1:4
datatemp=datasets(data,:);
loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)

temp_mat=MatAll_corrected2.(datatemp).(loom{i}).Mat;
temp_mat(isnan(temp_mat))=0;
P=participation_coef(temp_mat,Nodes2.Mod_clust_Fhab_pooled(keep));
[Gpos,~]=gateway_coef_sign(MatAll_corrected2.(datatemp).(loom{i}).Mat,Nodes2.Mod_clust_Fhab_pooled(keep),1);

MatAll_corrected2.(datatemp).(loom{i}).P2=P;
MatAll_corrected2.(datatemp).(loom{i}).Gpos2=Gpos;
end
end



%%% in the brain
figure;
counter=1;

for data=1:2
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected2.(datatemp));

for i=[1 2 3 4 5 10 11]%1:length(loom)
    
    low=0;high=0.8;
     
  subplot(2,7,counter);
  plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k'); xlim([300 1400]);
  hold on; scatter(Nodes2.Mod_loc(keep,1),Nodes2.Mod_loc(keep,2),25,MatAll_corrected2.(datatemp).(loom{i}).P2,'filled');colormap(inferno);view(-90,90);caxis([low high]);%colorbar; %
 
counter=counter+1;
end
end

%%% raster
figure;
counter=1;
for data=1:2
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected2.(datatemp));
temp=[];
for i=1:length(loom)

    temp(:,i)=(MatAll_corrected2.(datatemp).(loom{i}).P2)';
    

end 

subplot(1,2,counter);
imagesc(temp);caxis([0 0.8]);colormap(inferno); colorbar;
counter=counter+1;
end

%%% raster substraction
figure;
counter=1;
subs_Par_perLoom_f20_f60=[];
for data=1
    datatemp=datasets(data,:);
    loom=fieldnames(MatAll_corrected2.(datatemp));

for i=1:length(loom)

    subs_Par_perLoom_f20_f60(:,i)=(MatAll_corrected2.(datatemp).(loom{i}).P2-MatAll_corrected2.(datasets(data+1,:)).(loom{i}).P2)';
    

end 

imagesc(subs_Par_perLoom_f20_f60);caxis([-0.8 0.8]);colormap(RdBu); colorbar;
counter=counter+1;
end


