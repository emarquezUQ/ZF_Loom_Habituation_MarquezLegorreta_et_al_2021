
%%%%% this script is to make the circular graphs of the fmr1 correlation
%%%%% matrices.
%%%%% it's based on the functon from paul-kassebaum-mathworks-circularGraph
%%%%% but I made some edits (in folder: paul-kassebaum-mathworks-circularGraph-3a7926b_EML_edits)

%%%% Because of the changing of the original circulargraph script to add other
%%%% parameters like the edges colors. so some parts wont work anymore as
%%%% before

%%%% this is a continuation of Graph_analysis_with_BCT_forFmr1_figure.m 

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

groupnames=fieldnames(Nodes4.ROIs_idx);

%%

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


%%

g=3;
    group=groupnames{g,1};

loom=fieldnames(MatAll_corrected2.(group));
ctrl_11=MatAll_corrected2.(group).(loom{11}).Mat;
ctrl_11(find(isnan(ctrl_11)))= 0;

fmr1_11=MatAll_corrected2.(groupnames{g-1,1}).(loom{11}).Mat;
fmr1_11(find(isnan(fmr1_11)))= 0;



[~,~,bin] = histcounts(nonzeros(ctrl_11),20);
[Blues]=cbrewer('seq','Blues',20);
ctrl_11_blues=Blues(bin,:);

[BluesDark]=cbrewer('seq','Blues',40);
BluesDark=BluesDark(21:40,:);
ctrl_11_blues_dark=BluesDark(bin,:);

[~,~,bin] = histcounts(nonzeros(fmr1_11),20);
[Reds]=cbrewer('seq','Reds',20);
fmr1_11_reds=Reds(bin,:);

[RedsDark]=cbrewer('seq','Reds',40);
RedsDark=RedsDark(21:40,:);
fmr1_11_reds_dark=RedsDark(bin,:);



%%% creating a functional cluster colormap
funct_clust_colors=zeros(length(Nodes4.Mod_clust(keep)),3);
for i=unique(Nodes4.Mod_clust(keep))'
    
    temp_node=find(Nodes4.Mod_clust(keep)==i);
    
    if i<4    
    funct_clust_colors(temp_node,:)= repmat([0 1 0],size(temp_node,1),1);  
    elseif i==4      
     funct_clust_colors(temp_node,:)= repmat([0 0 1],size(temp_node,1),1);          
    elseif i==5
     funct_clust_colors(temp_node,:)= repmat([1 0 0],size(temp_node,1),1);      
    elseif i==6
     funct_clust_colors(temp_node,:)= repmat([1 0 1],size(temp_node,1),1);                     
    end   
    
end


black_nodes=zeros(length(Nodes4.Mod_clust(keep)),3);

figure;h=circularGraph(fmr1_11,'ColorMapNode',funct_clust_colors,'ColorMapEdges',fmr1_11_reds_dark);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');


%%%% making brain labels
RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};

RegionList_av={'Pal','Sp','Th','Hb','Pt','Tec','Tg','Cb','HB'};


Brain_label={};
for i=unique(Nodes4.Mod_brain(keep))'
    
    temp_node=find(Nodes4.Mod_brain(keep)==i);
    
     for j=temp_node'
    Brain_label{j}= RegionList{i}(1:4); 
     end
    
end

Brain_label2={};
for i=unique(Nodes4.Mod_brain(keep))'
    
    temp_node=find(Nodes4.Mod_brain(keep)==i);
    
     for j=temp_node'
    Brain_label2{j}= RegionList_av{i}; 
     end
    
end

figure;h=circularGraph(ctrl_11,'ColorMapNode',black_nodes,'ColorMapEdges',ctrl_11_blues_dark,'Label',Brain_label2);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');
colormap(BluesDark);
cbh=colorbar;
set(cbh,'XTickLabel',{'0.75','0.8','0.85','0.9','0.95','1'})


figure;h=circularGraph(fmr1_11,'ColorMapNode',black_nodes,'ColorMapEdges',fmr1_11_reds_dark,'Label',Brain_label2);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position;
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');
colormap(RedsDark);
cbh=colorbar;
set(cbh,'XTickLabel',{'0.75','0.8','0.85','0.9','0.95','1'})



%%
%%% now to order by brain regions. 
%%%% first I need to make it symetrical
%%% finding the midline. I should be around 310.
figure;
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([300 1350])
hold on;
gscatter(Nodes4.Mod_loc(keep,1),Nodes4.Mod_loc(keep,2),Nodes4.Mod_brain(keep));
view(-90,90)
plot([500:1200],(ones(length(500:1200))*310),'k');


%%% finding left and right nodes. 

leftH_nodes=find(Nodes4.Mod_loc(keep,2)>310);
rightH_nodes=find(Nodes4.Mod_loc(keep,2)<310);

%%% to see how many there are per brain region in each side and
%%% make it equal
for i=1:length(leftH_nodes)
leftH_nodes(i,2)=Nodes4.Mod_brain(keep(leftH_nodes(i)));
leftH_nodes(i,3)=Nodes4.Mod_clust(keep(leftH_nodes(i)));
end

for i=1:length(rightH_nodes)
rightH_nodes(i,2)=Nodes4.Mod_brain(keep(rightH_nodes(i)));
rightH_nodes(i,3)=Nodes4.Mod_clust(keep(rightH_nodes(i)));
end

Nb_nodes_perBrain=[];
for i=unique(Nodes4.Mod_brain(keep))'

    temp=length(find(leftH_nodes(:,2)==i));
    if isempty(temp)
        temp=0;
    end
    Nb_nodes_perBrain(1,i)=temp;
    
    
     temp=length(find(rightH_nodes(:,2)==i));
    if isempty(temp)
        temp=0;
    end
    Nb_nodes_perBrain(2,i)=temp;
    
end

 nodes_needed=Nb_nodes_perBrain(1,:)- Nb_nodes_perBrain(2,:);
 
 %%%  re-ordering them. left goes ascending and right goes descending
length(leftH_nodes)
[~,I]=sort(leftH_nodes(:,2));
leftH_nodes=leftH_nodes(I,:);


 length(rightH_nodes)
 [~,I]=sort(rightH_nodes(:,2),'descend');
 rightH_nodes=rightH_nodes(I,:);
 
 %%% adding empty nodes to the specific brain regions
new_rightH_nodes=[];
for i=sort(unique(Nodes4.Mod_brain(keep)),'descend')'
    
    add=nodes_needed(i);
    temp1=rightH_nodes(find(rightH_nodes(:,2)==i),:);  
           
    new_rightH_nodes=vertcat(new_rightH_nodes,temp1);
    
    temp2=repmat([0 i 10],add,1);

    new_rightH_nodes=vertcat(new_rightH_nodes,temp2);
    
end



new_brain_order=vertcat(leftH_nodes,new_rightH_nodes);

big_ctrl_11=zeros(length(new_brain_order),length(new_brain_order));
for i=1:length(new_brain_order)
    
   for j=1:length(new_brain_order) 
    
       if new_brain_order(i,1)==0 | new_brain_order(j,1)==0
       
           big_ctrl_11(i,j)=0;
       else
           
           big_ctrl_11(i,j)=ctrl_11(new_brain_order(i,1),new_brain_order(j,1));
           
       end
end

end


big_fmr1_11=zeros(length(new_brain_order),length(new_brain_order));
for i=1:length(new_brain_order)
    
   for j=1:length(new_brain_order) 
    
       if new_brain_order(i,1)==0 | new_brain_order(j,1)==0
       
           big_fmr1_11(i,j)=0;
       else
           
           big_fmr1_11(i,j)=fmr1_11(new_brain_order(i,1),new_brain_order(j,1));
           
       end
end

end


black_nodes_big=zeros(length(new_brain_order),3); %%% for the nodes edges

%%%% making brain labels
Brain_label_big={};
for i=unique(new_brain_order(:,2))'
    
    temp_node=find(new_brain_order(:,2)==i);
    
     for j=temp_node'
    Brain_label_big{j}= RegionList{i}(1:4); 
     end
    
end

%%% creating a functional cluster colormap
funct_clust_colors_big=zeros(length(new_brain_order),3);
for i=unique(new_brain_order(:,3))'
    
    temp_node=find(new_brain_order(:,3)==i);
    
    if i<4    
    funct_clust_colors_big(temp_node,:)= repmat([0 1 0],size(temp_node,1),1);  
    elseif i==4      
     funct_clust_colors_big(temp_node,:)= repmat([0 0 1],size(temp_node,1),1);          
    elseif i==5
     funct_clust_colors_big(temp_node,:)= repmat([1 0 0],size(temp_node,1),1);      
    elseif i==6
     funct_clust_colors_big(temp_node,:)= repmat([1 0 1],size(temp_node,1),1);  
    elseif i==10
     funct_clust_colors_big(temp_node,:)= repmat([0 0 0],size(temp_node,1),1);   
    end   
    
end

%%%% for the edges colormap
[D,edges,bin] = histcounts(nonzeros(big_ctrl_11),20);
[BluesDark]=cbrewer('seq','Blues',40);
BluesDark=BluesDark(21:40,:);
ctrl_11_blues_dark2=BluesDark(bin,:);

[~,~,bin] = histcounts(nonzeros(big_fmr1_11),20);
[RedsDark]=cbrewer('seq','Reds',40);
RedsDark=RedsDark(21:40,:);
fmr1_11_reds_dark2=RedsDark(bin,:);



figure;h=circularGraph(big_ctrl_11,'ColorMapNode',black_nodes_big,'ColorMapEdges',ctrl_11_blues_dark2,'Label',Brain_label_big);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position;
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors_big,'filled');



figure;h=circularGraph(big_fmr1_11,'ColorMapNode',black_nodes_big,'ColorMapEdges',fmr1_11_reds_dark2,'Label',Brain_label_big);

positions=[];
for i=1:length(h.Node)
positions(i,:)=h.Node(i).Position;
end

hold on;
scatter(positions(:,1),positions(:,2),30,funct_clust_colors_big,'filled');

%% making the circular graph figures brain location for looms 1-3,10 and 11


edges2=0.75:0.0125:1;
for L=[1 2 3 10 11]
    
    %%% for ctrls
   for g= [3 2] 
    group=groupnames{g,1}; 
    loom=fieldnames(MatAll_corrected2.(group));
    
    temp_mat=MatAll_corrected2.(group).(loom{L}).Mat;
    temp_mat(find(isnan(temp_mat)))= 0;
    
    big_temp_mat=zeros(length(new_brain_order),length(new_brain_order));
    for ii=1:length(new_brain_order)    
        for j=1:length(new_brain_order)     
            if new_brain_order(ii,1)==0 | new_brain_order(j,1)==0     
           big_temp_mat(ii,j)=0;
            else          
           big_temp_mat(ii,j)=temp_mat(new_brain_order(ii,1),new_brain_order(j,1));         
            end
        end
    end
    
    [~,~,bin] = histcounts(nonzeros(big_temp_mat),edges2);
    
    if g==3
    [TempColor]=cbrewer('seq','Blues',40);
    else
     [TempColor]=cbrewer('seq','Reds',40);   
    end
    
    TempColor=TempColor(21:40,:);
    TempColor_dark=[];
    TempColor_dark=TempColor(bin,:);
    
    %figure('Renderer', 'painters', 'Position', [100 100 650 650]);
    figure('Position', [100 100 650 650]);
    h=circularGraph(big_temp_mat,'ColorMapNode',black_nodes_big,'ColorMapEdges',TempColor_dark,'Label',Brain_label_big);

    positions=[];
        for p=1:length(h.Node)
        positions(p,:)=h.Node(p).Position;
        end

    hold on;
    scatter(positions(:,1),positions(:,2),30,funct_clust_colors_big,'filled');
    hold off;
    
     %set(gcf, 'Renderer','painters');
    %saveas(gcf,strcat('circularGraph_brain_',group,'_',loom{L},'.svg'));
    %saveas(gcf,strcat('circularGraph_brain_',group,'_',loom{L},'.emf'));
    
   end
end


%%

%%% for the graphs based on functional clustering 


edges2=0.75:0.0125:1;
for L=[2 3 10 11]%[1 2 3 10 11]
    
    %%% for ctrls
   for g= [3 2] 
    group=groupnames{g,1}; 
    loom=fieldnames(MatAll_corrected2.(group));
    
    temp_mat=MatAll_corrected2.(group).(loom{L}).Mat;
    temp_mat(find(isnan(temp_mat)))= 0;
         
    [~,~,bin] = histcounts(nonzeros(temp_mat),edges2);
    
    if g==3
    [TempColor]=cbrewer('seq','Blues',40);
    else
     [TempColor]=cbrewer('seq','Reds',40);   
    end
    
    TempColor=TempColor(21:40,:);
    TempColor_dark=[];
    TempColor_dark=TempColor(bin,:);
    
    figure('Renderer', 'painters', 'Position', [100 100 550 550]);
    %figure('Position', [100 100 550 550]);
    h=circularGraph(temp_mat,'ColorMapNode',black_nodes,'ColorMapEdges',TempColor_dark,'Label',Brain_label2);

    positions=[];
        for p=1:length(h.Node)
        positions(p,:)=h.Node(p).Position;
        end

    hold on;
    scatter(positions(:,1),positions(:,2),30,funct_clust_colors,'filled');
    hold off;
    
    %set(gcf, 'Renderer','painters');
    %saveas(gcf,strcat('circularGraph_cluster_',group,'_',loom{L},'_av','.svg'));
    %saveas(gcf,strcat('circularGraph_cluster_',group,'_',loom{L},'_av','.emf'));
   end
end








