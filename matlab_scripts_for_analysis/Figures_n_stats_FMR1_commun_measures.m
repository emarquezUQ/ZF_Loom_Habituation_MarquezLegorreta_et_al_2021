


%%%%% this script is to make figures and gather the data to pass to prism
%%%%% for the the community detection analysis. the data is based on the
%%%%% GnO_limits_N_commStats_FMR1.m script. 

%%%%% is for Figures 7a-f and Supplementary figure 10.

%%%%% NOTE: apart from loading some datasets it's also needed to run the
%%%%% parts of the "GnO_limits_N_commStats_FMR1.m" script to get the 
%%%%% optimization of gamma (g) and omega(o). 


load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;

S_cons_FMR1=load('S_cons_FMR1_G14O07.mat');
S_cons_FnS_All=load('S_cons_FnS_All_G14O07.mat');
load('S_cons_testing_GnO_measures6.mat');

cbrewer()

[RdBu]=cbrewer('div','RdBu',101);

[set1]=cbrewer('qual','Set1',40);

set1=flip(set1);

[spectral]=cbrewer('div','Spectral',40);
spectral=flip(spectral);


%%
%%%% with fasthab subtypes
oldclusttag=[4 2 5 6 1 7];
ClustID=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    ClustID(idx_temp)=i;
   
end
 
clustnames={'fasthab1','fasthab2','fasthab3','modhab','weakhab','inhib'};

%%%% with fasthab merged
oldclusttag=[4 2 5 6 1 7];
ClustID_CL4=zeros(size(Nodes.Nod_clustID(keep)));
for i=1:length(unique(Nodes.Nod_clustID(keep)))
    
    idx_temp=find(Nodes.Nod_clustID(keep)==oldclusttag(i));
    
    if i<4    
    ClustID_CL4(idx_temp)=1;
    elseif i==4
   ClustID_CL4(idx_temp)=2;
   elseif i==5
   ClustID_CL4(idx_temp)=3;
   elseif i==6
   ClustID_CL4(idx_temp)=4;
    end
end
 
clustnames_CL4={'fasthab','modhab','weakhab','inhib'};


%% dynamic network community detection


figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
sgtitle(strcat('g=',num2str(1.6),'/','o=',num2str(0.9)));
counter=1;
for group=[3 1 2]
S_good=Big_FMR1_OPT{16,9}.(groupnames{group}).S_cons;

%low=min(min(S_good));
%high=max(max(S_good));

 subplot(1,3,counter);imagesc(S_good);colormap(spectral);caxis([1 40]);%colorbar;
  title(groupnames{group}); 
 counter=counter+1;
    
 end   

saveas(gcf,'communities_fmr1_g16o09_2.svg')
close 


spectral=spectral(randperm(40,40),:);
figure;colormap(spectral);caxis([1 40]);colorbar;
saveas(gcf,'communities_fmr1_g16o09colorbar.svg')
close 



%% making a figure with the flex,cohe and prom matrices of all datasets


counter=1;
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])

for group=[3 1 2]
   subplot(3,3,counter);imagesc(FlexMat.(groupnames{group}));caxis([0 1]);colormap(inferno);%colorbar;
   tempMed=median(FlexMat.(groupnames{group})(:));
   tempSD=std(FlexMat.(groupnames{group})(:));
   title(strcat('flex/',groupnames{group},'/Med=',num2str(round(tempMed,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end

for group=[3 1 2]
   subplot(3,3,counter);imagesc(CoheMat.(groupnames{group}));caxis([0 1]);colormap(inferno);%colorbar;
   tempMed=median(CoheMat.(groupnames{group})(:));
   tempSD=std(CoheMat.(groupnames{group})(:));
   title(strcat('Cohe/',groupnames{group},'/Med=',num2str(round(tempMed,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end

for group=[3 1 2]
   subplot(3,3,counter);imagesc(PromMat.(groupnames{group}));caxis([0 1]);colormap(inferno);%colorbar;
   tempMed=median(PromMat.(groupnames{group})(:));
   tempSD=std(PromMat.(groupnames{group})(:));
   title(strcat('prom/',groupnames{group},'/Med=',num2str(round(tempMed,2)),'/SD=',num2str(round(tempSD,2))));
   
   counter=counter+1;
end

saveas(gcf,'Measures_FMR1.svg');
%saveas(gcf,'Measures_FMR1_wcolorbar.svg');
close



%% making a structure for the stat results

STATS_FMR1_all=struct;

%% getting the values for flexibility plot 

%%%%% the selectiono of gammas and omegas is based on:
%%%%% GnO_limits_N_commStats_FMR1.m after filtering the good g and o

%%%% or load these:
%o=[5 5 5 7 7 7 7 7 7 8 9 9 10 10 10 11 11 12 12 13 13 13 13 14 14 14 15 15 15 15 16 16 16 17];
%g=[14 15 16 13 15 16 17 18 19 16 16 18 17 18 20 11 16 11 18 14 19 20 21 14 17 20 16 18 19 20 16 17 19 17];

counter=1;
test=struct;
for i=1:length(o)
%%%% comparing flexibility


    
    for group=1:3   
    
    test.(groupnames{group})(counter,:)=Big_FMR1_OPT{g(i),o(i)}.(groupnames{group}).flex;
    end
    counter=counter+1;
 %end
end


%%
groups_flexibilityTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    groups_flexibilityTest(1:length(o),counter)=mean(test.(groupnames{group})'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes. 
    %groups_flexibilityTest(1:90,counter)=mean(test.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end

figure;boxplot(groups_flexibilityTest);

STATS_FMR1_all.flex.general.groups=groups_flexibilityTest;


    
%     [p, tbl,
%     stats]=ranksum(groups_flexibilityTest(:,1),groups_flexibilityTest(:,3));
%      %%%%% the ranksum test also gives sig. 
%     p
    
    [p, tbl, stats]=friedman(groups_flexibilityTest);
    
    STATS_FMR1_all.flex.general.friedman.p=p;
    STATS_FMR1_all.flex.general.friedman.p=tbl;
    STATS_FMR1_all.flex.general.friedman.p=stats;
    
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
    STATS_FMR1_all.flex.general.multcomp.c=c;
    STATS_FMR1_all.flex.general.multcomp.m=m;
    
%%% then I took the groups_flexibility to graphpad to graph it


%% brain region and clusters comparison. For flexibility

%%% I will use graphpad for the plots

%%% 
Flex_perBrain=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
       
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Flex_perBrain.(RegionList{brain})=temp_groups;
    
    STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).groups=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    
    STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).friedman.p=p;
    STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).friedman.tbl=tbl;
    STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).friedman.stats=stats;
    
    figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni');
        
    STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).multcomp.c=c;
    STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).multcomp.m=m;
end

%%% 

%%%%%%%% now by cluster with the CL4
Flex_perClust4=struct;
for clust=unique(ClustID_CL4(keepFmr1))'
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_groups=[];
    
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Flex_perClust4.(clustnames_CL4{clust})=temp_groups;
    
    
    STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).groups=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    
    STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).friedman.p=p;
    STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).friedman.tbl=tbl;
    STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).friedman.stats=stats;
    
    [c, m]=multcompare(stats,'CType','bonferroni');
    
    STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).multcomp.c=c;
    STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).multcomp.m=m;
    
end


%% relative difference in flex plot

[RdBu]=cbrewer('div','RdBu',101);

    groups_flexibilityTest=NaN(90,3);
%groups_flexibilityTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    %groups_flexibilityTest(1:length(o),counter)=mean(test.(groupnames{group})'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes. 
    groups_flexibilityTest(1:90,counter)=mean(test.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end
    
flex_dif=groups_flexibilityTest(:,1)-groups_flexibilityTest(:,3);
%figure;histogram(flex_dif);

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
%   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,flex_dif,'filled');colormap(RdBu);caxis([-0.4 0.4]);%colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;

saveas(gcf,'Flex_WTvsfmr1_GnO_limits.svg')
%saveas(gcf,'Flex_WTvsfmr1_GnO_limits_wcolorbar.svg')
close 

%% now for cohesion

counter=1;
test2=struct;
for i=1:length(o)
% for g=7:10
% for o=8:10
    
    for group=1:3   
    
    test2.(groupnames{group})(counter,:)=Big_FMR1_OPT{g(i),o(i)}.(groupnames{group}).node_cohesion;
    end
    counter=counter+1;
% end
% end
end

%groups_CoheTest=NaN(90,3);
groups_CoheTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    groups_CoheTest(1:length(o),counter)=mean(test2.(groupnames{group})');
counter=counter+1;
end

figure;boxplot(groups_CoheTest);


STATS_FMR1_all.cohe.general.groups=groups_CoheTest;

    [p, tbl, stats]=friedman(groups_CoheTest);
    figure;
    
    STATS_FMR1_all.cohe.general.friedman.p=p;
    STATS_FMR1_all.cohe.general.friedman.p=tbl;
    STATS_FMR1_all.cohe.general.friedman.p=stats;
    
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
    STATS_FMR1_all.cohe.general.multcomp.c=c;
    STATS_FMR1_all.cohe.general.multcomp.m=m;

       
    
%%% 
Cohe_perBrain=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
       
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test2.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Cohe_perBrain.(RegionList{brain})=temp_groups;
       
    STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).groups=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    
    STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).friedman.p=p;
    STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).friedman.tbl=tbl;
    STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).friedman.stats=stats;
    
    figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni');
        
    STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).multcomp.c=c;
    STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).multcomp.m=m;
end

%%% 


%%%%%%%% now by cluster with the CL4
Cohe_perClust4=struct;
for clust=unique(ClustID_CL4(keepFmr1))'
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_groups=[];
    
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test2.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Cohe_perClust4.(clustnames_CL4{clust})=temp_groups;
    
    
    STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).groups=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    
    STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).friedman.p=p;
    STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).friedman.tbl=tbl;
    STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).friedman.stats=stats;
    
    [c, m]=multcompare(stats,'CType','bonferroni');
    
    STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).multcomp.c=c;
    STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).multcomp.m=m;
    
end


    %%%% making figures
    %%%% for this one I need to add the cbrewer
    
    groups_CoheTest=NaN(90,3);
%groups_flexibilityTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    %groups_flexibilityTest(1:length(o),counter)=mean(test.(groupnames{group})'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes. 
    groups_CoheTest(1:90,counter)=mean(test2.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end
    
cohe_dif=groups_CoheTest(:,1)-groups_CoheTest(:,3);
%figure;histogram(flex_dif);

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
%   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,cohe_dif,'filled');colormap(RdBu);caxis([-0.4 0.4]);%colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;


saveas(gcf,'Cohe_WTvsfmr1_GnO_limits.svg')
%saveas(gcf,'Cohe_WTvsfmr1_GnO_limits_wcolorbar.svg')

close


%% now for promiscuity

counter=1;
test3=struct;
for i=1:length(o)

    
    for group=1:3   
    
    test3.(groupnames{group})(counter,:)=Big_FMR1_OPT{g(i),o(i)}.(groupnames{group}).P;
    end
    counter=counter+1;
% end
end


groups_PromTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    groups_PromTest(1:length(o),counter)=mean(test3.(groupnames{group})');
counter=counter+1;
end

figure;boxplot(groups_PromTest);


STATS_FMR1_all.prom.general.groups=groups_PromTest;

    [p, tbl, stats]=friedman(groups_PromTest);
    
    
    STATS_FMR1_all.prom.general.friedman.p=p;
    STATS_FMR1_all.prom.general.friedman.p=tbl;
    STATS_FMR1_all.prom.general.friedman.p=stats;
    
    figure;
    [c, m]=multcompare(stats,'CType','bonferroni');

    
    STATS_FMR1_all.prom.general.multcomp.c=c;
    STATS_FMR1_all.prom.general.multcomp.m=m;   
    
       
%%% 
Prom_perBrain=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    temp_idx=find(NodesFmr1.Nod_brainID(keepFmr1)==brain);
    temp_groups=[];
       
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test3.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Prom_perBrain.(RegionList{brain})=temp_groups;
     
    STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).groups=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    
    STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).friedman.p=p;
    STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).friedman.tbl=tbl;
    STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).friedman.stats=stats;
    
    figure;title(RegionList{brain});
    [c, m]=multcompare(stats,'CType','bonferroni');
        
    STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).multcomp.c=c;
    STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).multcomp.m=m;
end

%%% 


%%%%%%%% now by cluster with the CL4
Prom_perClust4=struct;
for clust=unique(ClustID_CL4(keepFmr1))'
    temp_idx=find(ClustID_CL4(keepFmr1)==clust);
    temp_groups=[];
    
        counter=1;
        for group=[3 1 2]
            temp_groups(1:length(o),counter)=mean(test3.(groupnames{group})(:,temp_idx)'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes.
            
            counter=counter+1;
            
        end
    
    Prom_perClust4.(clustnames_CL4{clust})=temp_groups;
    
    
    
    STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).groups=temp_groups;
    
    figure;boxplot(temp_groups);
    [p, tbl, stats]=friedman(temp_groups);
    
    STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).friedman.p=p;
    STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).friedman.tbl=tbl;
    STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).friedman.stats=stats;
    
    [c, m]=multcompare(stats,'CType','bonferroni');
    
    STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).multcomp.c=c;
    STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).multcomp.m=m;
    
end

    %%%% making figures
    %%%% for this one I need to add the cbrewer
    
    groups_PromTest=NaN(90,3);
%groups_flexibilityTest=NaN(length(o),3);
counter=1;
for group=[3 1 2]
    %groups_flexibilityTest(1:length(o),counter)=mean(test.(groupnames{group})'); %%%% Gilles thinks the units have to be the gammaxomega combinations and not the nodes. 
    groups_PromTest(1:90,counter)=mean(test3.(groupnames{group}));
    %groups_flexibilityTest(1:90,counter)=test.(groupnames{group});
counter=counter+1;
end
    
Prom_dif=groups_PromTest(:,1)-groups_PromTest(:,3);
%figure;histogram(flex_dif);

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
plot(Zbrain_brainMask2D(:,1),Zbrain_brainMask2D(:,2),'k');xlim([400 1350])
%  hold on;
%  gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),NodesFmr1.Nod_brainID(keepFmr1,:),[],[],21); 
%   hold on;
%   gscatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),ClustID(keepFmr1),'gggbrm',[],21); 
hold on;
scatter(NodesFmr1.Nod_coor(keepFmr1,1),NodesFmr1.Nod_coor(keepFmr1,2),30,Prom_dif,'filled');colormap(RdBu);caxis([-0.2 0.2]);colorbar;
view(-90,90);
%title('top 1SD dif flex');
hold off;

 
%saveas(gcf,'Prom_WTvsfmr1_GnO_limits.svg')
saveas(gcf,'Prom_WTvsfmr1_GnO_limits_wcolorbar.svg')

close   

%%


%% making a brain region and clusters stat comparison heatmap 

%% for brain regions

Brain_stats_test=struct;
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
    
    
    %%% flex 
    temp_groups=[];
    temp_groups(:,1)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).groups(:,1);  %%WT
    temp_groups(:,2)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).groups(:,3);  %%fmr1
    Brain_stats_test.difs(brain,1)=median(temp_groups(:,1))-median(temp_groups(:,2));
    
    %STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).friedman.p;       
    Brain_stats_test.mult_p(brain,1)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).multcomp.c(2,6);
    
    if STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).multcomp.c(2,6)<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        Brain_stats_test.sig_corrected(brain,1)=1;
    else
        Brain_stats_test.sig_corrected(brain,1)=0;  
    end
    
    
    %%% cohe
    temp_groups=[];
    temp_groups(:,1)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).groups(:,1);  %%WT
    temp_groups(:,2)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).groups(:,3);  %%fmr1
    Brain_stats_test.difs(brain,2)=median(temp_groups(:,1))-median(temp_groups(:,2));
    
    %STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).friedman.p;       
    Brain_stats_test.mult_p(brain,2)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).multcomp.c(2,6);
    
    if STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).multcomp.c(2,6)<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        Brain_stats_test.sig_corrected(brain,2)=1;
    else
        Brain_stats_test.sig_corrected(brain,2)=0;  
    end
    
    
    %%% prom
    temp_groups=[];
    temp_groups(:,1)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).groups(:,1);  %%WT
    temp_groups(:,2)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).groups(:,3);  %%fmr1
    Brain_stats_test.difs(brain,3)=median(temp_groups(:,1))-median(temp_groups(:,2));
    
    %STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).friedman.p;     
    Brain_stats_test.mult_p(brain,3)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).multcomp.c(2,6);
    
    if STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).multcomp.c(2,6)<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        Brain_stats_test.sig_corrected(brain,3)=1;
    else
        Brain_stats_test.sig_corrected(brain,3)=0;  
    end
    
end

max(Brain_stats_test.difs(:))
min(Brain_stats_test.difs(:))
%%% a heat map from -0.25 to 0.25

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
imagesc(Brain_stats_test.difs);colormap(RdBu);caxis([-0.25 0.25]);colorbar;

%saveas(gcf,'FMR1_Brain_comm_measures_difs_heatmap.svg')
saveas(gcf,'FMR1_Brain_comm_measures_difs_heatmap_wcolorbar.svg')

close   

%% for clusters 

%%%%%%%% now by cluster with the CL4


Clust_stats_test=struct;


for clust=unique(ClustID_CL4(keepFmr1))'
    
    
    %%% flex 
    temp_groups=[];
    temp_groups(:,1)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).groups(:,1);  %%WT
    temp_groups(:,2)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).groups(:,3);  %%fmr1
    Clust_stats_test.difs(clust,1)=median(temp_groups(:,1))-median(temp_groups(:,2));
    
    %STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).friedman.p;       
    Clust_stats_test.mult_p(clust,1)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6);
    
    if STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6)<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        Clust_stats_test.sig_corrected(clust,1)=1;
    else
        Clust_stats_test.sig_corrected(clust,1)=0;  
    end
    
    
    %%% cohe
    temp_groups=[];
    temp_groups(:,1)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).groups(:,1);  %%WT
    temp_groups(:,2)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).groups(:,3);  %%fmr1
    Clust_stats_test.difs(clust,2)=median(temp_groups(:,1))-median(temp_groups(:,2));
    
    %STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).friedman.p;       
    Clust_stats_test.mult_p(clust,2)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6);
    
    if STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6)<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        Clust_stats_test.sig_corrected(clust,2)=1;
    else
        Clust_stats_test.sig_corrected(clust,2)=0;  
    end
    
    
    %%% prom
    temp_groups=[];
    temp_groups(:,1)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).groups(:,1);  %%WT
    temp_groups(:,2)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).groups(:,3);  %%fmr1
    Clust_stats_test.difs(clust,3)=median(temp_groups(:,1))-median(temp_groups(:,2));
    
    %STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).friedman.p;     
    Clust_stats_test.mult_p(clust,3)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6);
    
    if STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6)<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        Clust_stats_test.sig_corrected(clust,3)=1;
    else
        Clust_stats_test.sig_corrected(clust,3)=0;  
    end
        
end


max(Clust_stats_test.difs(:))
min(Clust_stats_test.difs(:))
%%% a heat map from -0.25 to 0.25

figure;set(gcf, 'Position',  [200, 200, 350, 450]);
imagesc(Clust_stats_test.difs);colormap(RdBu);caxis([-0.25 0.25]);colorbar;

%saveas(gcf,'FMR1_Clust_comm_measures_difs_heatmap.svg')
saveas(gcf,'FMR1_Clust_comm_measures_difs_heatmap_wcolorbar.svg')

close   

%%

save('STATS_FMR1_all.mat','STATS_FMR1_all');

%% exploring the WTvsFMR1 diff per nodes

figure;
scatter3(flex_dif,cohe_dif,Prom_dif,10,'filled');

commMeasures=[];

commMeasures(:,1)=flex_dif;
commMeasures(:,2)=cohe_dif;
commMeasures(:,3)=Prom_dif;

%%%% raster plot of the differences sorted by cluster
max(commMeasures(:))
min(commMeasures(:))

%%%%% making a raster plot figure the same size as the community detection
%%%%% examples 
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1);imagesc(commMeasures);colormap(RdBu);caxis([-0.4 0.4]);colorbar;
title('Relative differences'); 
%saveas(gcf,'WTvsfmr1_GnO_limits_RelDif.svg')
saveas(gcf,'WTvsfmr1_GnO_limits_RelDif_wcolorbar.svg')

%%%% and also a bar with the cluster code

clustcolorbar=zeros(90,3);
temp_idx=[];
for clust=1:4
   temp_idx=find(ClustID_CL4(keepFmr1)==clust);
   for i=1:length(temp_idx)
   if clust==1 
   clustcolorbar(temp_idx(i),:)=[0 1 0]; %%green
   elseif clust==2
    clustcolorbar(temp_idx(i),:)=[0 0 1]; %%blue  
   elseif clust==3
    clustcolorbar(temp_idx(i),:)=[1 0 0]; %%red   
   else
   clustcolorbar(temp_idx(i),:)=[1 0 1]; %% magenta
   end
end
end    
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(1,3,2);colormap(flip(clustcolorbar));colorbar;
saveas(gcf,'FMR1_clustcolorbar.svg')    
 

%% now for the supp figure about Q optimization and selection of Gamma and Omega


%figure;set(gcf,'units','normalized','outerposition',[0 0 1 1])
figure;set(gcf, 'Position',  [50, 50, 800, 750]);
counter=1;
for group=[3 1 2]
   subplot(3,3,counter);imagesc(FlexMat.(groupnames{group}));caxis([0 1]);
   title(strcat('flex/',groupnames{group},'/med=',num2str(median(FlexMat.(groupnames{group})(:)))));colormap(inferno);%colorbar;
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);imagesc(CoheMat.(groupnames{group}));caxis([0 1]);
   title(strcat('cohe/',groupnames{group},'/med=',num2str(median(CoheMat.(groupnames{group})(:)))));colormap(inferno);%colorbar;
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);imagesc(PromMat.(groupnames{group}));caxis([0 1]);
   title(strcat('prom/',groupnames{group},'/med=',num2str(median(PromMat.(groupnames{group})(:)))));colormap(inferno);%colorbar;
   
   counter=counter+1;
end

saveas(gcf,'Measures_FMR1_extended.svg')

figure;
caxis([0 1]);colormap(inferno);colorbar;
saveas(gcf,'infernocolorbar.svg')


%% for Q optimization


%%%% getting an example of Q opt

group=3; %%% WTs

temp_dif=QoptFMR1_All.(groupnames{group}).meanQopt2-QoptFMR1_All.(groupnames{group}).meanQopt2_null2;
    
    figure;imagesc(temp_dif);
           
temp_relative_var=-(QoptFMR1_All.(groupnames{group}).varQopt2-max(QoptFMR1_All.(groupnames{group}).varQopt2(:)));
figure;imagesc(temp_relative_var);

thegood=temp_dif.*temp_relative_var;

%figure;imagesc(thegood);

thebest=thegood;thebest(thegood<0)=0;

figure;imagesc(thebest);


figure;set(gcf, 'Position',  [50, 50, 1050, 750]);
subplot(3,3,1);imagesc(temp_dif);colormap(inferno);colorbar; %%% WT example of Q-Qt
subplot(3,3,2);imagesc(temp_relative_var);colormap(inferno);colorbar; %%% WT example of relative Variance
subplot(3,3,3);imagesc(thebest);colormap(inferno);colorbar; %%% mixing variance and opt Q and setting negatives to 0
%title(strcat('mean=',num2str(mean(Opt_new2(:)))));
%saveas(gcf,'Qopt_WTnGnO_fmr1_1.svg');
%figure;set(gcf, 'Position',  [50, 50, 1450, 400]);
subplot(3,3,4);imagesc(Opt);colormap(inferno);colorbar; %%% average of the 3 datasets
subplot(3,3,5);imagesc(Opt_new2);colormap(inferno);colorbar; %%% after limiting GnO
subplot(3,3,6);imagesc(Opt_best2);colormap(inferno);colorbar; %%% selected sample
saveas(gcf,'Qopt_WTnGnO_fmr1.svg');



