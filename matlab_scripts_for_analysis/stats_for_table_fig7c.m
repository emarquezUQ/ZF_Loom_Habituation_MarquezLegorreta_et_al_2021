
%%%%% this script is to continue figures_n_stats_FMR1_commun_measures.m
%%%%% making extended data table with the stats of figure 7c
%%% it will have: 
% median diff; Chi2(2df); p(friedman); meanrank dif; 95%CIL; 95%CIH; p(bon)


%% just checking friedmans results. 

test_Mat=[];
for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
        
    %STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).friedman.p       
    
    if STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).friedman.p<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        test_Mat(brain,1)=1;
    else
        test_Mat=0;  
    end
    
    %STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).friedman.p;       
    
    if STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).friedman.p<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        test_Mat(brain,2)=1;
    else
        test_Mat(brain,2)=0;  
    end
    
    
    %%% prom
    
    %STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).friedman.p;     
    
    if STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).friedman.p<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        test_Mat(brain,3)=1;
    else
        test_Mat(brain,3)=0;  
    end
    
end


test_Mat2=[];
for clust=unique(ClustID_CL4(keepFmr1))'
    
    
    %STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).friedman.p;       
   
    if STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).friedman.p<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
       test_Mat2(clust,1)=1;
    else
        test_Mat2(clust,1)=0;  
    end
    
    %%% cohe
    
    %STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).friedman.p;       
    
    if STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).friedman.p<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        test_Mat2(clust,2)=1;
    else
        test_Mat2(clust,2)=0;  
    end
    
    
    %%% prom
    
    %STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).friedman.p;     
    
    if STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).friedman.p<0.05/13 %%% bonferroni correction 14 tests (friedman+4clust+9brainregions)
        test_Mat2(clust,3)=1;
    else
        test_Mat2(clust,3)=0;  
    end
        
end

%% adding flex data to table

Table=[];

for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
     
    Table(brain,1)=Brain_stats_test.difs(brain,1); %%% median difference i calculated
    Table(brain,2)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).friedman.tbl{2,5}; %% friedman's test chisquare. all have 2df
    Table(brain,3)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).friedman.p; %%% friedman's test p value
          
    Table(brain,4)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).multcomp.c(2,4);%%mean rank diff
    Table(brain,5)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).multcomp.c(2,3);%% 95% CI lower
    Table(brain,6)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).multcomp.c(2,5);%% 95% CI upper      
    Table(brain,7)=STATS_FMR1_all.flex.Flex_perBrain.(RegionList{brain}).multcomp.c(2,6);%% p value corrected for multiple comparisons (3 groups)
    
end

for clust=unique(ClustID_CL4(keepFmr1))'
    
    Table(9+clust,1)=Clust_stats_test.difs(clust,1); %%% median difference i calculated
    Table(9+clust,2)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).friedman.tbl{2,5}; %% friedman's test chisquare. all have 2df
    Table(9+clust,3)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).friedman.p; %%% friedman's test p value
          
    Table(9+clust,4)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).multcomp.c(2,4);%%mean rank diff
    Table(9+clust,5)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).multcomp.c(2,3);%% 95% CI Lower
    Table(9+clust,6)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).multcomp.c(2,5);%% 95% CI upper      
    Table(9+clust,7)=STATS_FMR1_all.flex.Flex_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6);%% p value corrected for multiple comparisons (3 groups)
        
end

%% adding coehsion data to table

Table2=[];

for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
     
    Table2(brain,1)=Brain_stats_test.difs(brain,2); %%% median difference i calculated
    Table2(brain,2)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).friedman.tbl{2,5}; %% friedman's test chisquare. all have 2df
    Table2(brain,3)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).friedman.p; %%% friedman's test p value
          
    Table2(brain,4)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).multcomp.c(2,4);%%mean rank diff
    Table2(brain,5)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).multcomp.c(2,3);%% 95% CI lower
    Table2(brain,6)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).multcomp.c(2,5);%% 95% CI upper      
    Table2(brain,7)=STATS_FMR1_all.cohe.Cohe_perBrain.(RegionList{brain}).multcomp.c(2,6);%% p value corrected for multiple comparisons (3 groups)
    
end

for clust=unique(ClustID_CL4(keepFmr1))'
    
    Table2(9+clust,1)=Clust_stats_test.difs(clust,2); %%% median difference i calculated
    Table2(9+clust,2)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).friedman.tbl{2,5}; %% friedman's test chisquare. all have 2df
    Table2(9+clust,3)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).friedman.p; %%% friedman's test p value
          
    Table2(9+clust,4)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).multcomp.c(2,4);%%mean rank diff
    Table2(9+clust,5)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).multcomp.c(2,3);%% 95% CI Lower
    Table2(9+clust,6)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).multcomp.c(2,5);%% 95% CI upper      
    Table2(9+clust,7)=STATS_FMR1_all.cohe.Cohe_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6);%% p value corrected for multiple comparisons (3 groups)
        
end

%% adding promiscuity data to table

Table3=[];

for brain=unique(NodesFmr1.Nod_brainID(keepFmr1))'
     
    Table3(brain,1)=Brain_stats_test.difs(brain,3); %%% median difference i calculated
    Table3(brain,2)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).friedman.tbl{2,5}; %% friedman's test chisquare. all have 2df
    Table3(brain,3)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).friedman.p; %%% friedman's test p value
          
    Table3(brain,4)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).multcomp.c(2,4);%%mean rank diff
    Table3(brain,5)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).multcomp.c(2,3);%% 95% CI lower
    Table3(brain,6)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).multcomp.c(2,5);%% 95% CI upper      
    Table3(brain,7)=STATS_FMR1_all.prom.Prom_perBrain.(RegionList{brain}).multcomp.c(2,6);%% p value corrected for multiple comparisons (3 groups)
    
end

for clust=unique(ClustID_CL4(keepFmr1))'
    
    Table3(9+clust,1)=Clust_stats_test.difs(clust,3); %%% median difference i calculated
    Table3(9+clust,2)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).friedman.tbl{2,5}; %% friedman's test chisquare. all have 2df
    Table3(9+clust,3)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).friedman.p; %%% friedman's test p value
          
    Table3(9+clust,4)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).multcomp.c(2,4);%%mean rank diff
    Table3(9+clust,5)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).multcomp.c(2,3);%% 95% CI Lower
    Table3(9+clust,6)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).multcomp.c(2,5);%% 95% CI upper      
    Table3(9+clust,7)=STATS_FMR1_all.prom.Prom_perClust4.(clustnames_CL4{clust}).multcomp.c(2,6);%% p value corrected for multiple comparisons (3 groups)
        
end

