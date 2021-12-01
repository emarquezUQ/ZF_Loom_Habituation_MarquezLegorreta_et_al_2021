
%%%% this script is to analyze the 20ISI stimulus train for the fmr1 experiment for Marquez-Legorreta et al 2021   

%%% IMPORTANT: we need to have the data in excel sorted by fish or location,
%%% not by bin. 

[~,~,SUMMARY] = xlsread('FMRSUMMARYDOC_ORGANISEDFINALaftertroubleshooting.xlsx');


%%%%now for s20 fish

larct_s20=zeros(length(SUMMARY)-1,1);


for i=1:length(larct_s20)

    larct_s20(i,1)=SUMMARY{i+1,11}; %%%the + 1 is to get rid of the first row wich has the names
end

Tot_sum_s20=sum(larct_s20);
Tot_mean_s20=mean(larct_s20);


%%%this is to get rid of the extra fast movements (>30mm/s)  that happend on each timepoints (of 1s)
%%% this means that the data I am analyzing from now on is taking into
%%% acount the times a fast movement happend in a second, not the total
%%% number of fast movements. 
larct_s20_new=zeros(length(SUMMARY)-1,1);
for i=1:length(larct_s20)
    if larct_s20(i,1)>1
    larct_s20_new(i,1)=1;
    
    else
        larct_s20_new(i,1)=larct_s20(i,1);
    
    end
end

Tot_sum_s20_new=sum(larct_s20_new);
Tot_mean_s20_new=mean(larct_s20_new);



%%%his is to stack them by fish. I need to take into account the number of
%%%bins. 
larct_s20 = vec2mat(larct_s20_new, max([SUMMARY{:,2}]));
larct_s20 = larct_s20';


%%% to try to analyse it by loom

loom_times=[0,22,40,60,78,100,120,140,158,180]; %%% this is the times to add at the beggining of each periods to have the loom timepoints

%%%first loom moment
loom1st_moment_s20_perloom=[];
for i=1:size(larct_s20,2)
    start=300;
    for k=1:10
    loom1st_moment_s20_temp=sum(larct_s20(start+loom_times(k):start+loom_times(k)+5,i)); %%%% this is to get to each loom time point and sum the responses in that happend on the 5s of the loom stimuli 
    
    %%%and this loop is to only count as 1 response to the looms all the
    %%% fast movements
    if loom1st_moment_s20_temp>1
    loom1st_moment_s20_perloom(k,i)=1;
    else
    loom1st_moment_s20_perloom(k,i)=loom1st_moment_s20_temp;
    end
    
    end
end

loom1st_moment_perloom_inverted=loom1st_moment_s20_perloom';

%%%second loom moment
loom2nd_moment_s20_perloom=[];
for i=1:size(larct_s20,2)
    start=794;
    for k=1:10
    loom2nd_moment_s20_temp=sum(larct_s20(start+loom_times(k):start+loom_times(k)+5,i));
    if loom2nd_moment_s20_temp>1
    loom2nd_moment_s20_perloom(k,i)=1;
    else
    loom2nd_moment_s20_perloom(k,i)=loom2nd_moment_s20_temp;
    end
    
    
    end
end

loom2nd_moment_perloom_inverted=loom2nd_moment_s20_perloom';


sumS20=zeros(10,3);
for i=1:10
    sumS20(i,1)=sum(loom1st_moment_s20_perloom(i,:));
    sumS20(i,2)=sum(loom2nd_moment_s20_perloom(i,:));
    %sumS20(i,3)=sum(loom3rd_moment_s20_perloom(i,:));
end
ratioS20=sumS20/size(larct_s20,2);



All_perloom_s20=vertcat(loom1st_moment_s20_perloom,loom2nd_moment_s20_perloom);

%%



%%% a bit of general analysis

%%%per loom
figure;
hold on;
plot(mean(All_perloom_s20'));

hold off;


%%%per block

blocks_s20=horzcat(mean(mean(loom1st_moment_s20_perloom)),mean(mean(loom2nd_moment_s20_perloom)));

figure;
hold on;
plot(blocks_s20);

hold off;


%%

%%% to try to get the data of particular fish
%%% example with fish 3,6 and 11. 

%%% it works!

fish_names=unique(SUMMARY(2:end,1)); %%% to get the name of all the fish

fish_list=['c1-003';'c1-006';'c1-011']; %%% to select the fish we want

idx_fish_list=[];
for i=1:size(fish_list,1)
idx_fish_list(i)=find(contains(fish_names,fish_list(i,:)));
end

fish_list2=['c1-002';'c1-008';'c1-012']; %%% to select the fish we want

idx_fish_list2=[];
for i=1:size(fish_list2,1)
idx_fish_list2(i)=find(contains(fish_names,fish_list2(i,:)));
end


%%%per loom
figure;
hold on;
plot(mean(All_perloom_s20(:,idx_fish_list)'));
plot(mean(All_perloom_s20(:,idx_fish_list2)'));
hold off;

