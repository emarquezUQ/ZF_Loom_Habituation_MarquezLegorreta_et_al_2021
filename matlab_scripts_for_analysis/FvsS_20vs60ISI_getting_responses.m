%%%% this script is to analyze the fast (f) vs slow (s) and the 20ISI vs 60ISI data from Marquez-Legorreta et al, 2021 

flooms_var20ISI={};
flooms_var60ISI={};
slooms_var20ISI={};
slooms_var60ISI={};

%%%%this is to read the specific excel spread sheet and get both the text
%%%%and numeric data:
%[~,~,slooms_var60ISI] = xlsread('slooms_var60ISI.xlsx','together');

%%%% first I will save the data that I copied from excel
 %save ('raw_FnSlooms_var20n60ISIdata.mat','flooms_var20ISI','flooms_var60ISI','slooms_var20ISI','slooms_var60ISI');

 
load raw_FnSlooms_var20n60ISIdata.mat



%%%%% first with the f20 fish

larct_f20=zeros(length(flooms_var20ISI)-1,1);


for i=1:length(larct_f20)

    larct_f20(i,1)=flooms_var20ISI{i+1,11}; %%%the + 1 is to get rid of the first row wich has the names
end

Tot_sum_f20=sum(larct_f20);
Tot_mean_f20=mean(larct_f20);


%%%this is to get rid of the extra fast movements (>30mm/s)  that happend on each timepoints (of 1s)
%%% this means that the data I am analyzing from now on is taking into
%%% acount the times a fast movement happend in a second, not the total
%%% number of fast movements. 
larct_f20_new=zeros(length(flooms_var20ISI)-1,1);
for i=1:length(larct_f20)
    if larct_f20(i,1)>1
    larct_f20_new(i,1)=1;
    
    else
        larct_f20_new(i,1)=larct_f20(i,1);
    
    end
end

Tot_sum_f20_new=sum(larct_f20_new);
Tot_mean_f20_new=mean(larct_f20_new);



%%%his is to stack them by fish. I need to take into account the number of
%%%bins. 
larct_f20 = vec2mat(larct_f20_new, max([flooms_var20ISI{:,2}]));
larct_f20 = larct_f20';






%%%% now I will try to analyze how many fast responses where in the pre
%%%% and rest moments. And in the 3 loom periods. 

%%% the pre moment

pre_moment_f20=[];
for i=1:size(larct_f20,2)
pre_moment_f20_temp=sum(larct_f20(1:299,i));
pre_moment_f20(1,i)=pre_moment_f20_temp;
end

mean_f20_pre=mean(pre_moment_f20);

%%% the rest moments

rest1_moment_f20=[];
for i=1:size(larct_f20,2)
rest1_moment_f20_temp=sum(larct_f20(494:793,i));
rest1_moment_f20(1,i)=rest1_moment_f20_temp;
end
mean_f20_rest1=mean(rest1_moment_f20);

%%% the rest moment

rest2_moment_f20=[];
for i=1:size(larct_f20,2)
rest2_moment_f20_temp=sum(larct_f20(988:1287,i));
rest2_moment_f20(1,i)=rest2_moment_f20_temp;
end
mean_f20_rest2=mean(rest2_moment_f20);



%%%first loom moment
loom1st_moment_f20=[];
for i=1:size(larct_f20,2)
loom1st_moment_f20_temp=sum(larct_f20(300:493,i));
loom1st_moment_f20(1,i)=loom1st_moment_f20_temp;
end

mean_f20_loom1st=mean(loom1st_moment_f20);

%%%second loom moment
loom2nd_moment_f20=[];
for i=1:size(larct_f20,2)
loom2nd_moment_f20_temp=sum(larct_f20(794:987,i));
loom2nd_moment_f20(1,i)=loom2nd_moment_f20_temp;
end

mean_f20_loom2nd=mean(loom2nd_moment_f20);

%%%third loom moment
loom3rd_moment_f20=[];
for i=1:size(larct_f20,2)
loom3rd_moment_f20_temp=sum(larct_f20(1288:1482,i));
loom3rd_moment_f20(1,i)=loom3rd_moment_f20_temp;
end
mean_f20_loom3rd=mean(loom3rd_moment_f20);


base_and_loom_resp_f20=[];

base_and_loom_resp_f20(1,1)=sum(pre_moment_f20);
base_and_loom_resp_f20(2,1)=sum(loom1st_moment_f20);
base_and_loom_resp_f20(3,1)=sum(rest1_moment_f20);
base_and_loom_resp_f20(4,1)=sum(loom2nd_moment_f20);
base_and_loom_resp_f20(5,1)=sum(rest2_moment_f20);
base_and_loom_resp_f20(6,1)=sum(loom3rd_moment_f20);




loom_resp_ratio_f20=(base_and_loom_resp_f20/size(larct_f20,2));



figure;
plot(loom_resp_ratio_f20);
%hold on;
%plot(loom_resp_ratio_ctrl);


%%%% now to try to analyse it by loom

loom_times=[0,22,40,60,78,100,120,140,158,180]; %%% this is the times to add at the beggining of each periods to have the loom timepoints

%%%first with fmr1 mutants

%%%first loom moment
loom1st_moment_f20_perloom=[];
for i=1:size(larct_f20,2)
    start=300;
    for k=1:10
    loom1st_moment_f20_temp=sum(larct_f20(start+loom_times(k):start+loom_times(k)+3,i)); %%%% this is to get to each loom time point and sum the responses in that happend on the 5s of the loom stimuli 
    
    %%%and this loop is to only count as 1 response to the looms all the
    %%% fast movements
    if loom1st_moment_f20_temp>1
    loom1st_moment_f20_perloom(k,i)=1;
    else
    loom1st_moment_f20_perloom(k,i)=loom1st_moment_f20_temp;
    end
    
    end
end

%%%second loom moment
loom2nd_moment_f20_perloom=[];
for i=1:size(larct_f20,2)
    start=794;
    for k=1:10
    loom2nd_moment_f20_temp=sum(larct_f20(start+loom_times(k):start+loom_times(k)+3,i));
    if loom2nd_moment_f20_temp>1
    loom2nd_moment_f20_perloom(k,i)=1;
    else
    loom2nd_moment_f20_perloom(k,i)=loom2nd_moment_f20_temp;
    end
    
    
    end
end


%%%third loom moment
loom3rd_moment_f20_perloom=[];
for i=1:size(larct_f20,2)
    start=1288;
    for k=1:10
    loom3rd_moment_f20_temp=sum(larct_f20(start+loom_times(k):start+loom_times(k)+3,i));
    if loom3rd_moment_f20_temp>1
    loom3rd_moment_f20_perloom(k,i)=1;
    else
    loom3rd_moment_f20_perloom(k,i)=loom3rd_moment_f20_temp;
    end
    
    
    end
end




sumF20=zeros(10,3);
for i=1:10
    sumF20(i,1)=sum(loom1st_moment_f20_perloom(i,:));
    sumF20(i,2)=sum(loom2nd_moment_f20_perloom(i,:));
    sumF20(i,3)=sum(loom3rd_moment_f20_perloom(i,:));
end
ratioF20=sumF20/size(larct_f20,2);






%%%% and to look for the sound response!


sound_f20=[];
sound_times=[0,5,10];
for i=1:size(larct_f20,2)
    start=1262;
    for k=1:3 %%%cause there are 3 sounds
    
    
    sound_f20_temp=sum(larct_f20(start+sound_times(k):start+sound_times(k)+2,i));
    if sound_f20_temp>1
    sound_f20(k,i)=1;
    else
    sound_f20(k,i)=sound_f20_temp;
    end
    end
end
mean_sound_f20=mean(sound_f20');



All_perloom_f20=vertcat(loom1st_moment_f20_perloom,loom2nd_moment_f20_perloom,loom3rd_moment_f20_perloom);



%%

%%%%% now with the f60 fish



larct_f60=zeros(length(flooms_var60ISI)-1,1);


for i=1:length(larct_f60)

    larct_f60(i,1)=flooms_var60ISI{i+1,11}; %%%the + 1 is to get rid of the first row wich has the names
end

Tot_sum_f60=sum(larct_f60);
Tot_mean_f60=mean(larct_f60);


%%%this is to get rid of the extra fast movements (>30mm/s)  that happend on each timepoints (of 1s)
%%% this means that the data I am analyzing from now on is taking into
%%% acount the times a fast movement happend in a second, not the total
%%% number of fast movements. 
larct_f60_new=zeros(length(flooms_var60ISI)-1,1);
for i=1:length(larct_f60)
    if larct_f60(i,1)>1
    larct_f60_new(i,1)=1;
    
    else
        larct_f60_new(i,1)=larct_f60(i,1);
    
    end
end

Tot_sum_f60_new=sum(larct_f60_new);
Tot_mean_f60_new=mean(larct_f60_new);



%%%his is to stack them by fish. I need to take into account the number of
%%%bins. 
larct_f60 = vec2mat(larct_f60_new, max([flooms_var60ISI{:,2}]));
larct_f60 = larct_f60';






%%%% now I will try to analyze how many fast responses where in the pre
%%%% and rest moments. And in the 3 loom periods. 

%%% the pre moment

pre_moment_f60=[];
for i=1:size(larct_f60,2)
pre_moment_f60_temp=sum(larct_f60(1:299,i));
pre_moment_f60(1,i)=pre_moment_f60_temp;
end

mean_f60_pre=mean(pre_moment_f60);

%%% the rest moments

rest1_moment_f60=[];
for i=1:size(larct_f60,2)
rest1_moment_f60_temp=sum(larct_f60(854:1153,i));
rest1_moment_f60(1,i)=rest1_moment_f60_temp;
end
mean_f60_rest1=mean(rest1_moment_f60);

%%% the rest moment

rest2_moment_f60=[];
for i=1:size(larct_f60,2)
rest2_moment_f60_temp=sum(larct_f60(1708:2007,i));
rest2_moment_f60(1,i)=rest2_moment_f60_temp;
end
mean_f60_rest2=mean(rest2_moment_f60);



%%%first loom moment
loom1st_moment_f60=[];
for i=1:size(larct_f60,2)
loom1st_moment_f60_temp=sum(larct_f60(300:853,i));
loom1st_moment_f60(1,i)=loom1st_moment_f60_temp;
end

mean_f60_loom1st=mean(loom1st_moment_f60);

%%%second loom moment
loom2nd_moment_f60=[];
for i=1:size(larct_f60,2)
loom2nd_moment_f60_temp=sum(larct_f60(1154:1707,i));
loom2nd_moment_f60(1,i)=loom2nd_moment_f60_temp;
end

mean_f60_loom2nd=mean(loom2nd_moment_f60);

%%%third loom moment
loom3rd_moment_f60=[];
for i=1:size(larct_f60,2)
loom3rd_moment_f60_temp=sum(larct_f60(2008:size(larct_f60,1),i));
loom3rd_moment_f60(1,i)=loom3rd_moment_f60_temp;
end
mean_f60_loom3rd=mean(loom3rd_moment_f60);


base_and_loom_resp_f60=[];

base_and_loom_resp_f60(1,1)=sum(pre_moment_f60);
base_and_loom_resp_f60(2,1)=sum(loom1st_moment_f60);
base_and_loom_resp_f60(3,1)=sum(rest1_moment_f60);
base_and_loom_resp_f60(4,1)=sum(loom2nd_moment_f60);
base_and_loom_resp_f60(5,1)=sum(rest2_moment_f60);
base_and_loom_resp_f60(6,1)=sum(loom3rd_moment_f60);




loom_resp_ratio_f60=(base_and_loom_resp_f60/size(larct_f60,2));



figure;
plot(loom_resp_ratio_f60);
%hold on;
%plot(loom_resp_ratio_ctrl);


%%%% now to try to analyse it by loom

loom_times=[0,66,120,180,234,300,360,420,474,540]; %%% this is the times to add at the beggining of each periods to have the loom timepoints



%%%first loom moment
loom1st_moment_f60_perloom=[];
for i=1:size(larct_f60,2)
    start=300;
    for k=1:10
    loom1st_moment_f60_temp=sum(larct_f60(start+loom_times(k):start+loom_times(k)+3,i)); %%%% this is to get to each loom time point and sum the responses in that happend on the 5s of the loom stimuli 
    
    %%%and this loop is to only count as 1 response to the looms all the
    %%% fast movements
    if loom1st_moment_f60_temp>1
    loom1st_moment_f60_perloom(k,i)=1;
    else
    loom1st_moment_f60_perloom(k,i)=loom1st_moment_f60_temp;
    end
    
    end
end

%%%second loom moment
loom2nd_moment_f60_perloom=[];
for i=1:size(larct_f60,2)
    start=1154;
    for k=1:10
    loom2nd_moment_f60_temp=sum(larct_f60(start+loom_times(k):start+loom_times(k)+3,i));
    if loom2nd_moment_f60_temp>1
    loom2nd_moment_f60_perloom(k,i)=1;
    else
    loom2nd_moment_f60_perloom(k,i)=loom2nd_moment_f60_temp;
    end
    
    
    end
end


%%%third loom moment
loom3rd_moment_f60_perloom=[];
for i=1:size(larct_f60,2)
    start=2008;
    for k=1:10
    loom3rd_moment_f60_temp=sum(larct_f60(start+loom_times(k):start+loom_times(k)+3,i));
    if loom3rd_moment_f60_temp>1
    loom3rd_moment_f60_perloom(k,i)=1;
    else
    loom3rd_moment_f60_perloom(k,i)=loom3rd_moment_f60_temp;
    end
    
    
    end
end




sumF60=zeros(10,3);
for i=1:10
    sumF60(i,1)=sum(loom1st_moment_f60_perloom(i,:));
    sumF60(i,2)=sum(loom2nd_moment_f60_perloom(i,:));
    sumF60(i,3)=sum(loom3rd_moment_f60_perloom(i,:));
end
ratioF60=sumF60/size(larct_f60,2);






%%%% and to look for the sound response!


sound_f60=[];
sound_times=[0,5,10];
for i=1:size(larct_f60,2)
    start=1982;
    for k=1:3 %%%cause there are 3 sounds
    
    
    sound_f60_temp=sum(larct_f60(start+sound_times(k):start+sound_times(k)+2,i));
    if sound_f60_temp>1
    sound_f60(k,i)=1;
    else
    sound_f60(k,i)=sound_f60_temp;
    end
    end
end
mean_sound_f60=mean(sound_f60');



All_perloom_f60=vertcat(loom1st_moment_f60_perloom,loom2nd_moment_f60_perloom,loom3rd_moment_f60_perloom);


%%


%%%%now for s20 fish

larct_s20=zeros(length(slooms_var20ISI)-1,1);


for i=1:length(larct_s20)

    larct_s20(i,1)=slooms_var20ISI{i+1,11}; %%%the + 1 is to get rid of the first row wich has the names
end

Tot_sum_s20=sum(larct_s20);
Tot_mean_s20=mean(larct_s20);


%%%this is to get rid of the extra fast movements (>30mm/s)  that happend on each timepoints (of 1s)
%%% this means that the data I am analyzing from now on is taking into
%%% acount the times a fast movement happend in a second, not the total
%%% number of fast movements. 
larct_s20_new=zeros(length(slooms_var20ISI)-1,1);
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
larct_s20 = vec2mat(larct_s20_new, max([slooms_var20ISI{:,2}]));
larct_s20 = larct_s20';






%%%% now I will try to analyze how many fast responses where in the pre
%%%% and rest moments. And in the 3 loom periods. 

%%% the pre moment

pre_moment_s20=[];
for i=1:size(larct_s20,2)
pre_moment_s20_temp=sum(larct_s20(1:299,i));
pre_moment_s20(1,i)=pre_moment_s20_temp;
end

mean_s20_pre=mean(pre_moment_s20);

%%% the rest moments

rest1_moment_s20=[];
for i=1:size(larct_s20,2)
rest1_moment_s20_temp=sum(larct_s20(496:795,i));
rest1_moment_s20(1,i)=rest1_moment_s20_temp;
end
mean_s20_rest1=mean(rest1_moment_s20);

%%% the rest moment

rest2_moment_s20=[];
for i=1:size(larct_s20,2)
rest2_moment_s20_temp=sum(larct_s20(990:1289,i));
rest2_moment_s20(1,i)=rest2_moment_s20_temp;
end
mean_s20_rest2=mean(rest2_moment_s20);



%%%first loom moment
loom1st_moment_s20=[];
for i=1:size(larct_s20,2)
loom1st_moment_s20_temp=sum(larct_s20(300:495,i));
loom1st_moment_s20(1,i)=loom1st_moment_s20_temp;
end

mean_s20_loom1st=mean(loom1st_moment_s20);

%%%second loom moment
loom2nd_moment_s20=[];
for i=1:size(larct_s20,2)
loom2nd_moment_s20_temp=sum(larct_s20(796:989,i));
loom2nd_moment_s20(1,i)=loom2nd_moment_s20_temp;
end

mean_s20_loom2nd=mean(loom2nd_moment_s20);

%%%third loom moment
loom3rd_moment_s20=[];
for i=1:size(larct_s20,2)
loom3rd_moment_s20_temp=sum(larct_s20(1290:1484,i));
loom3rd_moment_s20(1,i)=loom3rd_moment_s20_temp;
end
mean_s20_loom3rd=mean(loom3rd_moment_s20);


base_and_loom_resp_s20=[];

base_and_loom_resp_s20(1,1)=sum(pre_moment_s20);
base_and_loom_resp_s20(2,1)=sum(loom1st_moment_s20);
base_and_loom_resp_s20(3,1)=sum(rest1_moment_s20);
base_and_loom_resp_s20(4,1)=sum(loom2nd_moment_s20);
base_and_loom_resp_s20(5,1)=sum(rest2_moment_s20);
base_and_loom_resp_s20(6,1)=sum(loom3rd_moment_s20);




loom_resp_ratio_s20=(base_and_loom_resp_s20/size(larct_s20,2));



figure;
plot(loom_resp_ratio_s20);
%hold on;
%plot(loom_resp_ratio_ctrl);


%%%% now to try to analyse it by loom

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

%%%second loom moment
loom2nd_moment_s20_perloom=[];
for i=1:size(larct_s20,2)
    start=796;
    for k=1:10
    loom2nd_moment_s20_temp=sum(larct_s20(start+loom_times(k):start+loom_times(k)+5,i));
    if loom2nd_moment_s20_temp>1
    loom2nd_moment_s20_perloom(k,i)=1;
    else
    loom2nd_moment_s20_perloom(k,i)=loom2nd_moment_s20_temp;
    end
    
    
    end
end


%%%third loom moment
loom3rd_moment_s20_perloom=[];
for i=1:size(larct_s20,2)
    start=1290;
    for k=1:10
    loom3rd_moment_s20_temp=sum(larct_s20(start+loom_times(k):start+loom_times(k)+5,i));
    if loom3rd_moment_s20_temp>1
    loom3rd_moment_s20_perloom(k,i)=1;
    else
    loom3rd_moment_s20_perloom(k,i)=loom3rd_moment_s20_temp;
    end
    
    
    end
end




sumS20=zeros(10,3);
for i=1:10
    sumS20(i,1)=sum(loom1st_moment_s20_perloom(i,:));
    sumS20(i,2)=sum(loom2nd_moment_s20_perloom(i,:));
    sumS20(i,3)=sum(loom3rd_moment_s20_perloom(i,:));
end
ratioS20=sumS20/size(larct_s20,2);



%%%% and to look for the sound response!


sound_s20=[];
sound_times=[0,5,10];
for i=1:size(larct_s20,2)
    start=1264;
    for k=1:3 %%%cause there are 3 sounds
    
    
    sound_s20_temp=sum(larct_s20(start+sound_times(k):start+sound_times(k)+2,i));
    if sound_s20_temp>1
    sound_s20(k,i)=1;
    else
    sound_s20(k,i)=sound_s20_temp;
    end
    end
end
mean_sound_s20=mean(sound_s20');



All_perloom_s20=vertcat(loom1st_moment_s20_perloom,loom2nd_moment_s20_perloom,loom3rd_moment_s20_perloom);


%%

%%%%now the s60 fish


larct_s60=zeros(length(slooms_var60ISI)-1,1);


for i=1:length(larct_s60)

    larct_s60(i,1)=slooms_var60ISI{i+1,11}; %%%the + 1 is to get rid of the first row wich has the names
end

Tot_sum_s60=sum(larct_s60);
Tot_mean_s60=mean(larct_s60);


%%%this is to get rid of the extra fast movements (>30mm/s)  that happend on each timepoints (of 1s)
%%% this means that the data I am analyzing from now on is taking into
%%% acount the times a fast movement happend in a second, not the total
%%% number of fast movements. 
larct_s60_new=zeros(length(slooms_var60ISI)-1,1);
for i=1:length(larct_s60)
    if larct_s60(i,1)>1
    larct_s60_new(i,1)=1;
    
    else
        larct_s60_new(i,1)=larct_s60(i,1);
    
    end
end

Tot_sum_s60_new=sum(larct_s60_new);
Tot_mean_s60_new=mean(larct_s60_new);



%%%his is to stack them by fish. I need to take into account the number of
%%%bins. 
larct_s60 = vec2mat(larct_s60_new, max([slooms_var60ISI{:,2}]));
larct_s60 = larct_s60';






%%%% now I will try to analyze how many fast responses where in the pre
%%%% and rest moments. And in the 3 loom periods. 

%%% the pre moment

pre_moment_s60=[];
for i=1:size(larct_s60,2)
pre_moment_s60_temp=sum(larct_s60(1:299,i));
pre_moment_s60(1,i)=pre_moment_s60_temp;
end

mean_s60_pre=mean(pre_moment_s60);

%%% the rest moments

rest1_moment_s60=[];
for i=1:size(larct_s60,2)
rest1_moment_s60_temp=sum(larct_s60(856:1155,i));
rest1_moment_s60(1,i)=rest1_moment_s60_temp;
end
mean_s60_rest1=mean(rest1_moment_s60);

%%% the rest moment

rest2_moment_s60=[];
for i=1:size(larct_s60,2)
rest2_moment_s60_temp=sum(larct_s60(1710:2009,i));
rest2_moment_s60(1,i)=rest2_moment_s60_temp;
end
mean_s60_rest2=mean(rest2_moment_s60);



%%%first loom moment
loom1st_moment_s60=[];
for i=1:size(larct_s60,2)
loom1st_moment_s60_temp=sum(larct_s60(300:855,i));
loom1st_moment_s60(1,i)=loom1st_moment_s60_temp;
end

mean_s60_loom1st=mean(loom1st_moment_s60);

%%%second loom moment
loom2nd_moment_s60=[];
for i=1:size(larct_s60,2)
loom2nd_moment_s60_temp=sum(larct_s60(1156:1709,i));
loom2nd_moment_s60(1,i)=loom2nd_moment_s60_temp;
end

mean_s60_loom2nd=mean(loom2nd_moment_s60);

%%%third loom moment
loom3rd_moment_s60=[];
for i=1:size(larct_s60,2)
loom3rd_moment_s60_temp=sum(larct_s60(2010:size(larct_s60,1),i));
loom3rd_moment_s60(1,i)=loom3rd_moment_s60_temp;
end
mean_s60_loom3rd=mean(loom3rd_moment_s60);


base_and_loom_resp_s60=[];

base_and_loom_resp_s60(1,1)=sum(pre_moment_s60);
base_and_loom_resp_s60(2,1)=sum(loom1st_moment_s60);
base_and_loom_resp_s60(3,1)=sum(rest1_moment_s60);
base_and_loom_resp_s60(4,1)=sum(loom2nd_moment_s60);
base_and_loom_resp_s60(5,1)=sum(rest2_moment_s60);
base_and_loom_resp_s60(6,1)=sum(loom3rd_moment_s60);




loom_resp_ratio_s60=(base_and_loom_resp_s60/size(larct_s60,2));



figure;
plot(loom_resp_ratio_s60);
%hold on;
%plot(loom_resp_ratio_ctrl);


%%%% now to try to analyse it by loom

loom_times=[0,66,120,180,234,300,360,420,474,540]; %%% this is the times to add at the beggining of each periods to have the loom timepoints



%%%first loom moment
loom1st_moment_s60_perloom=[];
for i=1:size(larct_s60,2)
    start=300;
    for k=1:10
    loom1st_moment_s60_temp=sum(larct_s60(start+loom_times(k):start+loom_times(k)+5,i)); %%%% this is to get to each loom time point and sum the responses in that happend on the 5s of the loom stimuli 
    
    %%%and this loop is to only count as 1 response to the looms all the
    %%% fast movements
    if loom1st_moment_s60_temp>1
    loom1st_moment_s60_perloom(k,i)=1;
    else
    loom1st_moment_s60_perloom(k,i)=loom1st_moment_s60_temp;
    end
    
    end
end

%%%second loom moment
loom2nd_moment_s60_perloom=[];
for i=1:size(larct_s60,2)
    start=1156;
    for k=1:10
    loom2nd_moment_s60_temp=sum(larct_s60(start+loom_times(k):start+loom_times(k)+5,i));
    if loom2nd_moment_s60_temp>1
    loom2nd_moment_s60_perloom(k,i)=1;
    else
    loom2nd_moment_s60_perloom(k,i)=loom2nd_moment_s60_temp;
    end
    
    
    end
end


%%%third loom moment
loom3rd_moment_s60_perloom=[];
for i=1:size(larct_s60,2)
    start=2010;
    for k=1:10
    loom3rd_moment_s60_temp=sum(larct_s60(start+loom_times(k):start+loom_times(k)+5,i));
    if loom3rd_moment_s60_temp>1
    loom3rd_moment_s60_perloom(k,i)=1;
    else
    loom3rd_moment_s60_perloom(k,i)=loom3rd_moment_s60_temp;
    end
    
    
    end
end




sumS60=zeros(10,3);
for i=1:10
    sumS60(i,1)=sum(loom1st_moment_s60_perloom(i,:));
    sumS60(i,2)=sum(loom2nd_moment_s60_perloom(i,:));
    sumS60(i,3)=sum(loom3rd_moment_s60_perloom(i,:));
end
ratioF60=sumS60/size(larct_s60,2);






%%%% and to look for the sound response!


sound_s60=[];
sound_times=[0,5,10];
for i=1:size(larct_s60,2)
    start=1984;
    for k=1:3 %%%cause there are 3 sounds
    
    
    sound_s60_temp=sum(larct_s60(start+sound_times(k):start+sound_times(k)+2,i));
    if sound_s60_temp>1
    sound_s60(k,i)=1;
    else
    sound_s60(k,i)=sound_s60_temp;
    end
    end
end
mean_sound_s60=mean(sound_s60');



All_perloom_s60=vertcat(loom1st_moment_s60_perloom,loom2nd_moment_s60_perloom,loom3rd_moment_s60_perloom);


%%



%%% a bit of general analysis

%%%total responses per periods
figure;
hold on;
plot(loom_resp_ratio_f20);
plot(loom_resp_ratio_f60);
plot(loom_resp_ratio_s20);
plot(loom_resp_ratio_s60);
hold off;

%%%per loom
figure;
hold on;
plot(mean(All_perloom_f20'));
plot(mean(All_perloom_f60'));
plot(mean(All_perloom_s20'));
plot(mean(All_perloom_s60'));
hold off;


%%%per block

blocks_f20=horzcat(mean(mean(loom1st_moment_f20_perloom)),mean(mean(loom2nd_moment_f20_perloom)),mean(mean(loom3rd_moment_f20_perloom)));
blocks_f60=horzcat(mean(mean(loom1st_moment_f60_perloom)),mean(mean(loom2nd_moment_f60_perloom)),mean(mean(loom3rd_moment_f60_perloom)));
blocks_s20=horzcat(mean(mean(loom1st_moment_s20_perloom)),mean(mean(loom2nd_moment_s20_perloom)),mean(mean(loom3rd_moment_s20_perloom)));
blocks_s60=horzcat(mean(mean(loom1st_moment_s60_perloom)),mean(mean(loom2nd_moment_s60_perloom)),mean(mean(loom3rd_moment_s60_perloom)));

figure;
hold on;
plot(blocks_f20);
plot(blocks_f60);
plot(blocks_s20);
plot(blocks_s60);
hold off;


%%%the sound
figure;
hold on;
plot(mean(sound_f20'));
plot(mean(sound_f60'));
plot(mean(sound_s20'));
plot(mean(sound_s60'));
hold off;


