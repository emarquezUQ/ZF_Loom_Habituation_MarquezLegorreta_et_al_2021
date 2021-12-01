


%%%% this script is to do the cross-validation test. I will take one fish
%%%% out and regenerate the mean corr-matrices and show that the general
%%%% results still hold. now with the fmr1 dataset

%%% part of the script is based on how gilles did it for the fmr1 dataset

%%
load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

datasets=['f20'; 'f60'; 's20'; 's60'];


keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;



%%

cross_val_fmr1=struct;
for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));
CorrMatrices_mean2=zeros(21,length(names)-1,90,90);
for loom=1:21        
    for fish_rem_nb=1:length(names)-1
        temp=nan(length(names)-2,90,90);    
        counter=1;
        for fish_nb=1:length(names)-1
            if fish_nb ~= fish_rem_nb
                fish_name=names(fish_nb);
                temp(counter,:,:)=Data_corrMat4.(groupnames{group}).(fish_name{1}).loomsR{1,loom}(keepFmr1,keepFmr1);
                counter=counter+1;            
            end
        end
        CorrMatrices_mean2(loom,fish_rem_nb,:,:)=squeeze(nanmean(temp,1));
        
    end       
end
cross_val_fmr1.(groupnames{group}).CorrMatrices_mean2=CorrMatrices_mean2;
end


for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));

%%% setting the diagonal to 0. 
%%% looms 1,2,3,10 and 11
MatAll_corrected_crossval=[];

for fish=1:length(names)-1

    moment=[1 2 3 4 5 6 11 12]; %%% looms 1,2,3,4,5,10 and 11 (cause 1 is pre loom)
    %loom=[1 2 3 4 5 10 11];
     for m=1:length(moment)
     Mat = threshold_absolute(abs(squeeze(squeeze(cross_val_fmr1.(groupnames{group}).CorrMatrices_mean2(moment(m),fish,:,:)))),0.75);
     MatAll_corrected_crossval(m,fish,:,:)=Mat;
     end

end
cross_val_fmr1.(groupnames{group}).MatAll_corrected_crossval=MatAll_corrected_crossval;

end

%%


%%%% getting the number of high corr >0.75

for group=1:3

names = fieldnames(Data_corrMat4.(groupnames{group}));

crossval_highcorr=[];
for fish=1:length(names)-1
for m=1:length(moment)

temp_highcorr=length(find(squeeze(squeeze(cross_val_fmr1.(groupnames{group}).MatAll_corrected_crossval(m,fish,:,:)))));

crossval_highcorr(m,fish)=temp_highcorr;
end
end
cross_val_fmr1.(groupnames{group}).crossval_highcorrNum=crossval_highcorr;
end


%%%% then I took the results to Prism to plot it and do stats. 
