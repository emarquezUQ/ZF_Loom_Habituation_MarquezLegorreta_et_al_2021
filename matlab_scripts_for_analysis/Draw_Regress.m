function [SelectedSamples] = Draw_Regress(ModelResults,DataSet,idxKmeans,rsq_limit)

x = linspace(1,size(DataSet,2),size(DataSet,2));
SelectedSamples=find([ModelResults.rsquared]>rsq_limit);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(length(SelectedSamples)));yplot=ceil(length(SelectedSamples)/xplot);
for i=SelectedSamples    
    NumberOfCells=length(find(idxKmeans==i));
    subplot(xplot,yplot,counter);plot(x,DataSet(i,:),x,ModelResults(i).Fitted);title(num2str(NumberOfCells))
    xlim([0 size(DataSet,2)])
    counter=counter+1;
end