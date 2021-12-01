


%%%%% this script is to do a temporal null model for the whole fmr1 dataset
%%%%% with a temporal null model, but I am doing for values of gamma up to
%%%%% 2.5 and omega up to 2

load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;



%%
QoptFMR1_All=struct;
randOrder2 = randperm(21,21);
for group=1:3
A={};

for loom=1:21
  
   temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,loom}(keepFmr1,keepFmr1);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp

%%
Qopt2=[];

gamma = 0.1:0.1:2.5;
omega = 0.1:0.1:2;

for test=1:100
Qmat2=[];
for i=1:length(gamma)   
%parfor j=1:length(omega)
for j=1:length(omega)
    
N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma(i)*k'*k/twom;
end
twomu=twomu+2*omega(j)*N*(T-1);
B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
%S = reshape(S,N,T);

Qmat2(i,j)=Q;
   
end
end

Qopt2=cat(3,Qopt2,Qmat2);

end

meanQopt2=mean(Qopt2,3);
varQopt2=var(Qopt2,0,3);

figure;
subplot(1,2,1);imagesc(meanQopt2); colorbar;
subplot(1,2,2);imagesc(varQopt2); colorbar;


%% temporal null model




%unique(randOrder2)


At_null={};

for loom=1:21
    temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,randOrder2(loom)}(keepFmr1,keepFmr1);  
   temp(isnan(temp))=0;
    At_null{loom}=temp;  
end
clear temp



Qopt2_null2=[];

gamma = 0.1:0.1:2.5;
omega = 0.1:0.1:2;

for test=1:100
Qmat2_null2=[];
for i=1:length(gamma)   
%parfor j=1:length(omega)
for j=1:length(omega)

N=length(At_null{1});
T=length(At_null);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(At_null{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=At_null{s}-gamma(i)*k'*k/twom;
end
twomu=twomu+2*omega(j)*N*(T-1);
B = B + omega(j)*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
%S = reshape(S,N,T);

Qmat2_null2(i,j)=Q;
   
end
end

Qopt2_null2=cat(3,Qopt2_null2,Qmat2_null2);

end

meanQopt2_null2=mean(Qopt2_null2,3);
varQopt2_null2=var(Qopt2_null2,0,3);

figure;
subplot(1,2,1);imagesc(meanQopt2_null2); colorbar;
subplot(1,2,2);imagesc(varQopt2_null2); colorbar;


figure;imagesc(meanQopt2-meanQopt2_null2);colorbar;


QoptFMR1_All

QoptFMR1_All.(groupnames{group}).A=A;
QoptFMR1_All.(groupnames{group}).Qopt2=Qopt2;
QoptFMR1_All.(groupnames{group}).meanQopt2=meanQopt2;
QoptFMR1_All.(groupnames{group}).varQopt2=varQopt2;

QoptFMR1_All.(groupnames{group}).At_null=At_null;
QoptFMR1_All.(groupnames{group}).Qopt2_null2=Qopt2_null2;
QoptFMR1_All.(groupnames{group}).meanQopt2_null2=meanQopt2_null2;
QoptFMR1_All.(groupnames{group}).varQopt2_null2=varQopt2_null2;
end


save('QoptFMR1_All2.mat','QoptFMR1_All');



