

%%%%% this script is to get a range of gammas and omegas results. 


%%%%% range of gammas (0-2.5) and omegas (0-2) and the
%%%%% repetitions of the community detection will be 100. 

load('Nodes_WTnFmr1Exp.mat')
load('Nodes_FnS_all_samples.mat')
load('graphs_fmr1Exp.mat')
load('graphs_FnS_all.mat')

keepFmr1=load('graphs_fmr1Exp.mat','keep');
keepFmr1=keepFmr1.keep;

%%

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



%%

Big_FMR1_OPT={};
FlexMat=struct;
CoheMat=struct;
PromMat=struct;


%%
Allgamma = 0.1:0.1:2.5;
Allomega = 0.1:0.1:2;


%for g=1%:16
for g=1:25
for o=1:20

gamma = Allgamma(g);
omega = Allomega(o);   
    
%% making 100 iterations per group to calculate a representative modularity 
S_cons=struct;
%%
for group=1:3
A={};

for loom=1:21
  
   temp=Data_corrMat4.(groupnames{group}).Mean_corrMat{1,loom}(keepFmr1,keepFmr1);
   temp(isnan(temp))=0;
    A{loom}=temp;  
end
clear temp

%%

%%%% making a multidimensional structure where to store things
S_test=[];

for test=1:100
     

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B,10000,0);
%[S,Q,nb_it] = iterated_genlouvain(B);
Q = Q/twomu;
S = reshape(S,N,T);

S_test=cat(3,S_test,S);
  
  
end


%% trying consensus_iterative 


S_good=[];
for i=1:21
   C=S_test(:,i,:); 
   C=squeeze(C);
   [S2, Q2, X_new3, qpc] = consensus_iterative(C');
    S_good(:,i)=S2(i,:);
end

S_cons.(groupnames{group}).S_test=S_test;
S_cons.(groupnames{group}).S_cons=S_good;



%% testing flexibility

flex=flexibility(S_cons.(groupnames{group}).S_cons','temp'); %%% need to mind the orientation of the matrix to do it properly

S_cons.(groupnames{group}).flex=flex;

%% cohesion strength and related

options.figureFlag	= 0;
options.colormap	= 'jet';

[Cij,node_cohesion,node_disjoint,node_flexibility,strength_cohesion,commChanges,commCohesion,commDisjoint,commIndex] = calc_node_cohesion(S_cons.(groupnames{group}).S_cons,options);

S_cons.(groupnames{group}).node_cohesion=node_cohesion;

%% promiscuity

P = promiscuity(S_cons.(groupnames{group}).S_cons');  %%% need to mind the orientation of the matrix to do it properly

S_cons.(groupnames{group}).P=P;

end

%%
groups_flexibility=[];
counter=1;
for group=[3 1 2]

    groups_flexibility(:,counter)=S_cons.(groupnames{group}).flex;

counter=counter+1;
end

groups_cohesion=[];
counter=1;
for group=[3 1 2]

    groups_cohesion(:,counter)=S_cons.(groupnames{group}).node_cohesion;

counter=counter+1;
end

groups_prom=[];
counter=1;
for group=[3 1 2]

    groups_prom(:,counter)=S_cons.(groupnames{group}).P;

counter=counter+1;
end


%% relative change


%%%%%% flexibility

flex_dif=S_cons.(groupnames{3}).flex-S_cons.(groupnames{2}).flex;

%%%%%% cohesion

cohe_dif=S_cons.(groupnames{3}).node_cohesion-S_cons.(groupnames{2}).node_cohesion;

%%%%%% promiscuity
P_dif=S_cons.(groupnames{3}).P-S_cons.(groupnames{2}).P;

%%%%% 

FlexMat.(groupnames{3})(g,o)=mean(S_cons.(groupnames{3}).flex);
FlexMat.(groupnames{1})(g,o)=mean(S_cons.(groupnames{1}).flex);
FlexMat.(groupnames{2})(g,o)=mean(S_cons.(groupnames{2}).flex);

CoheMat.(groupnames{3})(g,o)=mean(S_cons.(groupnames{3}).node_cohesion);
CoheMat.(groupnames{1})(g,o)=mean(S_cons.(groupnames{1}).node_cohesion);
CoheMat.(groupnames{2})(g,o)=mean(S_cons.(groupnames{2}).node_cohesion);

PromMat.(groupnames{3})(g,o)=mean(S_cons.(groupnames{3}).P);
PromMat.(groupnames{1})(g,o)=mean(S_cons.(groupnames{1}).P);
PromMat.(groupnames{2})(g,o)=mean(S_cons.(groupnames{2}).P);


%%
name=strcat('g',num2str(10*gamma),'_o',num2str(10*omega));

Big_FMR1_OPT{g,o}=S_cons;


end
end


%% saving

save('S_cons_testing_GnO_measures6.mat','Big_FMR1_OPT','FlexMat','CoheMat','PromMat','Allgamma','Allomega');



%%

counter=1;
figure;
for group=[3 1 2]
   subplot(3,3,counter);imagesc(FlexMat.(groupnames{group}));caxis([0 1]);
   title(strcat('flex/',groupnames{group}));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);imagesc(CoheMat.(groupnames{group}));caxis([0 1]);
   title(strcat('cohe/',groupnames{group}));
   
   counter=counter+1;
end
for group=[3 1 2]
   subplot(3,3,counter);imagesc(PromMat.(groupnames{group}));caxis([0 1]);
   title(strcat('prom/',groupnames{group}));
   
   counter=counter+1;
end




