
%%%% Making a dataset sample 

RegionList={'Pallium','Subpallium','Thalamus','Habenula','Pretectum','Tectum','Tegmentum','Cerebellum','Hindbrain'};


load('means_F20_CL4n7.mat','mean_CL4_f20','mean_CL7_f20');
load('means_S20_CL4n7.mat','mean_CL4_s20','mean_CL7_s20');
clustersF=fieldnames(mean_CL7_f20);
clustersS=fieldnames(mean_CL7_s20);


%%% getting the nodes data for f60 and s20
load('Nodes_N_means_alldatasets2.mat')

Nodes.Nod_coor=Nodes2.Mod_loc;
Nodes.Nod_clustID=Nodes2.Mod_clust;
Nodes.Nod_brainID=Nodes2.Mod_brain;

Nodes.f60.NodeMats=Nodes2.f60.mean_matrix_K;
Nodes.s20.NodeMats=Nodes2.s20.mean_matrix_K;


%%% nodes for the fmr1 dataset
load('NodesNgraphFmr1Loomhab4.mat')

NodesFmr1.Nod_coor=Nodes4.Mod_loc;
NodesFmr1.Nod_clustID=Nodes4.Mod_clust;
NodesFmr1.Nod_brainID=Nodes4.Mod_brain;

NodesFmr1.fmr1.NodeMats=Nodes4.mean_matrix.fmr1;
NodesFmr1.wt.NodeMats=Nodes4.mean_matrix.control;


save('Nodes_WTnFmr1Exp.mat','Nodes','NodesFmr1','Zbrain_brainMask2D','RegionList');

