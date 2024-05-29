
clear;clc;
%% Figure4: connectivity distance calculation
addpath(genpath('/media/shuang/data/Covariance/packages/BrainSpace-0.1.2'));
addpath(genpath('/media/shuang/data/Covariance/packages/cifti-matlab-master'));
addpath(genpath('/media/shuang/data/Covariance/packages/BrewerMap-master'));
addpath(genpath('/media/shuang/data/Covariance/packages/brainstat_matlab'));
addpath(genpath('/media/shuang/data/Covariance/packages/ENIGMA_matlab'));
%% add surfaces and parcellations for plotting
[surf_lh, surf_rh] = load_conte69();
path_wb_command = ['/media/shuang/data/Covariance/packages/workbench/bin_linux64/wb_command'];
labeling2 = ciftiopen('/media/shuang/data/Covariance/template/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
labeling = labeling2.cdata;

%% load data
load('/media/shuang/data/Covariance/results/conn_HC_cov.mat');
load('/media/shuang/data/Covariance/results/conn_EOS_cov.mat');
load('/media/shuang/data/Covariance/data/GD400.mat');
load('/media/shuang/data/Covariance/results/Gradients.mat'); % own 400*10*2(EOS;HC) gradients
GD = GD400;
GD(201:end,1:200) = (GD400(1:200,1:200)+GD400(201:end,201:end))/2;
GD(1:200,201:end) = GD(201:end,1:200); % average l+r distance for interhemi
figure; heatmap(GD,'ColorLimits',[0 235],'GridVisible','off'); colormap(brewermap([],"BuPu"))

%------------plot node-wise geodesic distance degree-------------
obj=plot_hemispheres(mean(GD)',{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','Geodesic distance degree');
obj.colormaps(brewermap([],"BuPu"))

%% correlation with sc gradients
[gdr,gdp]=corr(mean(GD)',Gradient(:,2,2))
[rr,pp]=corr(reshape(conn_HC_cov,400*400,1),reshape(GD,400*400,1));
[spin_p,~]=spin_test(mean(GD)',Gradient(:,1,2),'surface_name','fsa5',...
    'parcellation_name','schaefer_400','n_rot',10000,'type','Pearson');

nperm = 10000;
nlen = 400*400;
x = reshape(conn_HC_cov,400*400,1); y = reshape(GD,400*400,1);
rr_null = zeros(nperm,1);
for m = 1:nperm % for each permutation
    clear x_rand
    x_rand = x(randperm(nlen));
    rr_null(m,1) = corr(x_rand,y); % null loadings
end
sperm = (1+(nnz(find(abs(rr_null)>=abs(rr)))))/(1+nperm); % permutation test

%% Connectivity distance: Threshold covariance matrix and sum geodesic connectivity
thre11 = ones(400,400);thre12 = ones(400,400);
[~,index1(:,:,1)] = sort(conn_HC_cov,'descend');[~,index1(:,:,2)] = sort(conn_EOS_cov,'descend');
per = 10; num = 400.*per/100; % localscn is relative stable for this parameter: TOP 10%
for i = 1:size(GD,1)
    thre11(index1(num+1:end,i,1),i) = 0;
    thre12(index1(num+1:end,i,2),i) = 0;
    CD_HC(1,i) = mean(GD(logical(thre11(:,i)),i)); % connectivity distance
    CD_EOS(1,i) = mean(GD(logical(thre12(:,i)),i));
end
CD=[CD_HC;CD_EOS]';

obj=plot_hemispheres(CD(:,1),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','HC connectivity distance');
obj.colormaps(brewermap([],"BuPu"))
obj.colorlimits([55, 133]) % HC
obj=plot_hemispheres(CD(:,2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','EOS connectivity distance');
obj.colormaps(brewermap([],"BuPu"))
obj.colorlimits([55, 133]) % EOS

%% z-test for covariance distance
z=zeros(400,1);zp=zeros(400,1);
CD_HC=CD(:,1);CD_EOS=CD(:,2);
CD_HC=zscore(CD_HC);CD_EOS=zscore(CD_EOS);
for i = 1:400
    clear t_r1 t_r2
    t_r1 = CD_EOS(i);
    t_r2 = CD_HC(i);
    z(i) = t_r1-t_r2;
    zp(i) = (1-normcdf(abs(z(i)),0,1))*2;
end
fdrz = mafdr(zp,'BHFDR',true);

%%%%% save as dlabel.nii mask file
mask=ciftiopen('/media/shuang/data/Covariance/template/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
for i = 1:400 % fs_LR 32k space for saving dlabel.nii file
    mask.cdata(find(mask.cdata==i))=logical(zp(i)<(0.05));
end
ciftisavereset(mask,'/media/shuang/data/Covariance/results/CDg2z05.dlabel.nii',path_wb_command);

mask=ciftiopen('/media/shuang/data/Covariance/template/SCON003_400.pscalar.nii',path_wb_command);
mask.cdata=z;
ciftisavereset(mask,'/media/shuang/data/Covariance/results/CDg2z.pscalar.nii',path_wb_command);
%%% end

%% statistical analyses
clc;fprintf('HC CD = %f +- %f, EOS CD = %f +- %f', mean(CD_HC), std(CD_HC), mean(CD_EOS), std(CD_EOS));
BIN=10; % divide into 10 bins
[~,index1]=sort(Gradient(:,2,2)); % HC small to big index; inverse G2
[~,rank1]=sort(index1); % position for each element
BinsHC=reshape(index1,400/BIN,BIN); % index for 10 bins
maskHC=zeros(400,1);
clear cdp cdt
for i = 1:BIN
    maskHC(BinsHC(:,i))=i;
    [~,cdp(i),~,cdt{i}]=ttest(CD_EOS(BinsHC(:,i)),CD_HC(BinsHC(:,i)));
end
FDRvalue = mafdr(cdp,'BHFDR',true);
