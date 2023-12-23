
clear;clc;
%% Figure4: connectivity distance
addpath(genpath('/media/shuang/data/Covariance/packages/BrainSpace-0.1.2'));
addpath(genpath('/media/shuang/data/Covariance/packages/cifti-matlab-master'));
addpath(genpath('/media/shuang/data/Covariance/packages/BrewerMap-master'));
addpath(genpath('/media/shuang/data/Covariance/packages/brainstat_matlab'));

%% add surfaces and parcellations for plotting
[surf_lh, surf_rh] = load_conte69();
path_wb_command = ['/media/shuang/data/repository/matlab_packages/workbench/bin_linux64/wb_command'];
labeling2 = ciftiopen('/media/shuang/data/repository/templates/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
labeling = labeling2.cdata;

%% load data
load('/media/shuang/data/Covariance/results/conn_HC_cov.mat');
load('/media/shuang/data/Covariance/results/conn_EOS_cov.mat');
load('/media/shuang/data/Covariance/data/GD400.mat');
GD = GD400;
GD(201:end,1:200) = (GD400(1:200,1:200)+GD400(201:end,201:end))/2;
GD(1:200,201:end) = GD(201:end,1:200); % average l+r distance for interhemi
figure; heatmap(GD,'ColorLimits',[0 235],'GridVisible','off'); colormap(brewermap([],"BuPu"))

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
             'labeltext','HC connectivityD');
obj.colormaps(brewermap([],"BuPu"))
obj.colorlimits([55, 133]) % HC

obj=plot_hemispheres(CD(:,2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','EOS connectivityD');
obj.colormaps(brewermap([],"BuPu"))
obj.colorlimits([55, 133]) % HC

% mask=ciftiopen('/media/shuang/data/Covariance/template/SCON003_400.pscalar.nii',path_wb_command);
% mask.cdata=CD_HC';
% ciftisavereset(mask,'/media/shuang/data/Covariance/results/HCCD.pscalar.nii',path_wb_command);
% mask.cdata=CD_EOS';
% ciftisavereset(mask,'/media/shuang/data/Covariance/results/EOSCD.pscalar.nii',path_wb_command);

%% statistical analyses
fprintf('HC CD = %f +- %f, EOS CD = %f +- %f', mean(CD_HC), std(CD_HC), mean(CD_EOS), std(CD_EOS));
load('/media/shuang/data/Covariance/results/GradientEOSHC.mat'); % own 400*10*2(EOS;HC) gradients
BIN=10; % divide into 10 bins
[~,index1]=sort(-Gradient(:,1,2)); % small to big index; inverse G1
[~,rank1]=sort(index1); % position for each element
BinsHC=reshape(index1,400/BIN,BIN); % index for 10 bins
maskHC=zeros(400,1);
clear cdp cdt
for i = 1:BIN
    maskHC(BinsHC(:,i))=i;
    [~,cdp(i),~,cdt{i}]=ttest2(CD_EOS(BinsHC(:,i)),CD_HC(BinsHC(:,i)));
end
% maskHC1=zeros(400,1);
% maskHC1(find(maskHC==10))=10;
% maskHC1(find(maskHC==1))=1;
% maskHC1(find(maskHC==2))=2;
% maskHC1=logical(maskHC1);
% mask=ciftiopen('/media/shuang/data/p_02658/AAAAAscripts/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
% for i = 1:400 % fs_LR 32k space for saving dlabel.nii file
%     mask.cdata(find(mask.cdata==i))=maskHC1(i);
% end
% ciftisavereset(mask,'/media/shuang/data/Covariance/results/CDBin.dlabel.nii',path_wb_command);
