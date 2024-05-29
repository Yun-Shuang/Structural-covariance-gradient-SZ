
clear;clc;
%% Figure1: structural covariance gradients calculation
addpath(genpath('/media/shuang/data/Covariance/packages/BrainSpace-0.1.2'));
addpath(genpath('/media/shuang/data/Covariance/packages/cifti-matlab-master'));
addpath(genpath('/media/shuang/data/Covariance/packages/BrewerMap-master'));
addpath(genpath('/media/shuang/data/Covariance/packages/brainstat_matlab'));

%% add surfaces and parcellations for plotting
[surf_lh, surf_rh] = load_conte69();
path_wb_command = ['/media/shuang/data/Covariance/packages/workbench/bin_linux64/wb_command'];
labeling2 = ciftiopen('/media/shuang/data/Covariance/template/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
labeling = labeling2.cdata;

%% load cortical thickness and covariates
load('/media/shuang/data/Covariance/data/CT400.mat'); % cortical thickness
load('/media/shuang/data/Covariance/data/regressors.mat'); %1:95 EOS; 96:end HC
eos = [1:95]; hc = [96:194];
HC = CT400(:,96:end); EOS = CT400(:,1:95);
% regressors(:,3)=[]; % discard mean cortical thickness
clc;fprintf('HC meanCT = %f +- %f, EOS meanCT = %f +- %f ', mean(HC(:)), std(HC(:)), mean(EOS(:)), std(EOS(:)));
[~,p,~,t]=ttest2(mean(HC,1),mean(EOS,1));

%% show an example of structural covariance
figure; plot(CT400(149,hc),CT400(91,hc),'k.','markersize',25);ylim([1.4 3.4]);xlim([2.2 3.8]);
hold on; plot(CT400(149,hc),polyval(polyfit(CT400(149,hc),CT400(91,hc),1),CT400(149,hc)),'r')
figure; plot(CT400(149,eos),CT400(91,eos),'k.','markersize',25);ylim([1.4 3.4]);xlim([2.2 3.8]);
hold on; plot(CT400(149,eos),polyval(polyfit(CT400(149,eos),CT400(91,eos),1),CT400(149,eos)),'r')

%% structural covariance matrix
[conn_HC_cov,conn_p(:,:,1)] = partialcorr(HC',regressors(hc,:)); % 
[conn_EOS_cov,conn_p(:,:,2)] = partialcorr(EOS',regressors(eos,:));
% ---------------- plot no threshold covariance matrix ---------------
for i=1:400
    for j=1:400
        if i~=j
    conn_HC_cov(i,j) = 0.5*log((1+conn_HC_cov(i,j))/(1-conn_HC_cov(i,j))); % fisher z-transform
    conn_EOS_cov(i,j) = 0.5*log((1+conn_EOS_cov(i,j))/(1-conn_EOS_cov(i,j))); % fisher z-transform
        else
    conn_HC_cov(i,j) = 0; 
    conn_EOS_cov(i,j) = 0; 
        end
    end
end
figure; heatmap(conn_HC_cov,'ColorLimits',[-0.5 0.5],'GridVisible','off'); colormap(flipud(brewermap([],"RdBu")))
figure; heatmap(conn_EOS_cov,'ColorLimits',[-0.5 0.5],'GridVisible','off'); colormap(flipud(brewermap([],"RdBu")))

% ---------------- plot thresholded covariance matrix ---------------
zz=conn_HC_cov;maskz=zeros(400,400);
[~,index] = sort(zz,'descend');
for i=1:400
    maskz(index(1:40,i),i)=1;
end
figure; heatmap(maskz,'ColorLimits',[0 1],'GridVisible','off'); colormap(brewermap([],"Greys")) % plot 10% threshold
zz=conn_EOS_cov;maskz=zeros(400,400);
[~,index] = sort(zz,'descend');
for i=1:400
    maskz(index(1:40,i),i)=1;
end
figure; heatmap(maskz,'ColorLimits',[0 1],'GridVisible','off'); colormap(brewermap([],"Greys"))
% ------------------ end ----------------------------------

%% Diffusion gradients

% ------------------- gradient calculating -----------------
Galign = GradientMaps('approach','diffusion embedding','alignment','pa','n_components',10);
Galign_HC = Galign.fit(conn_HC_cov);maskG=Galign_HC.gradients{1};
Galign_EOS = Galign.fit(conn_EOS_cov,'reference',maskG);
Gradient(:,:,1)=Galign_EOS.aligned{1};Gradient(:,:,2)=Galign_HC.gradients{1};
HC_G1=Gradient(:,1,2);EOS_G1=Gradient(:,1,1);
HC_G2=Gradient(:,2,2);EOS_G2=Gradient(:,2,1); %Gradient 2
% --------------------------- end -----------------------------
lambda=Galign_HC.lambda{1};
scree_plot(lambda(1:10)); % explained variance
clear lambda; lambda=Galign_EOS.lambda{1};
scree_plot(lambda(1:10)); % explained variance
% print mean +- SD
clc;fprintf('HC gradient1 = %f +- %f, EOS gradient1 = %f +- %f ', mean(HC_G1), std(HC_G1), mean(EOS_G1), std(EOS_G1));
fprintf('HC gradient2 = %f +- %f, EOS gradient2 = %f +- %f ', mean(HC_G2), std(HC_G2), mean(EOS_G2), std(EOS_G2));

% --------------------Gradient1 ------------------------------
obj=plot_hemispheres(HC_G1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','MaskG1');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"Greys"))])
obj.colorlimits([-0.04, 0.04]) 

obj=plot_hemispheres(EOS_G1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','EOS');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([-0.04, 0.04]) 
% ------------------ Gradient2 ---------------------------
obj=plot_hemispheres(HC_G2,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','maskg2');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([-0.05, 0.05]) 

obj=plot_hemispheres(EOS_G2,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','EOSg2');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([-0.05, 0.05]) 

%% compare with HCP gradient
addpath(genpath('/media/shuang/data/repository/matlab_packages/ENIGMA_matlab'));
maskG=load('/media/shuang/data/Covariance/template/strcov_gradient.csv'); %HCP 400*10 gradient mask
obj=plot_hemispheres(-maskG(:,1),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','HCPG1');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([-0.16, 0.16]) 
obj=plot_hemispheres(maskG(:,2),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','HCPG2');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([-0.16, 0.16]) 
[r(1,1),p(1,1)]=corr(rescale(HC_G1,-1,1),rescale(maskG(:,2),-1,1));
[r(1,2),p(1,2)]=corr(HC_G1,maskG(:,2));
[r(2,1),p(2,1)]=corr(HC_G2,-maskG(:,1));
[r(2,2),p(2,2)]=corr(HC_G2,maskG(:,2))

[spin_p(1,1),~]=spin_test(HC_G1,-maskG(:,1),'surface_name','fsa5',...
    'parcellation_name','schaefer_400','n_rot',10000,'type','Pearson');
[spin_p(1,2),~]=spin_test(HC_G1,maskG(:,2),'surface_name','fsa5',...
    'parcellation_name','schaefer_400','n_rot',10000,'type','Pearson');
[spin_p(2,1),~]=spin_test(HC_G2,-maskG(:,1),'surface_name','fsa5',...
    'parcellation_name','schaefer_400','n_rot',10000,'type','Pearson');
[spin_p(2,2),~]=spin_test(HC_G2,maskG(:,2),'surface_name','fsa5',...
    'parcellation_name','schaefer_400','n_rot',10000,'type','Pearson');

