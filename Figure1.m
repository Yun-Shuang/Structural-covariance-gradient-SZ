
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

%% show an example of structural covariance
figure; plot(CT400(149,hc),CT400(91,hc),'k.','markersize',25);ylim([1.4 3.4]);xlim([2.2 3.8]);
hold on; plot(CT400(149,hc),polyval(polyfit(CT400(149,hc),CT400(91,hc),1),CT400(149,hc)),'r')
figure; plot(CT400(149,eos),CT400(91,eos),'k.','markersize',25);ylim([1.4 3.4]);xlim([2.2 3.8]);
hold on; plot(CT400(149,eos),polyval(polyfit(CT400(149,eos),CT400(91,eos),1),CT400(149,eos)),'r')

%% structural covariance matrix
[conn_HC_cov,conn_p(:,:,1)] = partialcorr(HC',regressors(hc,:)); % 
[conn_EOS_cov,conn_p(:,:,2)] = partialcorr(EOS',regressors(eos,:));
% ---------------- plot no/threshold covariance matrix ---------------
figure; heatmap(conn_HC_cov,'ColorLimits',[-0.5 0.5],'GridVisible','off'); colormap(flipud(brewermap([],"RdBu")))
figure; heatmap(conn_EOS_cov,'ColorLimits',[-0.5 0.5],'GridVisible','off'); colormap(flipud(brewermap([],"RdBu")))
zz=conn_HC_cov;maskz=zeros(400,400);
[~,index] = sort(zz,'descend');
for i=1:400
    maskz(index(1:40,i),i)=1;
end
figure; heatmap(maskz,'ColorLimits',[0 1],'GridVisible','off'); colormap(brewermap([],"Greys"))
zz=conn_EOS_cov;maskz=zeros(400,400);
[~,index] = sort(zz,'descend');
for i=1:400
    maskz(index(1:40,i),i)=1;
end
figure; heatmap(maskz,'ColorLimits',[0 1],'GridVisible','off'); colormap(brewermap([],"Greys"))
% ------------------ end ----------------------------------

%% Diffusion gradients
maskG=load('/media/shuang/data/Covariance/template/strcov_gradient.csv'); %HCP 400*10 gradient mask
% ------------------- gradient calculating -----------------
Galign = GradientMaps('approach','diffusion embedding','alignment','pa');
Galign_HC = Galign.fit(conn_HC_cov,'reference',maskG);
Galign_EOS = Galign.fit(conn_EOS_cov,'reference',maskG);
Gradient(:,:,1)=Galign_EOS.aligned{1};Gradient(:,:,2)=Galign_HC.aligned{1};
% --------------------------- end -----------------------------

scree_plot(Galign_HC.lambda{1}); % explained variance
scree_plot(Galign_EOS.lambda{1}); 

HC_G1=-Gradient(:,1,2);EOS_G1=-Gradient(:,1,1);
% print mean +- SD
fprintf('HC gradient = %f +- %f, EOS gradient = %f +- %f ', mean(HC_G1), std(HC_G1), mean(EOS_G1), std(EOS_G1));
% --------------------Gradient1 ------------------------------
obj=plot_hemispheres(-maskG(:,1),{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','HCP');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([-0.16, 0.16]) 

obj=plot_hemispheres(HC_G1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','HC');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([-0.05, 0.05]) 

obj=plot_hemispheres(EOS_G1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext','EOS');
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([-0.05, 0.05]) 
% --------------- end ---------------------------

%% 7 YEO functional networks
load('/media/shuang/data/Covariance/data/N7index.mat')
clear g7 test p7 t7
for i=1:7
    g7(i,1)=-mean(Gradient(find(N7==i),1,1)); % EOS patients
    test{i,1}=-(Gradient(find(N7==i),1,1));
    g7(i,2)=-mean(Gradient(find(N7==i),1,2)); % controls
    test{i,2}=-(Gradient(find(N7==i),1,2));
    [~,p7(i),~,t7(i)]=ttest2(test{i,1},test{i,2}); % 2-t tests for each networks
end
% addpath(genpath('/media/shuang/data/p_02658/software/spider_plot-master'));
% figure;
% spider_plot(squeeze(g7)',...
%     'AxesLabels', {'VIS', 'SMN', 'DAN', 'VAN', 'LMB', 'FPN', 'DMN'},...
%     'AxesLimits', [-0.02*ones(1,7); 0.02*ones(1,7)],... % [min axes limits; max axes limits]
%     'AxesPrecision', 2*ones(1,7));
% legend('G1EOS', 'G1HC', 'Location', 'southoutside');

%% 5 cytoarchitectural networks
economo=load('/media/shuang/data/Covariance/template/economo_koskinas_conte69.csv');
mask=economo;mask(~logical(labeling))=0;
for i=1:400
    [mask1(i),f(i)]=mode(economo(find(labeling==i)));
    f(i)=f(i)/length(economo(find(labeling==i)));
end
dis=logical(f<0.5);dis=dis|~logical(mask1);
mask1(dis)=0; % discard 15 parcels
mask1=mask1';
%%% save as dlabel.nii mask file
% mask=ciftiopen('/media/shuang/data/p_02658/AAAAAscripts/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
% for i = 1:400 % fs_LR 32k space for saving dlabel.nii file
%     mask.cdata(find(mask.cdata==i))=mask1(i);
% end
% ciftisavereset(mask,'/media/shuang/data/Covariance/results/ECO5.dlabel.nii',path_wb_command);
clear e5 test p5 t5
for i=1:5
    e5(i,1)=-mean(Gradient(find(mask1==i),1,1)); % EOS patients
    test{i,1}=-(Gradient(find(mask1==i),1,1));
    e5(i,2)=-mean(Gradient(find(mask1==i),1,2)); % controls
    test{i,2}=-(Gradient(find(mask1==i),1,2));
    [~,p5(i),~,t5(i)]=ttest2(test{i,1},test{i,2}); % 2-t tests for each networks
end
C = [124,40,123;48,92,135;160,200,134;187,156,16;225,218,25]./255;
% for i = 1:5
%     test=mask1;test(find(test~=i))=0;
%     obj=plot_hemispheres(test,{surf_lh,surf_rh}, ...
%              'parcellation', labeling);
%     obj.colormaps([0.7 0.7 0.7; C(i,:)])
% end
obj=plot_hemispheres(mask1,{surf_lh,surf_rh},'parcellation', labeling);
obj.colormaps([0.7 0.7 0.7; C]) % cytoarchitectural networks
% figure; % spider plot
% spider_plot(squeeze(e5)',...
%     'AxesLabels', {'Agranular', 'Frontal', 'Parietal', 'Polar', 'Granular'},...
%     'AxesLimits', [-0.03*ones(1,5); 0.03*ones(1,5)],... % [min axes limits; max axes limits]
%     'AxesPrecision', 2*ones(1,5));
% legend('G1EOS', 'G1HC', 'Location', 'southoutside');