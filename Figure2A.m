
clear;clc
%% Figure2: structural covariance gradients comparisons
addpath(genpath('/media/shuang/data/Covariance/packages/cifti-matlab-master'));

%% load gradient data
load('/media/shuang/data/Covariance/results/Gradients.mat');
HC_G1=Gradient(:,1,2);EOS_G1=Gradient(:,1,1);
HC_G2=Gradient(:,2,2);EOS_G2=Gradient(:,2,1); %Gradient 2

%% ks-stat for gradient distribution
[h,p,ksstat]=kstest2(HC_G1,EOS_G1)
[h,p,ksstat]=kstest2(HC_G2,EOS_G2)

%% z-test for gradient values
HC_G1=zscore(HC_G1);EOS_G1=zscore(EOS_G1); % normalized
z=zeros(400,1);zp=zeros(400,1);
for i = 1:400
    clear t_r1 t_r2
    t_r1=EOS_G1(i);t_r2=HC_G1(i);
    z(i) = t_r1-t_r2;
    zp(i) = (1-normcdf(abs(z(i)),0,1))*2;
end
fdrz = mafdr(zp,'BHFDR',true);

%%%%% save as dlabel.nii mask file
path_wb_command = ['/media/shuang/data/Covariance/packages/workbench/bin_linux64/wb_command'];
mask=ciftiopen('/media/shuang/data/Covariance/template/Schaefer2018_400Parcels_7Networks_order.dlabel.nii',path_wb_command);
for i = 1:400 % fs_LR 32k space for saving dlabel.nii file
    mask.cdata(find(mask.cdata==i))=logical(zp(i)<(0.05));
end
ciftisavereset(mask,'/media/shuang/data/Covariance/results/g2z05.dlabel.nii',path_wb_command);

mask=ciftiopen('/media/shuang/data/Covariance/template/SCON003_400.pscalar.nii',path_wb_command);
mask.cdata=z;
ciftisavereset(mask,'/media/shuang/data/Covariance/results/g2z.pscalar.nii',path_wb_command);

%% 7 YEO functional networks
load('/media/shuang/data/Covariance/data/N7index.mat')
comp=1;
clear g7 test p7 t7
for i=1:7
    g7{i,1}=Gradient(find(N7==i),comp,1); % EOS patients
    test{i,1}=(Gradient(find(N7==i),comp,1));
    g7{i,2}=Gradient(find(N7==i),comp,2); % controls
    test{i,2}=(Gradient(find(N7==i),comp,2));
    [~,p7(i),~,t7(i)]=ttest(test{i,1},test{i,2}); % paired-t tests for each networks
end
fdr1 = mafdr(p7,'BHFDR',true);
