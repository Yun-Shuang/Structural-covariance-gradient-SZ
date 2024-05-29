
clear;clc;
%% Figure3: Projecting binned gradients axis onto covariance matrix
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

%% Gradient BINs
load('/media/shuang/data/Covariance/results/Gradients.mat'); % own 400*10*2(EOS;HC) gradients
load('/media/shuang/data/Covariance/results/conn_HC_cov.mat'); % covariance matrix
load('/media/shuang/data/Covariance/results/conn_EOS_cov.mat');

BIN=10; % divide into 10 bins
[~,index1]=sort(Gradient(:,2,2)); % HC small to big index; G2: anterior-to-posterior
[~,rank1]=sort(index1); % position for each element
BinsHC=reshape(index1,400/BIN,BIN); % index for 10 bins
maskHC=zeros(400,1);
for i = 1:BIN
    maskHC(BinsHC(:,i))=i;
    CTbinsHC(i,:)=mean(CT400(BinsHC(:,i),hc));
    CTbinsEOS(i,:)=mean(CT400(BinsHC(:,i),eos));
end

%-----------------plot binned maps---------------
obj=plot_hemispheres(maskHC,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'HC_10bins'});
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([10],"BrBg"))])
obj.colorlimits([0.9, 10]) 

test=partialcorr(CTbinsHC',regressors(hc,:));  % divide covariance matrix to 10 bins
figure;heatmap(test); caxis([-1 1]); colormap(flipud(brewermap([],"RdBu")));
test=partialcorr(CTbinsEOS',regressors(eos,:)); 
figure;heatmap(test); caxis([-1 1]); colormap(flipud(brewermap([],"RdBu")));
%-----------------------end-----------------------

%% Bins-based general linear model
%----------------set variables and covariates-----------
age=regressors(:,1);
gender = regressors(:,2);
for i =1:194
    if gender(i)==1
        gender1{i,1}= 'M';
    elseif gender(i)==2
        gender1{i,1}='F';
    end
end
for i=1:95
    diagnose{i,1}='S';
end
for i=96:length(age)
    diagnose{i,1}='H';
end
for i = 1:BIN
    seed(i,:)=mean(CT400(BinsHC(:,i),:));
end
demographics = table(age,diagnose,gender1,regressors(:,3),...
    seed(1,:)',seed(2,:)',seed(3,:)',seed(4,:)',seed(5,:)',...
    seed(6,:)',seed(7,:)',seed(8,:)',seed(9,:)',seed(10,:)',...
    'VariableNames',{'age','diag','sex','meanCT',...
    'bin1','bin2','bin3','bin4','bin5','bin6','bin7','bin8','bin9','bin10',});
term_meanCT = FixedEffect(demographics.meanCT, 'meanCT');
term_sex = FixedEffect(demographics.sex);
term_age = FixedEffect(demographics.age, 'Age');
term_diag = FixedEffect(demographics.diag);
%-------------------------end--------------------------

% Diagnose*seed effect
for k=1:10
    clear tempo Seed term_seed
    tempo=demographics(:,4+k);
    Seed = table2array(tempo);
    term_seed = FixedEffect(Seed, 'seed');
    model_diag = term_meanCT + term_age + term_sex + ...
     term_seed + term_diag + term_seed * term_diag;
    contrast_symp = (Seed .* ...
                    (demographics.diag == "S")) - ...
                   (Seed .* ...
                   (demographics.diag == "H"));
    slm_diag = SLM( ...
    model_diag, ...
    contrast_symp, ...
    'correction',{'fdr'});
    slm_diag.fit(seed');
    tt(k,:)=slm_diag.t;
    p(k,:)=tcdf(-abs(slm_diag.t),slm_diag.df);
    FDRvalue(k,:) = mafdr(p(k,:),'BHFDR',true);
end
p(find(p>0.05))=0;
FDRvalue(find(FDRvalue>0.05))=0;
tt = tt - diag(diag(tt)); % remove diagonal elements
figure;heatmap(tt.*logical(FDRvalue));  caxis([-2,2]); colormap(flipud(brewermap([],"RdBu")));

%% plot FDR results: bin1-9; bin7-10
source=1;target=9;
source=7;target=10;

% plot source and target brain regions
maskHC1=maskHC;
maskHC1(find(maskHC1~=source))=0;
obj=plot_hemispheres(maskHC1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'HC_source'});
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([0.9, 10])
maskHC1=maskHC;
maskHC1(find(maskHC1~=target))=0;
obj=plot_hemispheres(maskHC1,{surf_lh,surf_rh}, ...
             'parcellation', labeling, ...
             'labeltext',{'HC_target'});
obj.colormaps([0.7 0.7 0.7; flipud(brewermap([],"BrBg"))])
obj.colorlimits([0.9, 10]) 

% plot correlations
CTbins1=regress_cov([CTbinsEOS';CTbinsHC'],regressors)';
CTbinsEOS1=CTbins1(:,1:95);CTbinsHC1=CTbins1(:,96:end);
figure; plot(CTbinsHC1(source,:),CTbinsHC1(target,:),'.','MarkerEdge',[77,77,77]/255,'markersize',25);
hold on;plot(CTbinsHC1(source,:),polyval(polyfit(CTbinsHC1(source,:),CTbinsHC1(target,:),1),CTbinsHC1(source,:)),'k')
hold on;plot(CTbinsEOS1(source,:),CTbinsEOS1(target,:),'.','MarkerEdge',[51,160,44]/255,'markersize',25);
hold on;plot(CTbinsEOS1(source,:),polyval(polyfit(CTbinsEOS1(source,:),CTbinsEOS1(target,:),1),CTbinsEOS1(source,:)))

%% Association with symptoms
%%% load panss
load('/media/shuang/data/Covariance/data/PANSS.mat');
clc;fprintf('PANSS-P = %f +- %f, PANSS-N = %f +- %f ',...
    mean(cell2mat(celltypes(:,6))), std(cell2mat(celltypes(:,6))),...
    mean(cell2mat(celltypes(:,7))), std(cell2mat(celltypes(:,7))));
fprintf('PANSS-G = %f +- %f, PANSS-T = %f +- %f ',...
    mean(cell2mat(celltypes(:,8))), std(cell2mat(celltypes(:,8))),...
    mean(cell2mat(celltypes(:,10))), std(cell2mat(celltypes(:,10))));
list = '/media/shuang/data/p_02658/Lists/localSCN9599.txt'; % delete SCH108 for missing of rh131 parcel
fileID = fopen(list);
sbj = textscan(fileID,'%s'); sbj = sbj{1};
fclose(fileID);
index = zeros(size(sbj));
for i = 1:length(index)
    index(i) = ~isempty(find(strcmp(celltypes,sbj{i})));
end
sum(index)
index=logical(index);
demographics=demographics(index,:);
demographics=[demographics,cell2table(celltypes(:,6:10),"VariableNames",["panssp","panssn","panssg","pansss","pansst"]);];
term_meanCT = FixedEffect(demographics.meanCT, 'meanCT');
term_sex = FixedEffect(demographics.sex);
term_age = FixedEffect(demographics.age, 'Age');
term_diag = FixedEffect(demographics.diag);

symp = demographics.panssp; % select a symptom domain
term_panss = FixedEffect(symp);
%-------------------------end--------------------------

% Symptom*seed effect
for k=1:10
    clear tempo Seed term_seed
    tempo=demographics(:,4+k);
    Seed = table2array(tempo);
    term_seed = FixedEffect(Seed, 'seed');
    model_panss = term_meanCT + term_age + term_sex + ...
     term_seed + term_panss + term_seed * term_panss;
    contrast_symp = Seed .*symp;
    slm_diag = SLM( ...
    model_panss, ...
    contrast_symp, ...
    'correction',{'fdr'});
    slm_diag.fit(seed(:,index)');
    tt(k,:)=slm_diag.t;
    p(k,:)=tcdf(-abs(slm_diag.t),slm_diag.df); 
    FDRvalue(k,:) = mafdr(p(k,:),'BHFDR',true);
end
p(find(p>0.05))=0;
FDRvalue(find(FDRvalue>0.05))=0;
tt = tt - diag(diag(tt)); % remove diagonal elements
figure;heatmap(tt.*logical(FDRvalue));  caxis([-2,2]); colormap(flipud(brewermap([],"RdBu")));
