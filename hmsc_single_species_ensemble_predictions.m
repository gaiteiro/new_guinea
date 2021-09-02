
% Fit single species models for data with explicit spatial structure

% code adapted from: %
% Norberg, A., Abrego, N., Blanchet, F. G., Adler, F. R., Anderson, B. J., Anttila, J., Araújo, M. B.,
% Dallas, T., Dunson, D., Elith, J., Foster, S. D., Fox, R., Franklin, J., Godsoe, W., Guisan, A., 
% O'Hara, B., Hill, N. A., Holt, R. D., Hui, F. K. C., . Ovaskainen, O. (2019). 
% A comprehensive evaluation of predictive performance of 33 species distribution models at species and 
% community levels. Ecological Monographs, 89, e01370. https://doi.org/10.1002/ecm.1370

addpath('C:/Users/smith/Documents/Norberg comparison/aminorberg-SDM-comparison-a8028bf/bakeoff/pipeline/MODELS/HMSC class');
wdpath="C:/Users/smith/Documents/Norberg comparison/aminorberg-SDM-comparison-a8028bf/bakeoff/pipeline/";
cd("C:/Users/smith/Documents/Norberg comparison/aminorberg-SDM-comparison-a8028bf/bakeoff/pipeline");

hmscPath=fullfile(wdpath,'MODELS','HMSC class');
niters=100;
thin1=10;
nadaptrounds=3;
nrounds=9;
thin2=30;
MCMCcut=6;
MCMCsaveRam=true;
MCMCsaveFile=false;
intXs=false;
commSP=false;
MCMC2=false;
dsz=1;    % data size
typ=2; % spatial
     
cd(hmscPath);
Sets = {'plant'};
dSizes=[146,293];
nsets=size(Sets,2);
predN = 100;
set_no=Sets{1};      % data set
dSz=dSizes(dsz);    % data size
cd(hmscPath);

dSz=dSizes(dsz);
set_no=Sets{1};

mtypes = 1:2;
dTyp=1;
folder = fullfile(wdpath,'FITS',set_no,'ssHMSC');
folderData = fullfile(wdpath,'DATA');
folderPred = fullfile(wdpath,'PREDICTIONS',set_no,'ssHMSC');

dTyp=1;
%train
file=fullfile(folderData,strcat('Yt_', num2str(dTyp), '_', num2str(set_no),'.csv'));
Y_t=importdata(file);
%valid
file=fullfile(folderData,strcat('Yv_', num2str(dTyp), '_', num2str(set_no),'.csv'));
Y_v=importdata(file);

%train
file=fullfile(folderData,strcat('Xt_', num2str(dTyp), '_', num2str(set_no),'.csv'));
X_t=importdata(file);
X_t=[ones(size(X_t,1),1),X_t,X_t.^2];
%valid
file=fullfile(folderData,strcat('Xv_', num2str(dTyp), '_', num2str(set_no),'.csv'));
X_v=importdata(file);
X_v=[ones(size(X_v,1),1),X_v,X_v.^2];

%train
file=fullfile(folderData,strcat('St_', num2str(dTyp), '_', num2str(set_no),'.csv'));
S_t=importdata(file);
%valid
file=fullfile(folderData,strcat('Sv_', num2str(dTyp), '_', num2str(set_no),'.csv'));
S_v=importdata(file);

Xs_t=[X_t,S_t];
Xs_v=[X_v,S_v];

nsp=size(Y_t,2);
nsites=size(Y_t,1);

compTime=0;

for sp=1:nsp

    piCell = cellfun(@num2str, num2cell((1:size(Y_t,1))'), 'UniformOutput', false); % cell array numbered 1 to 146
    xyCell = [piCell, num2cell(S_t)]; % col 1 = ID, cols 2 & 3 are coords
    m = Hmsc(folder, false, false, [true]);
    m.setData(Y_t(:,sp),'probit',X_t,piCell,{xyCell},[],[]);
  
    m.setPriorsDefault();
    covScaleFlag=ones(1,m.nc);
    covScaleFlag(1)=2;
    m.setCovScaling(covScaleFlag);

    m.setMCMCOptions(niters, thin1);
    m.setMCMCAdapt([nadaptrounds,0], true);
    m.setMCMCSaveOptions(MCMCsaveRam, MCMCsaveFile);

    filebody=strcat('sp',num2str(sp),'_',num2str(set_no),'_hmsc',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    filebodyPred=strcat('sp',num2str(sp),'_',num2str(set_no),'_pred_hmsc',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    filebodyPredCsv=strcat('preds_ss_',num2str(set_no),'_hmsc',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    filebodyTime=strcat('compTime_',num2str(set_no),'_hmsc',num2str(typ),'_',num2str(dSz));

    try
        tic;
        m.sampleMCMC(nrounds, false, [], 3);
        compT=toc;
    catch
        fprintf(strcat('MCMC sampling failed for sp', num2str(sp)));
        save(fullfile(folder,strcat('sp',num2str(sp),'model',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz),'_NULL','.mat')),'m');
    continue;
end


m.setPostThinning(MCMCcut:m.repN, thin2);

save(fullfile(folder,strcat(filebody,'.mat')),'m');

m.postRamClear();
load(fullfile(folder,strcat(filebody,'.mat')));

if typ==1
    predList = m.predict(predN, X_v, [], [], false);
end
if typ==2
    piCell = cellfun(@num2str, num2cell((1:(m.ny+size(Y_v,1)))'), 'UniformOutput', false);
    xyCell = [piCell, num2cell([S_t;S_v])];
    predList = m.predict(predN, [m.X;X_v], piCell, {xyCell}, false);
end

save(fullfile(folderPred,strcat(filebodyPred,'.mat')),'predList');
compTime=compTime+compT;
end


if dTyp==1
save(fullfile(folder,strcat(filebodyTime,'.mat')),'compTime');
end
        
        % subsamples of data and species
        %file=fullfile(folderData,strcat('siteSamps_',num2str(set_no),'.mat'));
        %siteSamps=importdata(file);
        %siteSamps=struct2cell(siteSamps);
        %samp=siteSamps{dsz}';
        %file=fullfile(folderData,strcat('spSel_',num2str(set_no),'.csv'));
        %spSel=importdata(file);

        %Y_t=Y_t(samp,spSel);
        %Y_v=Y_v(:,spSel);
        %X_t=X_t(samp,:);
        %S_t=S_t(samp,:);
        %Xs_t=Xs_t(samp,:);


% modify predictions to .csv
cd(hmscPath);
dSz=dSizes(dsz);
set_no=Sets{1};

file=fullfile(folderData,strcat('Yv_', num2str(dTyp), '_', num2str(set_no),'.csv'));
Y_v=importdata(file);

% subsample of species
%file=fullfile(folderData,strcat('spSel_',num2str(set_no),'.csv'));
%spSel=importdata(file);

%Y_v=Y_v(:,spSel);

nsp=size(Y_v,2);
predsM=[];

filebodyPredCsv=strcat('preds_ss_',num2str(set_no),'_hmsc',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));


for sp=1:nsp
    filebodyPred=strcat('sp',num2str(sp),'_',num2str(set_no),'_pred_hmsc',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    load(fullfile(folderPred,strcat(filebodyPred,'.mat')),'predList');
    predsM=[predsM,cell2mat(predList)];
end

csvwrite(fullfile(folderPred, strcat(filebodyPredCsv,'.csv')),predsM);

% in reverse

%train
file=fullfile(folderData,strcat('Yv_', num2str(dTyp), '_', num2str(set_no),'.csv'));
Y_t=importdata(file);
%valid
file=fullfile(folderData,strcat('Yt_', num2str(dTyp), '_', num2str(set_no),'.csv'));
Y_v=importdata(file);

%train
file=fullfile(folderData,strcat('Xv_', num2str(dTyp), '_', num2str(set_no),'.csv'));
X_t=importdata(file);
X_t=[ones(size(X_t,1),1),X_t,X_t.^2];
%valid
file=fullfile(folderData,strcat('Xt_', num2str(dTyp), '_', num2str(set_no),'.csv'));
X_v=importdata(file);
X_v=[ones(size(X_v,1),1),X_v,X_v.^2];

%train
file=fullfile(folderData,strcat('Sv_', num2str(dTyp), '_', num2str(set_no),'.csv'));
S_t=importdata(file);
%valid
file=fullfile(folderData,strcat('St_', num2str(dTyp), '_', num2str(set_no),'.csv'));
S_v=importdata(file);

Xs_t=[X_t,S_t];
Xs_v=[X_v,S_v];

nsp=size(Y_t,2);
nsites=size(Y_t,1);

compTime=0;


for sp=1:nsp

    piCell = cellfun(@num2str, num2cell((1:size(Y_t,1))'), 'UniformOutput', false); % cell array numbered 1 to 146
    xyCell = [piCell, num2cell(S_t)]; % col 1 = ID, cols 2 & 3 are coords
    m = Hmsc(folder, false, false, [true]);
    m.setData(Y_t(:,sp),'probit',X_t,piCell,{xyCell},[],[]);
  
    m.setPriorsDefault();
    covScaleFlag=ones(1,m.nc);
    covScaleFlag(1)=2;
    m.setCovScaling(covScaleFlag);

    m.setMCMCOptions(niters, thin1);
    m.setMCMCAdapt([nadaptrounds,0], true);
    m.setMCMCSaveOptions(MCMCsaveRam, MCMCsaveFile);

    filebody=strcat('sp',num2str(sp),'_',num2str(set_no),'_hmscb',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    filebodyPred=strcat('sp',num2str(sp),'_',num2str(set_no),'_pred_hmscb',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    filebodyPredCsv=strcat('preds_ss_',num2str(set_no),'_hmscb',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    filebodyTime=strcat('compTime_',num2str(set_no),'_hmscb',num2str(typ),'_',num2str(dSz));

    try
        tic;
        m.sampleMCMC(nrounds, false, [], 3);
        compT=toc;
    catch
        fprintf(strcat('MCMC sampling failed for sp', num2str(sp)));
        save(fullfile(folder,strcat('sp',num2str(sp),'modelb',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz),'_NULL','.mat')),'m');
    continue;
end


m.setPostThinning(MCMCcut:m.repN, thin2);

save(fullfile(folder,strcat(filebody,'.mat')),'m');

m.postRamClear();
load(fullfile(folder,strcat(filebody,'.mat')));

if typ==1
    predList = m.predict(predN, X_v, [], [], false);
end
if typ==2
    piCell = cellfun(@num2str, num2cell((1:(m.ny+size(Y_v,1)))'), 'UniformOutput', false);
    xyCell = [piCell, num2cell([S_t;S_v])];
    predList = m.predict(predN, [m.X;X_v], piCell, {xyCell}, false);
end

save(fullfile(folderPred,strcat(filebodyPred,'.mat')),'predList');
compTime=compTime+compT;
end


if dTyp==1
save(fullfile(folder,strcat(filebodyTime,'.mat')),'compTime');
end
        
        % subsamples of data and species
        %file=fullfile(folderData,strcat('siteSamps_',num2str(set_no),'.mat'));
        %siteSamps=importdata(file);
        %siteSamps=struct2cell(siteSamps);
        %samp=siteSamps{dsz}';
        %file=fullfile(folderData,strcat('spSel_',num2str(set_no),'.csv'));
        %spSel=importdata(file);

        %Y_t=Y_t(samp,spSel);
        %Y_v=Y_v(:,spSel);
        %X_t=X_t(samp,:);
        %S_t=S_t(samp,:);
        %Xs_t=Xs_t(samp,:);


% modify predictions to .csv
cd(hmscPath);
dSz=dSizes(dsz);
set_no=Sets{1};

file=fullfile(folderData,strcat('Yt_', num2str(dTyp), '_', num2str(set_no),'.csv'));
Y_v=importdata(file);

% subsample of species
%file=fullfile(folderData,strcat('spSel_',num2str(set_no),'.csv'));
%spSel=importdata(file);

%Y_v=Y_v(:,spSel);

nsp=size(Y_v,2);
predsM=[];

filebodyPredCsv=strcat('preds_ss_',num2str(set_no),'_hmscb',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));


for sp=1:nsp
    filebodyPred=strcat('sp',num2str(sp),'_',num2str(set_no),'_pred_hmscb',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    load(fullfile(folderPred,strcat(filebodyPred,'.mat')),'predList');
    predsM=[predsM,cell2mat(predList)];
end

csvwrite(fullfile(folderPred, strcat(filebodyPredCsv,'.csv')),predsM);