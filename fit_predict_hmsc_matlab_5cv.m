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
dsz=1;    % data size
typ=3; 
% 1=without LF, 2==with LF; 3=with spatial structure as LF

cd(hmscPath);
Sets = {'plant'};
dSizes=[60,299];
nsets=size(Sets,2);
predN = 100;
set_no=Sets{1};      % data set
dSz=dSizes(dsz);    % data size

folder = fullfile(wdpath,'FITS',set_no,'HMSC');
folderData = fullfile(wdpath,'DATA');
folderPred = fullfile(wdpath,'PREDICTIONS',set_no,'HMSC');

for dTyp=1:5
    
    if dTyp==1
        a=1; b=2; c=3; d=4; e=5;
    end
    if dTyp==2
        a=2; b=3; c=4; d=5; e=1;
    end
    if dTyp==3
        a=3; b=4; c=5; d=1; e=2;
    end
    if dTyp==4
        a=4; b=5; c=1; d=2; e=3;
    end
    if dTyp==5
        a=5; b=1; c=2; d=3; e=4;
    end
    
    file=fullfile(folderData,strcat('Y_', num2str(a), '_', num2str(set_no),'.csv'));
    Y_1=importdata(file);
    file=fullfile(folderData,strcat('Y_', num2str(b), '_', num2str(set_no),'.csv'));
    Y_2=importdata(file);
    file=fullfile(folderData,strcat('Y_', num2str(c), '_', num2str(set_no),'.csv'));
    Y_3=importdata(file);
    file=fullfile(folderData,strcat('Y_', num2str(d), '_', num2str(set_no),'.csv'));
    Y_4=importdata(file);
    file=fullfile(folderData,strcat('Y_', num2str(e), '_', num2str(set_no),'.csv'));
    Y_5=importdata(file);
    %train
    Y_t1=[Y_1; Y_2; Y_3; Y_4];
    Y_t=[ones(size(Y_t1,1),1),Y_t1];
    %valid
    Y_v1=Y_5;
    Y_v=[ones(size(Y_v1,1),1),Y_v1];

    file=fullfile(folderData,strcat('X_', num2str(a), '_', num2str(set_no),'.csv'));
    X_1=importdata(file);
    file=fullfile(folderData,strcat('X_', num2str(b), '_', num2str(set_no),'.csv'));
    X_2=importdata(file);
    file=fullfile(folderData,strcat('X_', num2str(c), '_', num2str(set_no),'.csv'));
    X_3=importdata(file);
    file=fullfile(folderData,strcat('X_', num2str(d), '_', num2str(set_no),'.csv'));
    X_4=importdata(file);
    file=fullfile(folderData,strcat('X_', num2str(e), '_', num2str(set_no),'.csv'));
    X_5=importdata(file);
    %train
    X_t1=[X_1; X_2; X_3; X_4];
    X_t=X_t1;
    %X_t=[ones(size(X_t1,1),1),X_t1];
    %valid
    X_v1=X_5;
    %X_v=[ones(size(X_v1,1),1),X_v1];
    X_v=X_v1;

    file=fullfile(folderData,strcat('S_', num2str(a), '_', num2str(set_no),'.csv'));
    S_1=importdata(file);
    file=fullfile(folderData,strcat('S_', num2str(b), '_', num2str(set_no),'.csv'));
    S_2=importdata(file);
    file=fullfile(folderData,strcat('S_', num2str(c), '_', num2str(set_no),'.csv'));
    S_3=importdata(file);
    file=fullfile(folderData,strcat('S_', num2str(d), '_', num2str(set_no),'.csv'));
    S_4=importdata(file);
    file=fullfile(folderData,strcat('S_', num2str(e), '_', num2str(set_no),'.csv'));
    S_5=importdata(file);
    %train
    S_t1=[S_1; S_2; S_3; S_4];
    S_t=[ones(size(S_t1,1),1),S_t1];
    %valid
    S_v1=S_5;
    S_v=[ones(size(S_v1,1),1),S_v1];

    Xs_t=[X_t,S_t];
    Xs_v=[X_v,S_v];

    nsp=size(Y_t,2);
    nsites=size(Y_t,1);

    piCell = cellfun(@num2str, num2cell((1:size(Y_t,1))'), 'UniformOutput', false);
    xyCell = [piCell, num2cell(S_t)];
    m = Hmsc(folder, false, false, [true]);
    m.setData(Y_t,'probit',X_t,piCell,{xyCell},[],[]);
    m.setPriorsDefault();
    covScaleFlag=ones(1,m.nc);
    covScaleFlag(1)=2;
    m.setCovScaling(covScaleFlag);

    m.setMCMCOptions(niters, thin1);
    m.setMCMCAdapt([nadaptrounds,0], true);
    m.setMCMCSaveOptions(MCMCsaveRam, MCMCsaveFile);

    m.sampleMCMC(nrounds, false, [], 3);

    m.setPostThinning(MCMCcut:m.repN, thin2);

    filebody=strcat('hmsc_',num2str(set_no),'_',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));

    try
        save(fullfile(folder,strcat(filebody,'.mat')),'m');
    catch
        fprintf('File too large, using flag -v7.3')
        save(fullfile(folder,strcat(filebody,'.mat')),'m','-v7.3');
    end

    m.postRamClear();

    load(fullfile(folder,strcat(filebody,'.mat')));

    piCell = cellfun(@num2str, num2cell((1:(m.ny+size(Y_v,1)))'), 'UniformOutput', false);
    xyCell = [piCell, num2cell([S_t;S_v])];
    predList = m.predict(predN, [m.X;X_v], piCell, {xyCell}, false);

    filebodyPred=strcat('pred_hmsc_',num2str(set_no),'_',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));
    filebodyPredCsv=strcat('preds_',num2str(set_no),'_hmsc',num2str(typ),'_d',num2str(dTyp),'_',num2str(dSz));

    save(fullfile(folderPred,strcat(filebodyPred,'.mat')),'predList');
    predsM=cell2mat(predList);
    csvwrite(fullfile(folderPred,strcat(filebodyPredCsv,'.csv')),predsM);
end