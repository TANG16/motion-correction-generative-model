% here crossvalidation for 'regsc' setting
% fit the model after leaving one of the echos out and then assess the 
% probability of the left out echo according to the model. This could be 
% repeated for all echos from the various runs.  In principle though, the 
% values could be set to be very small so that it just gives the maximum 
% likelihood estimates.  I suspect that MTV regularisation would do a better 
% job than the current regularisation strategy
% 

clear 
close all
clc

addpath(fullfile(pwd,'utilities'))
addpath(fullfile(pwd,'utilities','auxiliary-functions-UCL'))


%% step 0: initializzation
% please for going back to full TP & full image check following settings:

do_gauss = true; % boolean for choosing the model to use(Gaussian or Rician)
timePoint_step=2; % used in the loop for looping over all TP or using step/increment 'timePoint_step'
makeMovie=1; % boolean for geneerating movies of iteration plots
doSave=1; % for saving all variables at the end of the main loop



%% --input
% input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reduceSize=0; % reduce size of the echoes maps using only lower part of images -> used data with reduced size in the appropiate folder

mainOutputFolderName='3-crossValidation';
outputFolderName_thisCase='1-test'; % to evaluate outputs from different tests/improvements/settings

% from John
% P = spm_select(22,'nifti','Select the raw MRI');
% T    = {[2.3 4.6 6.9 9.2 11.5 13.8 16.1 18.4],[2.3 4.6 6.9 9.2 11.5 13.8 16.1 18.4],[2.3 4.6 6.9 9.2 11.5 13.8]}; 

% here the same from MA and real T from images
[mainPath, folderName]=fileparts(pwd);
dataFolder='0-data_toUse'; % folder with reduced size images created in advance

if reduceSize
    filterStr='^reduced_subj.*';
else
    filterStr='^subj.*';
end

subjDataFolder=fullfile(mainPath,folderName,dataFolder,'PL01_Day00');
[pathData, subjFileName]=fileparts(subjDataFolder);
% output folder
outputFolder=fullfile(mainPath,folderName,mainOutputFolderName,outputFolderName_thisCase,subjFileName);



% read inputs:
% here automatic filled data, in the future from user in a batche format:
% here I read all echoes at once and then group usign path:
thisSubjFiles=spm_select('FPListRec',subjDataFolder,filterStr);
if isempty(thisSubjFiles)
    error('no selected data: please check path to data')
end


%% call the fcn
% some common computation needed for cross-validation:
mainMapPath=cellfun(@(tmp)fileparts(tmp ),cellstr( thisSubjFiles),'UniformOutput', false);
[mapIdx, ~]=grp2idx(mainMapPath);
TP_eachMap=arrayfun(@(tmp)sum(mapIdx==tmp),1:length(mapIdx));
TP_eachMap=TP_eachMap(~TP_eachMap==0);
[T ,fa,tr]=deal(cell(numel(TP_eachMap),1));    

hdr=spm_vol(thisSubjFiles);
% read echoes flip angles for each contrast
ii=0;
% read echo time values from header
for iMap=1:numel(TP_eachMap)
    thisMap_hdr=hdr(mapIdx==iMap);

    for iTP=1:TP_eachMap(iMap)
        ii = ii + 1;
      % from description read echo time info
      tmp = regexp(thisMap_hdr(iTP).descrip,...
         'TR=(?<tr>.+)ms/TE=(?<te>.+)ms/FA=(?<fa>.+)deg',...
         'names');
        T{iMap}(iTP)=str2double(tmp.te);
    end
        fa{iMap}=str2double(tmp.fa);
        tr{iMap}=str2double(tmp.tr);
end

T_all=cat(2,T{:});


% leave 1 out:
regsc=[ 1;  10;  10; 10];
tpm=arrayfun(@(tpm)[num2str(tpm) '_'],regsc,'UniformOutput',0);
regscStr=cat(2,tpm{:});

[MSD, MAD]= deal(zeros(size(thisSubjFiles,1),1));
for iLeaveOut=1:size(thisSubjFiles,1)
    
    this=thisSubjFiles;
    leftOut=this(iLeaveOut,:);
    this(iLeaveOut,:)=[];
    
    thisCross=['cross-validat_leave-1-out_elem_' num2str(iLeaveOut)];
        
    thisOutputFolder=fullfile([outputFolder '_regsc_' regscStr(1:end-1) ],thisCross);
    if ~exist(thisOutputFolder,'dir')
        sts = mkdir(thisOutputFolder);
        if ~sts, error('Error creating output directory "%s".',thisOutputFolder); end
    end 
    
  [ output_args ] = mapGeneration_fcn_vs20190606(this,thisOutputFolder,...
     'do_gauss',do_gauss,'timePoint_step',timePoint_step,'makeMovie',makeMovie,'doSave',doSave,'regsc',regsc);
 
    %  compute the left out from the model:
    mapID= regexp(leftOut,'fil_(\w+).*_siemens','tokens');
    switch char(mapID{1})   
        case {'mt'}
            thisMap='MT';
        case {'pd'}
            thisMap='PD';
        case {'t1'}
            thisMap='T1';                      
    end
    alpha=spm_select('FPListRec',thisOutputFolder,['alpha_' thisMap]);
    beta=spm_select('FPListRec',thisOutputFolder,'beta_R2s.nii');
    
%     beta=spm_vol([thisOutputFolder filesep 'beta_R2s.nii']);
    
    % compute model:
    formula=sprintf('exp(i1-i2.*%.2f)',T_all(iLeaveOut));
    iDataFit=[thisOutputFolder filesep 'fitData.nii'];
    modelThis=spm_imcalc({alpha,beta},iDataFit,formula);
    spm_check_registration(char({alpha;beta;iDataFit;leftOut}) )
    
    
    % comparison with the real acquisition using 
    % (1)mean squared difference:
    differenceMap=[thisOutputFolder filesep 'leftOut_min_modelFit.nii'];
    differenceMapComputed=  spm_imcalc({leftOut,iDataFit},differenceMap,'i1-i2');
    valVectorDiff=differenceMapComputed.private.dat(:);
    MSD_this=sum(valVectorDiff.^2)/prod(differenceMapComputed.dim);
    % (2)mean absolute difference:
    MAD_this=sum(abs(valVectorDiff))/prod(differenceMapComputed.dim);

    MSD(iLeaveOut)=MSD_this;
    MAD(iLeaveOut)=MAD_this;

    % print and save differnce computation:
    sprintf('mean squared difference: %.4f\nmean absolute difference: %.4f',MSD_this,MAD_this)
    save([thisOutputFolder filesep 'crossValidationsValues.mat'],'MSD_this','MAD_this')
end
save([fileparts(thisOutputFolder) filesep 'crossValidationsValues.mat'],'MSD','MAD')


