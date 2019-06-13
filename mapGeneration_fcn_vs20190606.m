function [ output_args ] = mapGeneration_fcn_vs20190606 (thisSubjFiles,outputFolder,varargin)
% the main function for map generation for MPM protocol
%
% INPUTS: 
%	thisSubjFiles:all echoes all contrast files for a subject
%   outputFolder: Folder where to save outputs
% 
% 
% OPTIONAL INPUTS: 
% optional inptuts seed to be called using string for each parameter i.e.
% [Ell ,phi, g_thisMap, h_thisMap]=mapGeneration_fcn_vs20190606(..,'paramName',paramValue);
% 
%   timePoint_step: increase step for each loop over Nii data
%       the loop goues 1:timePoint_step:numel(Nii)
%   do_gauss : for gaussian (do_gauss=1) or Rice (do_gauss=0) distribution
%   makeMovie: for saving plots at each iterations
%   doSave: for saving all variables at the end of the main loop
%   axesHdl: to plot in specific axes 
%   plotModel: for plotting the loop images note if axesHdl is given this is
%       set to 1 (default is 0, no plot is done)
%
% OUTPUTS:
%   
%   
% MA: 06-06-2019


    %% --input/output 
    % define contrast by grouping based on the folder path:
    mainMapPath=cellfun(@(tmp)fileparts(tmp ),cellstr( thisSubjFiles),'UniformOutput', false);
    [mapIdx, ~]=grp2idx(mainMapPath);
    TP_eachMap=arrayfun(@(tmp)sum(mapIdx==tmp),1:length(mapIdx));
    TP_eachMap=TP_eachMap(~TP_eachMap==0);
    [T ,fa,tr]=deal(cell(numel(TP_eachMap),1));    


    % default inputs:
    def.timePoint_step = 1;
    def.plotModel = 0;
    def.figureHdl=[]; % by default all kind of hv will be considered 
    def.axesHdl=[]; % by default all kind of hv will be considered 
    def.do_gauss=1;
    def.makeMovie=0;
    def.doSave=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Settings for spm_field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters are voxel sizes (x, y and z) followed by 3 regularisation settings
    % and then the number of iterations and multi-grid cycles.
    % If no smoothing is required, then this bit of code could be speeded up lots.
    % Some form of total variation regularisation could be used, in which case
    % the smoothing part would be useful.
    % reg      = [1 1 1  [0 1 0]];
    % % note cahnge the full regsc-> look if the smooth is too much, adding
    % % artifacts
    % regsc    = ones(numel(T)+1,1)*1000;
    % regsc(1) = 100000;

    % default optimization parameters:
    def.reg= [1 1 1  [0 1 0]];
    def.regsc    = ones(numel(T)+1,1)*1000;
    def.regsc(1) = 100000;
    

    % -- settings deformation 
    def.doDef= 0; % boolen to choose if the deformation part should be done or not
    def.vx          = [1 1 1]; % Voxel sizes
    def.settings_v  = [def.vx  0.00001 0.01 10 0.05 0.1]; % Regularisation for displacements
    def.alpha       = 0.9; % Update step (can be reduced if updates overshoot)

    
    % Parameters Validation
    ip = inputParser;
    ip.FunctionName = mfilename;
    for fn = fieldnames(def)'
        ip.addParameter(fn{1},def.(fn{1}));
    end
    ip.parse(varargin{:});
    options = ip.Results;
    % additional correction to the inputs
    if ~isempty( options.axesHdl)
        options.plotModel = 1;   
    end
    
    do_gauss=options.do_gauss;
    makeMovie=options.makeMovie;
    timePoint_step=options.timePoint_step;
    doSave=options.doSave;
    
    reg=options.reg;
    regsc=options.regsc;
    
    doDef= options.doDef;       % boolen to choose if the deformation part should be done or not
    vx          = options.vx;   % Voxel sizes
    settings_v  = options.settings_v; % Regularisation for displacements
    alpha       = options.alpha; % Update step (can be reduced if updates overshoot)

 if makeMovie
    writerObj = VideoWriter([outputFolder filesep 'movie.avi']); % Name it.
    writerObj.FrameRate = 1; % How many frames per second.
    open(writerObj); 
 end 

 
%% --read inputs
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


%% -- figures
% figure position and settings
% screenSize=get(0,'ScreenSize');
xCoor(1)=0; xCoor(2)=1/2;
yCoor=1/4; 
xSize=1/2; ySize=2/3;
tagNames={'model','def'};
% axes position
dd2=0.035; dxlin2=(1-(dd2*4))/3;
dd3=0.1; dxlin3=(1-(dd3*3))/2;
axesHdl=cell(1,numel(xCoor));

for iFig=1:numel(xCoor)
    figHdl(iFig)=figure('units','normalized','Tag',tagNames{iFig},'Position',[xCoor(iFig) yCoor xSize ySize ]);
%     figHdl(iFig)=figure('units','normalized','Tag',tagNames{iFig},'Position',[xCoor(iFig) yCoor xSize ySize ],'Visible','off');
    
    axesHdl{iFig}(1)=axes('Parent',figHdl(iFig),'Tag', 'alphas','Unit','Normalized','Position', [1/2-1/3 2/3 2/3 1/3-dd2]);
%     textHdl(1)=title('alpha_M_T - alpha_P_D -alpha_T_1');    
    axesHdl{iFig}(end+1)=axes('Parent',figHdl(iFig),'Tag', 'hisMT','Unit','Normalized','Position', [dd2               1/3 dxlin2 dxlin2]); % for plotting the joint histogram to PD
    axesHdl{iFig}(end+1)=axes('Parent',figHdl(iFig),'Tag', 'hisT1','Unit','Normalized','Position', [2*dd2+dxlin2      1/3 dxlin2 dxlin2]);
    axesHdl{iFig}(end+1)=axes('Parent',figHdl(iFig),'Tag', 'RGB','Unit','Normalized','Position',   [3*dd2+2*dxlin2    1/3 dxlin2 dxlin2]);
    axesHdl{iFig}(end+1)=axes('Parent',figHdl(iFig),'Tag', 'ldef','Unit','Normalized','Position', [dd3           0.04 dxlin3 1/3-dd3-0.05]);
    axesHdl{iFig}(end+1)=axes('Parent',figHdl(iFig),'Tag', 'ltot','Unit','Normalized','Position', [2*dd3+dxlin3  0.04 dxlin3 1/3-dd3-0.05]);

end

%% -- settings model (exp fit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sd    = spm_noise_estimate(thisSubjFiles); % here ask if it should be computed to each contrast or as 
% here use a common value
s     = mean(sd.^2);
Nii   = nifti(thisSubjFiles);
dm    = size(Nii(1).dat);
%range = {1:dm(1),1:dm(2),1:dm(3)};

% initialiaze
L_model = [];
L_def =[];
L_tot=[];
clear thisReadVol tmp
thisSubjFiles=cellstr(thisSubjFiles);

% reduce by a factor of 1:
titleSubplot={'beta', 'a MT', 'a PD' 'a T1'}; % title of subplots in the model loop 


%% -- settings deformation 
% doDef= 1; % boolen to choose if the deformation part should be done or not
% vx          = [1 1 1]; % Voxel sizes
% settings_v  = [vx  0.00001 0.01 10 0.05 0.1]; % Regularisation for displacements
% alpha       = 0.9; % Update step (can be reduced if updates overshoot)

% Construct an identity transform
% [x1,x2,x3] = ndgrid(single(1:dm(1)),single(1:dm(2)), single(1:dm(3)));
% x = cat(4,x1,x2,x3);

% Initial estimates of displacement fields
% at the end we should have 1 displacement field for each map, because we
% assume deformation is only between modalities and not/small inside
% modalities(within map the different echoes we assume no movement or very small)
V=cell(size(T));
for iMap=1:numel(V)% Loop over modalities
    ii=sum(TP_eachMap(1:iMap-1))+1;
    V{iMap} = zeros([Nii(ii).dat.dim 3],'single'); % ii ensures we are taking the size of the first echo of each map
end




%% --estimate Theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting estimates computed using the simple linear approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np    = 1+numel(T); % number of parameters
ii    = 0;
g     = zeros([dm, np],'single');         % Gradients
H     = zeros([dm,(1+np)*np/2],'single'); % Hessians

for iMap=1:timePoint_step:numel(T) % Loop over time series/ or maps
    msk = true(dm);
    for jTP=1:timePoint_step:numel(T{iMap}) % Loop over time points
        ii  = ii + timePoint_step;
        if jTP==1
            Mi  = Nii(ii).mat;
            M1  = Nii(1).mat;
            phi = AffMap(dm,M1\Mi);
            g_thisMap  = zeros([dm 2],'single');
            h_thisMap  = zeros([dm 3],'single');

        end
        y       = single(Nii(ii).dat(:,:,:,1));
        msk     = isfinite(y) & (y>0);
        g0      = zeros(dm,'single');
        h0      = zeros(dm,'single');
        g0(msk) = log(y(msk));
        h0(msk) = 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INCREMENTS
        t    = T{iMap}(jTP);

        % Increment gradients 
        g_thisMap(:,:,:,1) = g_thisMap(:,:,:,1) + g0*(-t);
        g_thisMap(:,:,:,2) = g_thisMap(:,:,:,2) + g0;

        % Increment Hessians 
        h_thisMap(:,:,:,1) = h_thisMap(:,:,:,1) + h0*(t^2);
        h_thisMap(:,:,:,2) = h_thisMap(:,:,:,2) + h0;
        h_thisMap(:,:,:,3) = h_thisMap(:,:,:,3) + h0*(-t);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('.');
    end

    % Push to common space
    spm_diffeo('boundary',1);            % Neumann boundary condition
    g_thisMap            = spm_diffeo('push',g_thisMap,phi);
    h_thisMap            = spm_diffeo('push',h_thisMap,phi);

    % Increment gradients 
    g(:,:,:,1)    = g(:,:,:,1)    + g_thisMap(:,:,:,1);
    g(:,:,:,iMap+1)  = g(:,:,:,iMap+1)  + g_thisMap(:,:,:,2);

    % Increment Hessians 
    H(:,:,:,1)    = H(:,:,:,1)    + h_thisMap(:,:,:,1);
    H(:,:,:,iMap+1)  = H(:,:,:,iMap+1)  + h_thisMap(:,:,:,2);
    H(:,:,:,iMap+np) = H(:,:,:,iMap+np) + h_thisMap(:,:,:,3);
end % end loop over maps

spm_field('boundary',1);                       % Neumann boundary condition
Theta = spm_field(H,g,[1 1 1  1e-6 0 0  1 1]); % Solve theta = H\g
Theta(:,:,:,1) = max(Theta(:,:,:,1),0);  % here set to 0 all what is smaller (= set to 0 neg values)
fprintf('\n');

clear g0 h0 g1 h1 H g Mi M1 phi iTP jTP ii 


%%  %%%%%%%%%% MAIN LOOP - GN %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The actual Gauss-Newton updates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nits = 5; % number of iterations
sliceIteraions=cell(1,nits); % for development purpose to keep one slice for each iteration
% Set boundary conditions 
%-----------------------------------------------------------------------
spm_field('boundary',1);        % Bias correction - Neumann
spm_diffeo('boundary',0);       % Diffeomorphism  - circulant

sliceToPlot = ceil(dm(3)/2);


%% rigid deformation 
% to ensure that the images are rigidly aligned before reconstructing the various parameter maps 
% here I am coregistering all echoes to PD ad refernece (as done in the new hMRI toolbox)
refIdx=2; % index 2 refers to PD
refImg=thisSubjFiles(find(mapIdx==2,1));

[mapIdx, mapPath]=grp2idx(mainMapPath);

job.ref=refImg;
job.eoptions=struct('cost_fun','nmi','sep',[4 2],'fwhm',[7 7]);
 
for iContrast=[1 3]
    this=mapIdx==iContrast;
    thisContrastToCoreg=thisSubjFiles(this);
    job.source=thisContrastToCoreg(1);
    job.other=thisContrastToCoreg(2:end);
    spm_run_coreg(job)
end

% main loop 
%-----------------------------------------------------------------------
for it=1:nits
   
    
    %% step 1- estimate Theta (model/gaussian model part)
    % model Estimation
    %=======================================================================
    % Recompute theta using GN
	%-----------------------------------------------------------------------
    % Ell  = 0;   % Log-likelihood -> % -log p(F|mu,V) = \sum_i log p(f_i | mu, v_i)
    
    [Ell,llr,Theta] = updateModelFit (Nii,mapIdx,s,Theta,T,Nii(1),V,reg,regsc,'axesHdl',axesHdl{1}(2:4),'timePoint_step',timePoint_step,'do_gauss',do_gauss);

    %% -- Some visualisation stuff
    fprintf('Model fit iter:%3d -  Ell:%g llr:%g Ell+llr:%g\n', it, Ell, llr, Ell+llr);
    L_model     = [L_model Ell+llr]; 
    thisAxes=axesHdl{1}(end-1);
    plot(thisAxes,L_model,'.-'); drawnow;

    if any(~isfinite(Theta(:))), crash; end
    thisAxes=axesHdl{1}(1);
    imagesc([Theta(:,:,sliceToPlot,2)',Theta(:,:,sliceToPlot,3)',Theta(:,:,sliceToPlot,4)'],...
        'parent',thisAxes ); 
    set(thisAxes,'DataAspectRatio' , [1 1 1],'DataAspectRatioMode','manual','Visible','off',...
        'XLimMode','auto','XDir','normal','YLimMode','auto','YDir','normal', 'ZLimMode','auto','PlotBoxAspectRatio',[300 80 1]);
%     thisAxes.Title.String='alpha_M_T - alpha_P_D -alpha_T_1';
    
    thisAxes=axesHdl{1}(end);
    imagesc(Theta(:,:,sliceToPlot,1)','parent',thisAxes ); 
    set(thisAxes,'DataAspectRatio' , [1 1 1],'DataAspectRatioMode','manual','Visible','off',...
        'XLimMode','auto','XDir','normal','YLimMode','auto','YDir','normal', 'ZLimMode','auto','PlotBoxAspectRatio',[100 50 1]);
    
    clear H g gr Hreg
%     for iSubplot=1:size(Theta,4)
%         imagesc(Theta(:,:,ceil(end/2),iSubplot)'); axis image xy off; colorbar
%         title(titleSubplot{iSubplot})
%     end
    drawnow
    
    %% make movie
    if makeMovie
        frame = getframe(figHdl(1)); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
%         saveas(figHdl,[ mainFolder filesep '3-figures' filesep 'boxplot_4Gpr' filesep parameters{iParam} '.eps'],'epsc')
        saveas(figHdl(1),[outputFolder filesep  'modelUpdate_iter_' num2str(it) '.png'])
    end
    

    %% step 2- update deformation
    if doDef
        % optionally deformation could also not be done
        % Deformations
        %=======================================================================
        % Minimises -log p(F,V,mu).  Note that some constants are
        % ignored as they don't feature in the objective function
        %-----------------------------------------------------------------------        
        [Ell, Ev,V] = updateDeformation (Nii,mapIdx,s,Theta,T,V,alpha,settings_v,'timePoint_step',timePoint_step)  ;      
        


        %% --Some visualisation stuff
        if it==1
            set(figHdl(2),'Visible','on')
        end
        fprintf('Deform Update iter:%3d  -  Ell:%g Ev:%g  Ell+Ev:%g\n', it, Ell, Ev, Ell+Ev);
                
%         thisAxes=findobj(figHdl(2),'Tag','ldef');
        thisAxes=axesHdl{2}(end-1);
        L_def     = [L_def Ell+Ev]; plot(thisAxes,L_def,'.-'); drawnow; 
        title(thisAxes,'Likelihood-Deformable Model')
        
%         thisAxes=findobj(figHdl(2),'Tag','ltot');
        thisAxes=axesHdl{2}(end);
        L_tot     = [L_tot Ell+Ev+llr]; plot(thisAxes,L_tot,'.-'); drawnow;
        title(thisAxes,'Likelihood-Total Model')

        % plot alphas
%         thisAxes=findobj(figHdl(2),'Tag','alphas');
        thisAxes=axesHdl{2}(1);
        imagesc([Theta(:,:,sliceToPlot,2)',Theta(:,:,sliceToPlot,3)',Theta(:,:,sliceToPlot,4)'],...
            'parent',thisAxes); 
        set(thisAxes,'DataAspectRatio' , [1 1 1],'DataAspectRatioMode','manual','Visible','off',...
            'XLimMode','auto','XDir','normal','YLimMode','auto','YDir','normal', 'ZLimMode','auto','PlotBoxAspectRatio',[300 80 1]);
        title(thisAxes,'a_M_T   a_P_D   a_T_1')
        
        % joint histogram all refered to PD (as done in unicort)
        MapIdx=[2 4]; % index for comparison with PD map
        for iHist=1:2
            thisAxes=axesHdl{2}(1+iHist);
            [N_iMap_pd,XEDGES_iMap,YEDGES_iMap] = histcounts2(Theta(:,:,sliceToPlot,MapIdx(iHist)),Theta(:,:,sliceToPlot,3));
            imagesc(N_iMap_pd,'parent',thisAxes)
    %         colormap gray
%             set(thisAxes,'DataAspectRatio' , [1 1 1],'DataAspectRatioMode','manual','Visible','off',...
%                 'XLimMode','auto','XDir','normal','YLimMode','auto','YDir','normal', 'ZLimMode','auto','PlotBoxAspectRatio',[1 1 1]);
            set(thisAxes,'Tag','RGB','DataAspectRatio' , [1 1 1],'Visible','off',...
            'XLimMode','auto', 'YLimMode','auto', 'ZLimMode','auto','DataAspectRatioMode','auto','YDir','normal'); 
        end
        
        
        % additional visual inspection for coregistration:
%         thisAxes=findobj(figHdl(2),'Tag','RGB');
        thisAxes=axesHdl{2}(4);
        rgb = cat(3,Theta(:,:,sliceToPlot,2)',Theta(:,:,sliceToPlot,3)',Theta(:,:,sliceToPlot,4)');
%         rgb = cat(3,Theta(:,:,sliceToPlot,2)',Theta(:,:,sliceToPlot,3),Theta(:,:,sliceToPlot,4));
        % keep one slice for offline comaprison between iterations
        sliceIteraions{it}=rgb;
        rgb = rgb - min(rgb(:));
        rgb = rgb/max(rgb(:));
        image(rgb,'parent',thisAxes);
%         axis image xy off;
%         title('check Registration-compare a_i')
        set(thisAxes,'Tag','RGB','DataAspectRatio' , [1 1 1],'Visible','off',...
            'XLimMode','auto', 'YLimMode','auto', 'ZLimMode','auto','DataAspectRatioMode','auto','YDir','normal'); 
                
        drawnow
        
        if makeMovie
%         frame = getframe(figHdl(1)); % 'gcf' can handle if you zoom in to take a movie.
%         writeVideo(writerObj, frame);
%         saveas(figHdl,[ mainFolder filesep '3-figures' filesep 'boxplot_4Gpr' filesep parameters{iParam} '.eps'],'epsc')
        saveas(figHdl(2),[outputFolder filesep  'deformationUpdate_iter_' num2str(it) '.png'])
        end

    else
        if it ==1; close(figHdl(2)); end
    end % end for if-> for doing deformation part
    
end % end main loop

% figure
% imshow3D(Theta(:,:,:,1))
% figure
% imshow3D(Theta(:,:,:,2))
% figure
% imshow3D(Theta(:,:,:,3))
% figure
% imshow3D(Theta(:,:,:,4))


% try to understand possible value for alpha and beta:
% beta=Theta(:,:,ceil(end/2),1);
% alpha1=Theta(:,:,ceil(end/2),2);
% alpha2=Theta(:,:,ceil(end/2),3);
% alpha3=Theta(:,:,ceil(end/2),4);


if makeMovie
    close(writerObj); % Saves the movie.
end
 
if doSave
    save([outputFolder filesep 'runned_MultiImage_code_' datestr(now,'yyyymmdd') '.mat'])
end

%%  %%%%%%%%%% SAVE MAPS %%%%%%%%%%

%% -- save Theta
mapsName={'beta_R2s','alpha_MT','alpha_PD','alpha_T1'};

for iMap=1:length(mapsName)
    thetaMaps= nifti;
    thetaMaps.mat = Nii(1).mat;
    thetaMaps.dat = file_array([outputFolder filesep mapsName{iMap} '.nii'], [size(Theta,1) size(Theta,2) size(Theta,3)], 'float32');
    % thetaMaps.dat = file_array([ mapsName{iMap} '.nii'], [size(Theta,1) size(Theta,2) size(Theta,3)], 'float32');
    create(thetaMaps);
    thetaMaps.dat(:,:,:) = Theta(:,:,:,iMap);
    thetaMaps_ALL{iMap}=thetaMaps;
    clear thetaMaps
end

%% -- save deformation 
% to be used as proof or evaluation of applied deformation

mapsName2={'def_MT','def_PD','def_T1'}; 

for iDef=1:length(mapsName2)
    defMaps= nifti;
    defMaps.mat = Nii(1).mat;
    defMaps.dat = file_array([outputFolder filesep mapsName2{iDef} '.nii'], size(V{iDef}), 'float32');
    % thetaMaps.dat = file_array([ mapsName{iMap} '.nii'], [size(Theta,1) size(Theta,2) size(Theta,3)], 'float32');
    create(defMaps);
    defMaps.dat(:,:,:,:) = V{iDef};

end


% optional save v here... to be done potentially useful for motion
% evaluation

disp('done saved maps in:' )
disp(outputFolder)


%% -- save maps
% create maps as in hMRI
% as in fcn:
% hmri_create_MTProt.m

% =======================================================================%
% Prepare output for R1, PD and MT maps
%=========================================================================%
mpm_params.output={'R1','A','MT'};
[path ,outbasename]=fileparts(Nii(1).dat.fname);
dm= Nii(1).dat.dim;

% define NIFTI objects for output images
Nmap    = nifti;
for ii=1:length(TP_eachMap) 
    %dm         = V_pdw(1).dim;
    Ni         = nifti;
    Ni.mat     = Nii(1).mat ; %V_pdw(1).mat;
    Ni.mat0    = Nii(1).mat ; %V_pdw(1).mat;
    Ni.descrip = 'generate map from generative model';
    Ni.dat     = file_array(fullfile(outputFolder,[outbasename '_' mpm_params.output{ii} '.nii']),dm,'float32');
    create(Ni);
    Nmap(ii) = Ni;
end

% =======================================================================%
% Map calculation continued (R1, PD, MT) 
%=========================================================================%


fprintf(1,'\n    -------- Map calculation continued (R1, PD, MT) --------\n');

M0 = Nii(1).mat;

for iMap=1:length(TP_eachMap) 
    fa_rad{iMap}=fa{iMap}* pi / 180;
end
% First calculate R1 & MTR

%  remember theta has at theta(:,:,:, 2:end)has in the following order: MT, Pd, T1
idx_MT=1;
idx_PD=2;
idx_T1=3;
mpm_params.interp=3;

fa_t1w_rad=fa_rad{idx_T1};
fa_pdw_rad=fa_rad{idx_PD};
fa_mtw_rad=fa_rad{idx_MT};
TR_t1w=tr{idx_T1};
TR_pdw=tr{idx_PD};
TR_mtw=tr{idx_MT};
% thresholds for maps
threshall.R1  =   2000;
threshall.A = 100000;
threshall.MT =5;

%% -- compute map R1

for p = 1:dm(3)
    M = M0*spm_matrix([0 0 p]);

    % PDw images are always available, so this bit is always loaded:
    PDw = spm_slice_vol(Theta(:,:,:,idx_PD),thetaMaps_ALL{idx_PD}.mat\M,dm(1:2),mpm_params.interp);
    
%     if ~isempty(V_trans)
%         f_T = spm_slice_vol(V_trans(2,:),V_trans(2,:).mat\M,dm(1:2),mpm_params.interp)/100; % divide by 100, since p.u. maps
%     else
        f_T = [];
%     end
    
%     % Standard magnetization transfer ratio (MTR) in percent units [p.u.]
%     % only if trpd = trmt and fapd = famt and if PDw and MTw images are
%     % available
%     if (MTwidx && PDwidx)
%         MTw = spm_slice_vol(Vavg(MTwidx),Vavg(MTwidx).mat\M,dm(1:2),mpm_params.interp);
%         if (TR_mtw == TR_pdw) && (fa_mtw == fa_pdw) % additional MTR image...
%             MTR = (PDw-MTw)./(PDw+eps) * 100;
%             % write MTR image
%             Nmap(mpm_params.qMTR).dat(:,:,p) = max(min(MTR,threshall.MTR),-threshall.MTR);
%         end          
%     end
    
    % T1 map and A/PD maps can only be calculated if T1w images are
    % available:
    if idx_T1

        T1w = spm_slice_vol(Theta(:,:,:,idx_T1),thetaMaps_ALL{idx_T1}.mat\M,dm(1:2),mpm_params.interp);


        if isempty(f_T)
            % semi-quantitative T1
            R1 = (((T1w * (fa_t1w_rad / 2 / TR_t1w)) - (PDw * fa_pdw_rad / 2 / TR_pdw)) ./ ...
                max(((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)),eps))*10^6;
        else
%             % Transmit bias corrected quantitative T1 values
%             % correct T1 for transmit bias f_T with fa_true = f_T * fa_nom
%             % T1corr = T1 / f_T / f_T
%             
%             if RFC.RFCorr
%                 % MFC: We do have P2_a and P2_b parameters for this sequence
%                 % => T1 = A(B1) + B(B1)*T1app (see Preibisch 2009)
%                 T1 = RFC.P2_a(1)*f_T.^2 + ...
%                     RFC.P2_a(2)*f_T + ...
%                     RFC.P2_a(3) + ...
%                     (RFC.P2_b(1)*f_T.^2+RFC.P2_b(2)*f_T+RFC.P2_b(3)) .* ...
%                     ((((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)+eps) ./ ...
%                     max((T1w * fa_t1w_rad / 2 / TR_t1w) - (PDw * fa_pdw_rad / 2 / TR_pdw),eps))./f_T.^2);
%             else
%                 % MFC: We do not have P2_a or P2_b parameters for this sequence
%                 % => T1 = T1app
%                 T1 = ((((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)+eps) ./ ...
%                     max((T1w * fa_t1w_rad / 2 / TR_t1w) - (PDw * fa_pdw_rad / 2 / TR_pdw),eps))./f_T.^2);
%             end
%             
%             R1 = 1./T1*10^6;
        end
        
        R1(R1<0) = 0;
        tmp      = R1;
        Nmap(1).dat(:,:,p) = min(max(tmp,-threshall.R1),threshall.R1)*0.001; % truncating images
                
    end
    spm_progress_bar('Set',p);
end

%% -- compute map MT & PD

for p = 1:dm(3)
    M = M0*spm_matrix([0 0 p]);

%     if ~isempty(V_trans)
%         f_T = spm_slice_vol(V_trans(2,:),V_trans(2,:).mat\M,dm(1:2),mpm_params.interp)/100; % divide by 100, since p.u. maps
%     elseif ~isempty(V_trans_unicort)
%         f_T = spm_slice_vol(V_trans_unicort(1,:),V_trans_unicort(1,:).mat\M,dm(1:2),mpm_params.interp)/100; % divide by 100, since p.u. maps        
%     else
        f_T = [];
%     end
%     
    % PDw images are always available, so this bit is always loaded:
    PDw = spm_slice_vol(Theta(:,:,:,idx_PD),thetaMaps_ALL{idx_PD}.mat\M,dm(1:2),mpm_params.interp);
    
    % T1 map and A/PD maps can only be calculated if T1w images are
    % available:
    if idx_T1
        V_R1=spm_vol(thetaMaps_ALL{idx_T1}.dat.fname) ; % or the masked one (R1_unicort)
        T1w = spm_slice_vol(Theta(:,:,:,idx_T1),thetaMaps_ALL{idx_T1}.mat\M,dm(1:2),mpm_params.interp);
        T1 = 10^3./spm_slice_vol(V_R1,V_R1.mat\M,dm(1:2),mpm_params.interp);
        
        % if "fullOLS" option enabled, use the OLS fit at TE=0 as
        % "T1w_forA"; otherwise use the average calculated earlier (by
        % default, corresponds to the first echo to reduce R2* bias)
%         if mpm_params.fullOLS
            T1w_forA = T1w;
%         else
%             T1w_forA = spm_slice_vol(VT1w_forA,VT1w_forA.mat\M,dm(1:2),mpm_params.interp);
%         end
%                 
        % A values proportional to PD
        % f_T correction is applied either if:
        % - f_T has been provided as separate B1 mapping measurement (not
        % UNICORT!) or
        % - f_T has been calculated using UNICORT *AND* the UNICORT.PD flag
        % is enabled (advanced user only! method not validated yet!)
%         if(~isempty(f_T))&&(~mpm_params.UNICORT.R1 || mpm_params.UNICORT.PD)
%             A = T1 .* (T1w_forA .*(fa_t1w_rad*f_T) / 2 / TR_t1w) + (T1w_forA ./ (fa_t1w_rad*f_T));
%         else
            % semi-quantitative A
            A = T1 .* (T1w_forA * fa_t1w_rad / 2 / TR_t1w) + (T1w_forA / fa_t1w_rad);
%         end
        
        tmp      = A;
        Nmap(2).dat(:,:,p) = max(min(tmp,threshall.A),-threshall.A);
        % dynamic range increased to 10^5 to accommodate phased-array coils and symmetrical for noise distribution

        % for MT maps calculation, one needs MTw images on top of the T1w
        % and PDw ones...
%         if MTwidx
            MTw = spm_slice_vol(Theta(:,:,:,idx_MT),thetaMaps_ALL{idx_MT}.mat\M,dm(1:2),3);
            T1_forMT = ((PDw / fa_pdw_rad) - (T1w / fa_t1w_rad)) ./ ...
                max((T1w * (fa_t1w_rad / 2 / TR_t1w)) - (PDw * fa_pdw_rad / 2 / TR_pdw),eps);
            A_forMT = T1_forMT .* (T1w * fa_t1w_rad / 2 / TR_t1w) + (T1w / fa_t1w_rad);
            
            % MT in [p.u.]; offset by - famt * famt / 2 * 100 where MT_w = 0 (outside mask)
            MT       = ( (A_forMT * fa_mtw_rad - MTw) ./ (MTw+eps) ./ (T1_forMT + eps) * TR_mtw - fa_mtw_rad^2 / 2 ) * 100;
            % f_T correction is applied either if:
            % - f_T has been provided as separate B1 mapping measurement (not
            % UNICORT!) or
            % - f_T has been calculated using UNICORT *AND* the UNICORT.MT flag
            % is enabled (advanced user only! method not validated yet!)
%             if (~isempty(f_T))&&(~mpm_params.UNICORT.R1 || mpm_params.UNICORT.MT)
%                 MT = MT .* (1 - 0.4) ./ (1 - 0.4 * f_T);
%             end
            
            tmp      = MT;
            Nmap(3).dat(:,:,p) = max(min(tmp,threshall.MT),-threshall.MT);
        end
end

  



output_args=[]% decide to the outputs, not done jet
end % end fcn 







