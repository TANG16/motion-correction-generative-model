% here I use only a portion of the intere FOV and less TP for speeding
% reason

clear 
close all
clc

%% step 0: initializzation
% please for going back to full TP & full image check following settings:
timePoint_step=2; % used in the loop for looping over all TP or using step/increment 'timePoint_step'
reduceSize=1; % reduce size of the echoes maps using only lower part of images -> used data with reduced size in the appropiate folder

addpath(fullfile(pwd,'utilities'))
addpath(fullfile(pwd,'utilities','auxiliary-functions-UCL'))

do_gauss = true; % boolean for choosing the model to use(Gaussian or Rician)

%% --input
% input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% from John
% P = spm_select(22,'nifti','Select the raw MRI');
% T    = {[2.3 4.6 6.9 9.2 11.5 13.8 16.1 18.4],[2.3 4.6 6.9 9.2 11.5 13.8 16.1 18.4],[2.3 4.6 6.9 9.2 11.5 13.8]}; 

% here the same from MA and real T from images
[mainPath, folderName]=fileparts(pwd);
if reduceSize
    dataFolder='0-data_ReducedSize'; % folder with reduced size images created in advance
    filterStr='^reduced_subj.*';
else
    dataFolder='0-data';
    filterStr='^subj.*';
end
subjDataFolder=fullfile(mainPath,folderName,dataFolder,'CL01_Day00');
[pathData, subjFileName]=fileparts(subjDataFolder);
% output folder
outputFolder=fullfile(mainPath,folderName,['1-outputCoregistr_Red3Regsc_' datestr(now,'yyyymmdd')],subjFileName);
if ~exist(outputFolder,'dir')
    sts = mkdir(outputFolder);
    if ~sts, error('Error creating output directory "%s".',outputFolder); end
end 

% read inputs:
% here automatic filled data, in the future from user in a batche format:

P=spm_select('FPListRec',subjDataFolder,filterStr);
if isempty(P)
    error('no selected data: please check path to data')
end
hdr=spm_vol(P);

% for grouping based on the folder path:
mainMapPath=cellfun(@(tmp)fileparts(tmp ),cellstr( P),'UniformOutput', false);
[mapIdx, mapPath]=grp2idx(mainMapPath);
TP_eachMap=arrayfun(@(tmp)sum(mapIdx==tmp),1:length(mapIdx));
TP_eachMap=TP_eachMap(~TP_eachMap==0);

ii=0;
T=cell(numel(TP_eachMap),1);    

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
end
%% -- figures
figHdl(1)=figure('Tag','model','Position',[43  296 709 472 ],'Visible','off');
figHdl(1)=figure('Tag','model','Position',[43  296 709 472 ]);
figHdl(2)=figure('Tag','def','Position',[754   295 709 472 ],'Visible','off');
% figHdl(2)=figure('Tag','def','Position',[754   295 709 472 ]);

% create nice position for each plot:
% settings for second line plot:
dd2=0.035; dxlin2=(1-(dd2*4))/3;
% settings for last line plot:
dd3=0.1; dxlin3=(1-(dd3*3))/2;

% figure1
% aaf1(1)=axes('Parent',figHdl(1),'Tag', 'alphas','Unit','Normalized','Position', [1/2-1/3 2/3 2/3 1/3-dd2]);
% 
% aaf1(2)=axes('Parent',figHdl(1),'Tag', 'hisMT','Unit','Normalized','Position', [dd2               1/3 dxlin2 dxlin2]); % for plotting the joint histogram to PD
% aaf1(3)=axes('Parent',figHdl(1),'Tag', 'hisT1','Unit','Normalized','Position', [2*dd2+dxlin2      1/3 dxlin2 dxlin2]);
% aaf1(4)=axes('Parent',figHdl(1),'Tag', 'RGB','Unit','Normalized','Position',   [3*dd2+2*dxlin2    1/3 dxlin2 dxlin2]);
% 
% aaf1(5)=axes('Parent',figHdl(1),'Tag', 'ldef','Unit','Normalized','Position', [dd3           0.04 dxlin3 1/3-dd3-0.05]);

%figure2
aa(1)=axes('Parent',figHdl(2),'Tag', 'alphas','Unit','Normalized','Position', [1/2-1/3 2/3 2/3 1/3-dd2]);

aa(2)=axes('Parent',figHdl(2),'Tag', 'hisMT','Unit','Normalized','Position', [dd2               1/3 dxlin2 dxlin2]); % for plotting the joint histogram to PD
aa(3)=axes('Parent',figHdl(2),'Tag', 'hisT1','Unit','Normalized','Position', [2*dd2+dxlin2      1/3 dxlin2 dxlin2]);
aa(4)=axes('Parent',figHdl(2),'Tag', 'RGB','Unit','Normalized','Position',   [3*dd2+2*dxlin2    1/3 dxlin2 dxlin2]);

aa(5)=axes('Parent',figHdl(2),'Tag', 'ldef','Unit','Normalized','Position', [dd3           0.04 dxlin3 1/3-dd3-0.05]);
aa(6)=axes('Parent',figHdl(2),'Tag', 'ltot','Unit','Normalized','Position', [2*dd3+dxlin3  0.04 dxlin3 1/3-dd3-0.05]);

%% -- settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sd    = spm_noise_estimate(P);
s     = mean(sd.^2);
Nii   = nifti(P);
dm    = size(Nii(1).dat);
%range = {1:dm(1),1:dm(2),1:dm(3)};

L_model = [];
L_def =[];
L_tot=[];
clear thisReadVol tmp P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for spm_field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters are voxel sizes (x, y and z) followed by 3 regularisation settings
% and then the number of iterations and multi-grid cycles.
% If no smoothing is required, then this bit of code could be speeded up lots.
% Some form of total variation regularisation could be used, in which case
% the smoothing part would be useful.
reg      = [1 1 1  [0 1 0]];
% note cahnge the full regsc-> look if the smooth is too much, adding
% artifacts
regsc    = ones(numel(T)+1,1)*1000; % intercept regulariztion
regsc(1) = 100000; % regulraz of decay map
% reduce by a factor of 1:
regsc=regsc/1000;
titleSubplot={'beta', 'a MT', 'a PD' 'a T1'}; % title of subplots in the model loop 

%% --estimate Theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting estimates computed using the simple linear approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np    = 1+numel(T); % number of parameters
ii    = 0;
g     = zeros([dm, np],'single');         % Gradients
H     = zeros([dm,(1+np)*np/2],'single'); % Hessians

for iMap=1:numel(T) % Loop over time series/ or maps
    msk = true(dm);
    for jTP=1:timePoint_step:numel(T{iMap}) % Loop over time points
        ii  = ii + 1;
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

%% -- deformation settings
doDef= 1; % boolen to choose if the deformation part should be done or not
vx          = [1 1 1]; % Voxel sizes
settings_v  = [vx  0.00001 0.01 10 0.05 0.1]; % Regularisation for displacements
alpha       = 0.9; % Update step (can be reduced if updates overshoot)

% Construct an identity transform
[x1,x2,x3] = ndgrid(single(1:dm(1)),single(1:dm(2)), single(1:dm(3)));
x = cat(4,x1,x2,x3);

% Initial estimates of displacement fields
% at the end we should have 1 displacement field for each map, because we
% assume deformation is only between modalities and not/small inside
% modalities(within map the different echoes we assume no movement or very small)
V=cell(size(T));
for iMap=1:numel(V)% Loop over modalities
    if iMap==1, ii=1; 
    else 
        ii=sum(TP_eachMap(1:iMap-1))+1; 
    end
    V{iMap} = zeros([Nii(ii).dat.dim 3],'single'); % ii ensures we are taking the size of the first echo of each map
end



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

% main loop 
%-----------------------------------------------------------------------
for it=1:nits
   
    %% step 1- estimate Theta (model/gaussian model part)
    % model Estimation
    %=======================================================================
    % Recompute theta using GN
	%-----------------------------------------------------------------------
    Ell  = 0;   % Log-likelihood -> % -log p(F|mu,V) = \sum_i log p(f_i | mu, v_i)

    
    g   = zeros([dm, np],'single');         % Gradients
    H   = zeros([dm,(np+1)*np/2],'single'); % Hessians
    ii  = 0;
    for iMap=1:numel(T) % Loop over time series
        v=V{iMap};
        Nii_this=Nii(mapIdx==iMap);
        timeVector=T{iMap};
        theta_this=cat(4,Theta(:,:,:,1),Theta(:,:,:,iMap+1));
        [Ell ,phi, g_thisMap, h_thisMap]=estimateGradHessEll(Nii_this,s,theta_this,timeVector,x+v,Nii(1));

        % Push to common space (the motion corrected one)
        g_thisMap            = spm_diffeo('push',g_thisMap,phi);
        h_thisMap            = spm_diffeo('push',h_thisMap,phi);

        % Increment gradients 
        g(:,:,:,1)    = g(:,:,:,1)    + g_thisMap(:,:,:,1);
        g(:,:,:,iMap+1)  = g(:,:,:,iMap+1)  + g_thisMap(:,:,:,2);

        % Increment Hessians 
        H(:,:,:,1)    = H(:,:,:,1)    + h_thisMap(:,:,:,1);
        H(:,:,:,iMap+1)  = H(:,:,:,iMap+1)  + h_thisMap(:,:,:,2);
        H(:,:,:,np+iMap) = H(:,:,:,np+iMap) + h_thisMap(:,:,:,3);
    end
    
    gr     = spm_field('vel2mom',Theta,reg,regsc); 

    llr    = -0.5*sum(sum(sum(sum(gr.*Theta))));

    %% -- Some visualisation stuff
    figure(figHdl(1))
    fprintf('Model fit iter:%3d -  Ell:%g llr:%g Ell+llr:%g\n', it, Ell, llr, Ell+llr);
    L_model     = [L_model Ell+llr]; subplot(4,2,7); plot(L_model,'.-'); drawnow;

    % Add something to the Hessians for stability
    Hreg    = sum(H(:,:,:,1:size(Theta,4)),4)*(sqrt(eps('single')));
    H(:,:,:,1:size(Theta,4)) = bsxfun(@plus,H(:,:,:,1:size(Theta,4)),Hreg);

    Theta = Theta - spm_field(H,g+gr,[reg 1 1],regsc);  % Solve theta = theta - H\g
    Theta(:,:,:,1) = max(Theta(:,:,:,1),0);             % set to 0 negative values

    if any(~isfinite(Theta(:))), crash; end
    clear H g gr Hreg
    for iSubplot=1:size(Theta,4)
        subplot(4,2,iSubplot+2);
        imagesc(Theta(:,:,ceil(end/2),iSubplot)'); axis image xy off; colorbar
%         title(titleSubplot{iSubplot})
    end
    drawnow

    %% step 2- update deformation
    if doDef
        % optionally deformation could also not be done

        % Deformations
        %=======================================================================
        % Minimises -log p(F,V,mu).  Note that some constants are
        % ignored as they don't feature in the objective function
        %-----------------------------------------------------------------------

        Ell = 0; % -log p(F|mu,V) = \sum_i log p(f_i | mu, v_i)
        Ev  = 0; % -log p(V)      = \sum_i log p(v_i)         

        % Loop over images, updating the displacements
        for iMap=1:numel(T) % loop over modalities

           %% -- compute g_map and H_map
            v=V{iMap};
            Nii_this=Nii(mapIdx==iMap);
            timeVector=T{iMap};
            theta_this=cat(4,Theta(:,:,:,1),Theta(:,:,:,iMap+1));       

            [Ell ,phi, g_thisMap, H_thisMap]=estimateGradHessEll(Nii_this,s,theta_this,timeVector,x+v,Nii(1)); % add input for not plotting


            %% -- deformation convert
           % here loop over TP not needed because the coregistration is between
           % modalities -> therefore the reference image is the first of each modality

            u   = spm_diffeo('vel2mom',v,settings_v);

            % Convert a velocity field to a momentum field by u = A*v, where
            % A is the large sparse matrix encoding some form of regularisation.
            % v and m are single precision floating point.

            Ev  = Ev  + 0.5*sum(sum(sum(sum(v.*u))));

            [~,gax,gay,gaz] = spm_diffeo('bsplins',Theta(:,:,:,iMap+1),phi,[1 1 1  1 1 1]);
            [~,gbx,gby,gbz] = spm_diffeo('bsplins',Theta(:,:,:,1),phi,[1 1 1  1 1 1]);

    %         including deformations
            % bsxfun(@times,dL/da, da/dv)-> with:
            %   dL/da = g_map(a_thisMap) 
            %   da/dv = [gax gay gaz] = cat(4,gax,gay,gaz) -> [dim 3]
            ga=cat(4,gax,gay,gaz);
            gb=cat(4,gbx,gby,gbz);
            gdef   = bsxfun(@times,g_thisMap(:,:,:,2),ga ) + bsxfun(@times,g_thisMap(:,:,:,1),gb );  

            Hdef = cat(4, H_thisMap(:,:,:,1).*gbx.*gbx + H_thisMap(:,:,:,2).*gax.*gax + H_thisMap(:,:,:,3).*(gbx.*gax + gax.*gbx), ...
                       H_thisMap(:,:,:,1).*gby.*gby + H_thisMap(:,:,:,2).*gay.*gay + H_thisMap(:,:,:,3).*(gby.*gay + gay.*gby), ...
                       H_thisMap(:,:,:,1).*gbz.*gbz + H_thisMap(:,:,:,2).*gaz.*gaz + H_thisMap(:,:,:,3).*(gbz.*gaz + gaz.*gbz), ...
                       H_thisMap(:,:,:,1).*gbx.*gby + H_thisMap(:,:,:,2).*gax.*gay + H_thisMap(:,:,:,3).*(gbx.*gay + gax.*gby), ...
                       H_thisMap(:,:,:,1).*gbx.*gbz + H_thisMap(:,:,:,2).*gax.*gaz + H_thisMap(:,:,:,3).*(gbx.*gaz + gax.*gbz), ...
                       H_thisMap(:,:,:,1).*gby.*gbz + H_thisMap(:,:,:,2).*gay.*gaz + H_thisMap(:,:,:,3).*(gby.*gaz + gay.*gbz));

%% -- apply affine rotation to gradients    
        % rotate gradient such that g1=M(1:3,1:3)'*g0
        Affine=Nii(1).mat\Nii_this(1).mat;
        Jz = reshape(Affine(1:3,1:3),[1 1 1 3 3]);
        
        g   = zeros([dm(1:3),3],'single');
        H   = zeros([dm(1:3),6],'single');
        lkp    = [1 4 5; 4 2 6; 5 6 3];

        for z=1:dm(3)
            % Rotate gradients, such that g1 = J'*g0;
            for d1=1:3
                tmp = 0;
                for d2=1:3
                    tmp = tmp + Jz(:,:,:,d2,d1).*gdef(:,:,z,d2);
                end
                g(:,:,z,d1) = tmp;
            end

            % Rotate Hessian, such that H2 = J'*H0*J
            % First do H1 = J'*H0
            RH  = zeros([dm(1:2),1,3,3],'single');
            for d1=1:3
                for d3=1:3
                    tmp = 0;
                    for d2=1:3
                    	tmp = tmp + Jz(:,:,:,d2,d1).*Hdef(:,:,z,lkp(d2,d3));
                    end
                    RH(:,:,:,d1,d3) = tmp;
                end
            end

            % Then do H2 = H1*J
            for d1=1:3
                for d3=d1:3 % Only need compute an upper or lower triangle
                    tmp = 0;
                    for d2=1:3
                        tmp = tmp + RH(:,:,:,d1,d2).*Jz(:,:,:,d2,d3);
                    end
                    H(:,:,z,lkp(d1,d3)) = tmp;
                end
            end
        end

        Hdef=H; gdef=g;
        clear H g
        
        %% --GN update                   
            % Gauss-Newton update. Note that it is not guaranteed to improve things
            v_old=v;
            dv = spm_diffeo(Hdef,gdef+u,[settings_v 3 3]);
            alphaThis=alpha;

            for icheck=1:6 
                % check if things are improving otherwise: 
                %  go for smaller alpha, if still not improve that do not
                %  update, it means we are already very close to ideal/optimum

                v = v_old -alpha*dv;

                %  try with smaller alpha
                [Ell_new]=estimateGradHessEll(Nii_this,s,theta_this,timeVector,x+v,Nii(1)); % add input for not plotting
                Ev_new  = Ev  + 0.5*sum(sum(sum(sum(v.*u))));
                L_thisModel=Ell_new+Ev_new;
                if L_thisModel > Ev+Ell
                    if icheck==6
                        v = v_old; % if any alpha is not improving our model, keep old value
                    end
                    break
                end
                alphaThis=alphaThis/2;
            end

            % here please chec for improvements otherwise go back to prev value
            % and try with a different alpha??
            V{iMap} = v;

    %         imagesc([f(:,:,ceil(end/2))' thisModel(:,:,ceil(end/2))' a(:,:,ceil(end/2))']); axis image xy off; drawnow
    %         title('updateDef: acquired Image - model - difference(acquired Image-model)')

        end % end loop over modalities  
        clear g_thisMap H_thisMap Mi M1 phi th1 th2 f1 f2 g0 h0

        %% --Some visualisation stuff
        sliceToPlot = ceil(dm(3)/2);
        if it==1
            set(figHdl(2),'Visible','on')
        end
        fprintf('Deform Update iter:%3d  -  Ell:%g Ev:%g  Ell+Ev:%g\n', it, Ell, Ev, Ell+Ev);
        thisAxes=findobj(figHdl(2),'Tag','ldef');
        L_def     = [L_def Ell+Ev]; plot(thisAxes,L_def,'.-'); drawnow; 
        thisAxes.Tag='ldef';
        title(thisAxes,'Likelihood-Deformable Model')
        
        thisAxes=findobj(figHdl(2),'Tag','ltot');
        L_tot     = [L_tot Ell+Ev+llr]; plot(thisAxes,L_tot,'.-'); drawnow;
        thisAxes.Tag='ltot';
        title(thisAxes,'Likelihood-Total Model')

        % plot alphas
        thisAxes=findobj(figHdl(2),'Tag','alphas');
        imagesc([Theta(:,:,sliceToPlot,2)',Theta(:,:,sliceToPlot,3)',Theta(:,:,sliceToPlot,4)'],...
            'parent',thisAxes); 
%         axis image xy off; 
        set(thisAxes,'Tag','alphas','DataAspectRatio' , [1 1 1],'Visible','off',...
            'XLimMode','auto', 'YLimMode','auto', 'ZLimMode','auto','DataAspectRatioMode','auto','YDir','normal'); %
        title(thisAxes,'a_M_T   a_P_D   a_T_1')
        
        % joint histogram all refered to PD (as done in unicort)
%         thisAxes=findobj(figHdl(2),'Tag','hisMT');
        thisAxes=findobj(figHdl(2),'Tag','hisMT');
        [N_mt_pd,XEDGES_mt,YEDGES_mt] = histcounts2(Theta(:,:,sliceToPlot,2),Theta(:,:,sliceToPlot,3));
        imagesc(N_mt_pd,'parent',thisAxes)
        colormap gray
        set(thisAxes,'Tag','hisMT','DataAspectRatio' , [1 1 1],'Visible','off',...
            'XLimMode','auto', 'YLimMode','auto', 'ZLimMode','auto','DataAspectRatioMode','auto','YDir','normal'); %
        
        thisAxes=findobj(figHdl(2),'Tag','hisT1');
        [N_t1_pd,XEDGES_t1,YEDGES_t1] = histcounts2(Theta(:,:,sliceToPlot,4),Theta(:,:,sliceToPlot,3));
        imagesc(N_t1_pd,'parent',thisAxes)
        set(thisAxes,'Tag','hisT1','DataAspectRatio' , [1 1 1],'Visible','off',...
            'XLimMode','auto', 'YLimMode','auto', 'ZLimMode','auto','DataAspectRatioMode','auto','YDir','normal'); %
        
        % additional visual inspection for coregistration:
        thisAxes=findobj(figHdl(2),'Tag','RGB');
        
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

    else
        close(figHdl(2))
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

disp('done')
save([pwd filesep 'runned_MultiImage_code_' datestr(now,'yyyymmdd') '.mat'])

% try to understand possible value for alpha and beta:
beta=Theta(:,:,ceil(end/2),1);
alpha1=Theta(:,:,ceil(end/2),2);
alpha2=Theta(:,:,ceil(end/2),3);
alpha3=Theta(:,:,ceil(end/2),4);





%% step3: save maps
mapsName={'beta_R2s','alpha_MT','alpha_PD','alpha_T1'};

for iMap=1:length(mapsName)
    thetaMaps= nifti;
    thetaMaps.mat = Nii(1).mat;
    thetaMaps.dat = file_array([outputFolder filesep mapsName{iMap} '.nii'], [size(Theta,1) size(Theta,2) size(Theta,3)], 'float32');
    % thetaMaps.dat = file_array([ mapsName{iMap} '.nii'], [size(Theta,1) size(Theta,2) size(Theta,3)], 'float32');
    create(thetaMaps);
    thetaMaps.dat(:,:,:) = Theta(:,:,:,iMap);

    clear thetaMaps
end









