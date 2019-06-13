function [Ell, llr,Theta] = updateModelFit (Nii,mapIdx,s,Theta,T,refNii,V,reg,regsc,varargin)
% the fcn will return the model fit parameters
%
% SYNOPSIS: []=updateModelFit (varagin)
%
% INPUT:
%	Nii: array of nifti object(i.e all echoes of a map)
%       do_gauss = boolean (1 fo gaussian model, 0 for Rice ) -> Rice model
%       need to be tested
%   mapIdx: grouping index for contrasts   
%   s: variance of the model
%   Theta: parameters map each map parameter along the 4th dimension
%       (i.e. parameterMaps(:,:,:,1)=beta; parameterMaps(:,:,:,2)=alpha_mt
%       etc...)
%   timeVector: vector with time points should be one for each Nii
%   T: transformation needed to warp the model to this time point map
%   V: velocity fields
%   refNii: nifti file to use as reference for affine
% 
% OPTIONAL INPUTS: 
% look 'estimateGradHessEll' fcn
% thisAxes axes for each map
%
% OUTPUT the model fit parameters
%   Ell: conditional probability based on the model
%   llr: probability of the likelihood
%   Theta: model fit parameters Theta(:,:,:,1)=beta;
%       Theta(:,:,:,2:4)=alphas for each contrast
%
% REMARKS
%
% SEE ALSO Affscale,estimateGradHessEll
%
% EXAMPLE
%

% created with MATLAB ver.: 9.1.0.441655 (R2016b) on Microsoft Windows 7 Professional  Version 6.1 (Build 7601: Service Pack 1)
%
% created by: mazzarito
% DATE: 06-Sep-2018
%
% Last revision $Rev: 3231 $ $Date: 2013-02-06 12:35:49 -0500 (Wed, 06 Feb 2013) $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% input/output 
    % default inputs:
    def.timePoint_step = 1;
    def.plotModel = 0;
    def.figureHdl=[]; % by default all kind of hv will be considered 
    def.axesHdl=[]; % by default all kind of hv will be considered 
    def.do_gauss=1;

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


    % initialization
    Ell=0;  % Log-likelihood -> % -log p(F|mu,V) = \sum_i log p(f_i | mu, v_i)
    dm=Nii(1).dat.dim; % assuming all images of same size
    np=numel(T)+1;

    %% intialize
    
    % Construct an identity transform
    [x1,x2,x3] = ndgrid(single(1:dm(1)),single(1:dm(2)), single(1:dm(3)));
    x = cat(4,x1,x2,x3);

    g   = zeros([dm, np],'single');         % Gradients
    H   = zeros([dm,(np+1)*np/2],'single'); % Hessians

    for iMap=1:numel(T) % Loop over time series
        v=V{iMap};
        Nii_this=Nii(mapIdx==iMap);
        timeVector=T{iMap};
        theta_this=cat(4,Theta(:,:,:,1),Theta(:,:,:,iMap+1));
        thisAxes=options.axesHdl(iMap);

        [Ell ,phi, g_thisMap, h_thisMap]=estimateGradHessEll(Nii_this,s,theta_this,timeVector,x+v,refNii,'timePoint_step',options.timePoint_step,'axesHdl',thisAxes,'do_gauss',do_gauss);

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
    
    % Add something to the Hessians for stability
    Hreg    = sum(H(:,:,:,1:size(Theta,4)),4)*(sqrt(eps('single')));
    H(:,:,:,1:size(Theta,4)) = bsxfun(@plus,H(:,:,:,1:size(Theta,4)),Hreg);

    Theta = Theta - spm_field(H,g+gr,[reg 1 1],regsc);  % Solve theta = theta - H\g
    Theta(:,:,:,1) = max(Theta(:,:,:,1),0);             % set to 0 negative values

    
end


















