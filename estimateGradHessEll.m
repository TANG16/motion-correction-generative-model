function [Ell,varargout]  = estimateGradHessEll(Nii,s,parameterMaps,timeVector,T,refNii,varargin)
% the function estimates the gradient Hessian and conditional probability
% using all maps included in Nii (nifti files). in this model parameters
% are a and b from an exponential fit
%
% INPUTS: 
%	Nii: array of nifti object(i.e all echoes of a map)
%       do_gauss = boolean (1 fo gaussian model, 0 for Rice ) -> Rice model
%       need to be optimized
%   s: variance of the model
%   parameterMaps: parameters map each map parameter along the 4th dimension
%       (i.e. parameterMaps(:,:,:,1)=b; parameterMaps(:,:,:,2)=a 
%       (note here only parameters related to this Nii series, i.e only parameters of this map)
%   time: vector with time points should be one for each Nii
%   T: transforamtion needed to parameterMaps to warp to this time point
%       map
%   refNii: nifti file to use as reference for affine
%
% OUTPUTS:
%   Ell: conditional probability based on the model 
%   
% OPTIONAL OUTPUTS:
%   phi: is the deformation
%   g: [Nii_dim 2] gradient or first derivative
%   H: [Nii_dim 3] Hessian or second derivative
%
% MA: 29-08-2018


%% input/output 
% input
timePoint_step=2; % should be an optional input fro varagin
plotModel=1;% for plotting staff 
plotModel=0;% for plotting staff 
figureHdl=[]; %figure handle where to plot, if empty new figure
do_gauss=1; % should be an optional input to change default 
% output

if nargout>1 % because the first 2 output are necesary
    isGradHess=true; % calculate gradient/hessian
else
    isGradHess=false; % NOT calculate gradient/hessian
end

% initialization
Ell=0;
dm=Nii(1).dat.dim; % assuming all images of same size

%% main loop
    for iTP=1:timePoint_step:numel(Nii) % Loop over time points
    	if iTP==1
            Mi  = Nii(iTP).mat;
%             M1  = Nii(1).mat;
            phi = AffScale(T,refNii.mat\Mi);
            spm_diffeo('boundary',1);            % Neumann boundary condition

            modelParam_a  = spm_diffeo('samp',parameterMaps(:,:,:,2),phi); % alpha of this map
            modelParam_b  = spm_diffeo('samp',parameterMaps(:,:,:,1),  phi); % common beta 
            if isGradHess
                g   = zeros([dm 2],'single');
                H   = zeros([dm 3],'single');
            end
        end
        y    = single(Nii(iTP).dat(:,:,:,1));
        msk  = uint8(isfinite(y) & (y~=0) & isfinite(modelParam_a) & isfinite(modelParam_b));
        t    = timeVector(iTP);
        f    = exp(modelParam_a - modelParam_b*t);
        f1   = y.*f/s;  % needed for besseli input
        f2   = f.^2;    % needed for H update (in the definintion of x for hessian)
        g0   = zeros(dm,'single');
        h0   = zeros(dm,'single');
            
        % Visualisation
        if plotModel
            if ~isempty( figureHdl)
                figure(figureHdl) % could be improved without popup
                subplot(4,1,1);
            else
                if iTP==1, figure , end % needcd only once and overwrite the same over the loop
            end
            imagesc([y(:,:,ceil(end/2))' f(:,:,ceil(end/2))' (y(:,:,ceil(end/2))'-f(:,:,ceil(end/2))').*single(msk(:,:,ceil(end/2)))']);
            title('ModelFit: acquired Image - model - difference(acquired Image-model)')
            axis image xy off;
            drawnow
        end
                    
        
        % Conditions for using Rice model. I1/I0 fails for single precision
        % values over about 91.9.
        % following line to be considered for Rice model
%         if and(it>nits-2,~do_gauss), msk(msk & (f1<91) & (f2<91*(4*s))) = 2; end % here we define in msk were we can apply rice distribution

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % note from the original version of John:
        % msk == 1  -> where we apply Gaussian distr
        % msk == 2  -> where we apply Rice distr

        % GAUSSIAN MODEL
        % (applied where msk ==1, even for Rice case, where Rice is not working)
        % Compute likelihood (ll), gradient (g0) and approximate Hessian (h0)
        % for voxels using Gaussian model
        if do_gauss
            Ell = Ell - 0.5*sum(log(s) + log(2*pi) + (y(msk==1)-f(msk==1)).^2/s);
        else
            % here Rice calculations -> to be checked
            I0       = besseli(0, f1(msk==2)); % These steps are slow
            I1       = besseli(1, f1(msk==2));
            Ell       = Ell + sum(log(y(msk==2)) - 0.5*y(msk==2).^2/s - 0.5*f(msk==2).^2/s + log(I0)); % Log-likelihood  -> please double check if it is not added twice from gauss model
        end
        
        if nargout==1
            % if only Ell is need go to the next iteration for speeding
            continue
        end
        
        g0(msk==1) = f2(msk==1)/s - f1(msk==1);
        h0(msk==1) = f2(msk==1)/s;
        % here additional optimization from last John's e-mail:
        % use the full Hessian for those points where
        %     2*exp(a-b*t) - y > 0
        % but only the approximate Hessian if not 
        % -> this should probably give a more stable optimisation than the prev version.
        %
        % from matlab the full Hessian has follwing form:
        % H(1,1)=f2/s +((f*(f-y))/s)
        % so it should be -> please uncomment the following
        isFullHessian=and(2*f - y > 0,msk==1);
        h0(isFullHessian)=2.*f2(isFullHessian)/s-f1(isFullHessian);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if ~do_gauss
            % RICIAN MODEL
            % Compute likelihood (ll), gradient (g0) and approximate Hessian (h0)
            % for voxels using Rician model
%             I0       = besseli(0, f1(msk==2)); % These steps are slow
%             I1       = besseli(1, f1(msk==2));
%             Ell       = Ell + sum(log(y(msk==2)) - 0.5*y(msk==2).^2/s - 0.5*f(msk==2).^2/s + log(I0)); % Log-likelihood  -> please double check if it is not added twice from gauss model

            % I_1(x)/I_0(x) could be roughly approximated by x./(alpha+sqrt(beta^2+x.^2)), where alpha=0.42
            % and beta=1.43.  See (eg) Hornik & Grun. "Amos-type bounds for modified Bessel function ratios".
            % J Math Anal Appl. 2013 Dec 1; 408(1): 91???101.  Note that brackets are needed around (I1./I0)
            % to prevent overflow of f1.*I1.
            g0(msk==2) = f2(msk==2)/s - f1(msk==2).*(I1./I0);

            % Use a Gaussian approximation of the Rice distribution for computing a stable Hessian
            x          = -f2(msk==2)/(2*s);
            Laguerre   = exp(x/2).*((1-x).*besseli(0,-x/2)-x.*besseli(1,-x/2));
            v          = 2*s + f2(msk==2) - pi*s/2*Laguerre.^2; % Variance of Rician distribution
            h0(msk==2) = f2(msk==2)./v;
           %Ey         = sqrt(pi*s/2).*Laguerre;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INCREMENTS
            % Increment gradients 
            g(:,:,:,1) = g(:,:,:,1) + g0*(-t);    % df/d_beta
            g(:,:,:,2) = g(:,:,:,2) + g0;         % df/d_alpha

            % Increment Hessians 
            H(:,:,:,1) = H(:,:,:,1) + h0*(t^2);     % d2f/d_betas^2
            H(:,:,:,2) = H(:,:,:,2) + h0;           % d2f/d_alpha^2
            H(:,:,:,3) = H(:,:,:,3) + h0*(-t);      % d2f/d_alpha_d_beta
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('.');
    end % end loop over all time series
    
    if nargout>1
        varargout{1}=phi;
        varargout{2}=g;
        varargout{3}=H;
    end
    
end % end of this fcn







