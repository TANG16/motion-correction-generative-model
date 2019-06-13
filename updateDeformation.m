function [Ell, Ev,V] = updateDeformation (Nii,mapIdx,s,Theta,T,V,alpha,settings_v,varargin)
% the fcn will return the model fit parameters
%
% SYNOPSIS: []=updateDeformation (varagin)
%
% INPUT:
%	Nii: array of nifti object(i.e all echoes of a map)
%       do_gauss = boolean (1 fo gaussian model, 0 for Rice ) -> Rice model
%       need to be optimized
%   mapIdx: grouping index for contrasts   
%   s: variance of the model
%   Theta: parameters map each map parameter along the 4th dimension
%       (i.e. parameterMaps(:,:,:,1)=b; parameterMaps(:,:,:,2)=a 
%       (note here only parameters related to this Nii series, i.e only parameters of this map)
%   T: transformation needed to parameterMaps to warp to this time point
%       map
%   V: velocity fields
%   alpha: used to weight the update
%   settings_v: settings for the deformation
% 
% OPTIONAL INPUTS: 
% look 'estimateGradHessEll' fcn
% thisAxes axes for each map
%
% OUTPUT the model fir parameters
%   Ell: probability of the likelihood
%   Ev: conditional prob of the deformation
%   updated V, deformation fields 
%
% REMARKS
% Deformations
%=======================================================================
% Minimises -log p(F,V,mu).  Note that some constants are
% ignored as they don't feature in the objective function
%-----------------------------------------------------------------------
%
% SEE ALSO Affscale,estimateGradHessEll
%
% EXAMPLE
%         [Ell, Ev,V] = updateDeformation (Nii,mapIdx,s,Theta,T,V,alpha,settings_v,'timePoint_step',timePoint_step)  ;      
%
% created with MATLAB ver.: 9.1.0.441655 (R2016b) on Microsoft Windows 7 Professional  Version 6.1 (Build 7601: Service Pack 1)
%
% created by: mazzarito
% DATE: 14-Sep-2018
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
    Ell=0;  % -log p(F|mu,V) = \sum_i log p(f_i | mu, v_i)
    Ev=0;   % -log p(V)      = \sum_i log p(v_i)         

    dm=Nii(1).dat.dim; % assuming all images of same size
    np=numel(T)+1;

    %% intialize
    
    % Construct an identity transform
    [x1,x2,x3] = ndgrid(single(1:dm(1)),single(1:dm(2)), single(1:dm(3)));
    x = cat(4,x1,x2,x3);

    g   = zeros([dm, np],'single');         % Gradients
    H   = zeros([dm,(np+1)*np/2],'single'); % Hessians

    
        for iMap=1:numel(T) % loop over modalities

           %% -- compute g_map and H_map
            v=V{iMap};
            Nii_this=Nii(mapIdx==iMap);
            timeVector=T{iMap};
            theta_this=cat(4,Theta(:,:,:,1),Theta(:,:,:,iMap+1));       

            [Ell ,phi, g_thisMap, H_thisMap]=estimateGradHessEll(Nii_this,s,theta_this,timeVector,x+v,Nii(1),'timePoint_step',options.timePoint_step,'do_gauss',do_gauss); % add input for not plotting

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
    
    
    
end


















