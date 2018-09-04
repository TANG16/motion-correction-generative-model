function Y = AffScale(T,M)
%Compute an affine mapping as AffMap in addition scales using specified
% transformation T
% FORMAT Y = AffMap(T,M)
% T = identity trasform (x) + velocity field (v)
% M affine (i.e. ratio between nii.mat of 2 images )
dm=size(T(:,:,:,1));
Y          = zeros([dm 3],'single');
Y(:,:,:,1) = M(1,1)*T(:,:,:,1) + M(1,2)*T(:,:,:,2) + M(1,3)*T(:,:,:,3) + M(1,4);
Y(:,:,:,2) = M(2,1)*T(:,:,:,1) + M(2,2)*T(:,:,:,2) + M(2,3)*T(:,:,:,3) + M(2,4);
Y(:,:,:,3) = M(3,1)*T(:,:,:,1) + M(3,2)*T(:,:,:,2) + M(3,3)*T(:,:,:,3) + M(3,4);

end

