function Y = AffMap(d,M)
% Compute an affine mapping
% FORMAT Y = AffMap(d,M)

if iscell(d)
    [x1,x2,x3] = ndgrid(single(d{1}),single(d{2}),single(d{3}));
    dm = [size(x1,1) size(x1,2) size(x1,3)];
else
    dm = [d 1 1 1];
    dm = dm(1:3);
    [x1,x2,x3] = ndgrid(single(1:dm(1)),single(1:dm(2)),single(1:dm(3)));
end

if nargin<2,
    M = single(eye(4));
else
    M = single(M);
end
Y          = zeros([dm 3],'single');
Y(:,:,:,1) = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
Y(:,:,:,2) = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
Y(:,:,:,3) = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);

