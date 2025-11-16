function p = fun_L2DownConv(x,gamma,Lamb,d,y)
%function p = fun_L2DownConv(x,gamma,Lamb,d,y)
%
% This procedure evaluates the following functional
%
% f(x) = gamma/2 ||S_d H x -y ||^2_2
%
% where S_d is a downsampling operator with factor d>0 and H is a convolution
% operator that is described by the Fourier kernel Lamb 
% [i.e. Hx = ifftn(Lamb.*fftn(x)) ]
%
%  INPUTS
% ========
%  x      - ND array
%  gamma  - positive scalar 
%  Lamb   - ND array with the same size as 'x' (Fourier kernel of the
%           convolution operator H)
%  d      - vector of positive integers (downsampling factor in each dimension of x)
%  y      - ND array with size : size(x)./d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (23-10-2019)
% Author  : Emmanuel Soubies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check input
if ~isscalar(gamma) ||  gamma<=0
    error('''gamma'' must be a positive scalar')
end
if size(x)~=size(Lamb)
    error('''Lamb'' must be an ND array with the same size as ''x''')
end
if size(y)~=size(x)./d
    error('''y'' must be an ND array with size : size(x)./d')
end
%-----%

% Precomputations
sz=size(x);
sel=cell(length(sz),1);
for ii=1:length(size(d))
    sel{ii}=1:d(ii):sz(ii);
end

% evaluate the function
Hx=real(ifftn(Lamb.*fftn(x)));
p = gamma/2*norm(Hx(sel{:}) - y,'fro')^2;

end


