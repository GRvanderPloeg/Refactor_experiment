function p = prox_L2DownConv(x,gamma,Lamb,d,y)
%function p = prox_L2DownConv(x,gamma,Lamb,d,y)
%
% This procedure computes the proximity operator of
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
%
%  REFERENCE
% ===========
% [1] Emmanuel Soubies and Michael Unser. "Computational Super-Sectioning for Single-Slice
%     Structured-Illumination Microscopy",  5-2 (2019), pp.~240--250.

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
yUp=zeros(sz);
yUp(sel{:})=y;

szPatch=sz./d;
for n=1:length(sz)
    sel2{n}=szPatch(n)*ones(1,d(n));
end
P= @(x) repmat(x,d);
H2H2t=Pt(abs(Lamb).^2,sel2,szPatch);   
                
% prox computation     
fftr=conj(Lamb).*fftn(yUp)+fftn(x/gamma);
p=real(ifftn(gamma*(fftr - conj(Lamb).*P(Pt(Lamb.*fftr,sel2,szPatch)./(prod(d)/gamma+H2H2t)))));         
end

function y=Pt(x,sel,szPatch)
y=zeros(szPatch);
tmp=mat2cell(x,sel{:});
for n=1:numel(tmp)
    y=y+tmp{n};
end
end
