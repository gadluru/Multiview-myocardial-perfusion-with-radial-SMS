function [G,kSpace_cart] = GROG_seperate_SMS( kSpace_radial, kx, ky, phase_mod, para)
%  This is a GROG interpolation function for 2D MRI k-space data

%  Reference:
%  Nicole Seiberlich, et al (2008) Magnetic Resonance in Medicine
%  Self-Calibrating GRAPPA Operator Gridding for Radial and Spiral
%  Trajectories. 59:930-935.

%  Inputs:
%         kSpace_radial: up to 5 dimentions
%         kx, ky       : k-space positions

%  Copyright. Ye Tian, U of U DiBella Group
%  phye1988@gmail.com


% This version does not support multi-slices inputs.
nSMS = para.Recon.nSMS;
[sx,nor,nof,nc,~,NSlice] = size(kSpace_radial);
core_size = para.core_size;
over_sampling = para.over_sampling;

if size(phase_mod,5) == 3
    phase_mod = phase_mod(:,:,:,:,2);
    idx = phase_mod==0;
    phase_mod = round(imag(phase_mod));
    temp = [0,1,-1];
    N = sum(idx,2)/3;
    for i = 1:size(idx,3)
        phase_mod(1,idx(1,:,i),i) = repmat(temp,[1,N(i)]);
    end
else
    phase_mod(phase_mod==2) = -1;
end
if sum(vec(phase_mod)) == numel(phase_mod)
    phase_mod(phase_mod==1) = 0;
end

SMS(1,:,:,1,1,1) = phase_mod==0;
SMS(1,:,:,1,1,2) = phase_mod==1;
SMS(1,:,:,1,1,3) = phase_mod==-1;
G = cell(1,nSMS);
for j=1:nSMS
    for i=1:NSlice
        kSpace_radial_temp = kSpace_radial(:,:,:,:,i);
        kSpace_radial_temp = kSpace_radial_temp(repmat(SMS(:,:,:,:,:,j),[sx 1 1 nc]));
        kx_temp = kx(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
        ky_temp = ky(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
        kSpace_radial_temp = reshape(kSpace_radial_temp,[sx nor/nSMS nof nc]);
        kx_temp = reshape(kx_temp,[sx nor/nSMS nof]);
        ky_temp = reshape(ky_temp,[sx nor/nSMS nof]);
        
        G{j} = GROG.init(kSpace_radial_temp,kx_temp,ky_temp,over_sampling,core_size);
        kSpace_cart(:,:,:,:,:,:,j) = GROG.rad2cart(kSpace_radial_temp,G{j});
    end
end

sx_over = size(kSpace_cart,1);
kSpace_cart = reshape(kSpace_cart,[sx_over,sx_over,nof,nc,1,NSlice,nSMS]);

if size(kSpace_cart,1)~=sx*over_sampling
    kSpace_cart([1,end],:,:,:,:,:,:) = [];
    kSpace_cart(:,[1,end],:,:,:,:,:) = [];
end
