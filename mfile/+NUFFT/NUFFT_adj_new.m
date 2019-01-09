function image = NUFFT_adj_new(kSpace_radial,N)

kSpace_radial = bsxfun(@times,N.W,kSpace_radial); % density compensation
%filter
kSpace_cart = NUFFT.rad2cart_new(kSpace_radial,N);

image = ifft2(kSpace_cart);

sx = N.siz(1);

if N.sx_over >= N.siz(1)
    image = image(1:sx,1:sx,:,:,:,:);
else
    image = fftshift2(image);
end

image = bsxfun(@times,image,N.kb_density_comp);