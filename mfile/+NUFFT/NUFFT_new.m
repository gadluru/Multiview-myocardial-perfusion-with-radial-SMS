function kSpace_radial = NUFFT_new(image,N)

sx_over = N.sx_over;

image = bsxfun(@times,image,N.kb_density_comp);

if N.sx_over < N.siz(1)
    image = fftshift2(image);
end
kSpace_cart = fft2(image,sx_over,sx_over);

kSpace_radial = NUFFT.cart2rad_new(kSpace_cart,N);