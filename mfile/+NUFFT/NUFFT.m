function kSpace_radial = NUFFT(image,N)

sx = N.siz(1);
sx_over = N.sx_over;
core_size = N.core_size;
nof = N.siz(3);
nc = size(image,4);
nSMS = size(image,5);
throw_away_x = (sx_over+core_size(1)-sx)/2;
throw_away_y = (sx_over+core_size(2)-sx)/2;
throw_away_x = round(throw_away_x);
throw_away_y = round(throw_away_y);
%image(1:2:end,:,:,:) = -image(1:2:end,:,:,:);
%image(:,1:2:end,:,:) = -image(:,1:2:end,:,:);
if mod(core_size(1),2) == 0
    image_over = zeros(sx_over+core_size(1),sx_over+core_size(2),nof,nc,nSMS);
else
    image_over = zeros(sx_over+core_size(1)+1,sx_over+core_size(2)+1,nof,nc,nSMS);
end
image = bsxfun(@times,image,N.kb_density_comp);
image_over(throw_away_x+1:throw_away_x+sx,throw_away_y+1:throw_away_y+sx,:,:,:) = image;

%image_over(1:2:end,:,:,:) = - image_over(1:2:end,:,:,:);
%image_over(:,1:2:end,:,:) = - image_over(:,1:2:end,:,:);

image_over = fftshift(image_over,1);
image_over = fftshift(image_over,2);

kSpace_cart = fft2(image_over);

kSpace_cart = fftshift(kSpace_cart,1);
kSpace_cart = fftshift(kSpace_cart,2);
%kSpace_cart(1:2:end,:,:,:) = - kSpace_cart(1:2:end,:,:,:);
%kSpace_cart(:,1:2:end,:,:) = - kSpace_cart(:,1:2:end,:,:);

kSpace_radial = NUFFT.cart2rad(kSpace_cart,N);