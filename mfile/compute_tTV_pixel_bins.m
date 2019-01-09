function tTV_update = compute_tTV_pixel_bins(Image,weight,beta_square,Motion,bins)
tTV_update = zeros(size(Image),'like',Image);

for i=1:size(bins,1)
    bin_temp = bins(i,:);
    Image_temp = Image(:,:,bin_temp,:,:);
    Motion_temp = Motion{i};
    tTV_update(:,:,bin_temp,:,:) = compute_tTV_pixel(Image_temp,weight,beta_square,Motion_temp);
    %tTV_update(:,:,bin_temp,:,:) = compute_tTV_yt(Image_temp,weight,beta_square);
end