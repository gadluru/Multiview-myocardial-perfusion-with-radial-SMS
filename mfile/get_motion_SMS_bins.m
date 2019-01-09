function Motion_bins = get_motion_SMS_bins(Image,para,step,noi)

bins = para.Recon.bins;
nb = size(bins,1);

for i=1:nb
    bin_temp = bins(i,:);
    Motion_bins{i} = get_motion_SMS(Image(:,:,bin_temp,:),step,noi);
end