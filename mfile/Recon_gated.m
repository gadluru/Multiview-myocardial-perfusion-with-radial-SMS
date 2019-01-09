function Recon_gated(Image_ref,Image_PD,Data,para)
%% self gating
para = auto_gating_SMS_yt(Image_ref,para);

%% motion estimation
for i=1:length(Data)
    para{i}.Recon.noi = 100;
    Data{i}.Motion = get_motion_SMS(Image_ref(:,:,:,i,:),1,10);
    Data{i}.Motion_bins = get_motion_SMS_bins(Image_ref(:,:,:,i,:),para{i},1,10);
end

%% pixel tracking STCR recon
for i=1:length(Data)
    para{i}.Recon.weight_tTV = para{i}.Recon.weight_tTV*3;
    [Image(:,:,:,i,:),para{i}] = STCR_pixel_tracking_bins(Data{i},para{i});
end

%% save image and parameters
disp('Saving image into Results...')
Image_PD = abs(crop_half_FOV(Image_PD));
Image = abs(crop_half_FOV(Image));
save([para{1}.dir.save_recon_img_mat_dir,para{1}.dir.save_recon_img_name],'Image','Image_PD','para','-v7.3');
disp('Reconstruction done');fprintf('\n')
