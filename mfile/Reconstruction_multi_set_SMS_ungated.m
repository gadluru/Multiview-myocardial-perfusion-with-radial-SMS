function [Image,para] = Reconstruction_multi_set_SMS_ungated(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');

%% load parameters from para
para_temp = prepare_para(para);
clear para

%% load k-space data
[kSpace_all,para_temp] = PCA_kSpace_unstreaking_yt(para_temp);

%% prepare Data
for i=1:size(kSpace_all,5)
    para{i} = para_temp;
    if size(para_temp.phase_mod,2) > 1
        para{i}.phase_mod = para{i}.phase_mod(:,i);
        para{i}.angle_mod = para{i}.angle_mod(:,i);
    end
    if length(para_temp.image_orintation)>1
        para{i}.image_orintation = para{i}.image_orintation(i);
    end
    [Data{i},para{i}] = prepare_Data(kSpace_all(:,:,:,:,i),para{i});
end

%% reconstructe proton density images
for i=1:length(Data)
    [Image_PD_ref(:,:,:,i,:),para{i}] = STCR(Data{i},para{i});
    para{i}.Recon.noi = 70;
    Data{i}.Motion = get_motion_SMS(Image_PD_ref(:,:,:,i,:),1,10);
    Data{i}.first_guess = Image_PD_ref(:,:,:,i,:);
    [Image_PD(:,:,:,i,:),para{i}] = STCR_pixel_tracking(Data{i},para{i});
    [Data{i},para{i}] = remove_PD(Data{i},para{i});
end

%% reconstruct preliminary reference images
for i=1:length(Data)
    para{i}.Recon.noi = para_temp.Recon.noi;
    para{i}.step_size = para_temp.step_size;
    para{i}.weight_tTV = para{i}.weight_tTV/3;keyboard
    [Image_ref(:,:,:,i,:),para{i}] = STCR(Data{i},para{i});
    Data{i}.first_guess = Image_ref(:,:,:,i,:);
end

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

Image = cat(3,Image_PD,Image);

if para_temp.Recon.crop_half_FOV == 1
    Image = abs(crop_half_FOV(Image));
else
    Image = abs(Image);
end

disp('Saving image into Results...')

Image = gather(Image);
save([para_temp.dir.save_recon_img_mat_dir,para_temp.dir.save_recon_img_name],'Image','para','-v7.3');

disp('Reconstruction done');fprintf('\n')

end