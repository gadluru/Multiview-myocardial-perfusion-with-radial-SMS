function Reconstruction_multi_set_SMS_ungated(para)

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
    para{i}.weight_tTV = para{i}.weight_tTV/3;
    [Image_ref(:,:,:,i,:),para{i}] = STCR(Data{i},para{i});
    Data{i}.first_guess = Image_ref(:,:,:,i,:);
end

%% pixel tracking reconstruction
if para_temp.gated
    Recon_gated(Image_ref,Image_PD,Data,para);
else
    Recon_ungated(Image_ref,Image_PD,Data,para);
end

end