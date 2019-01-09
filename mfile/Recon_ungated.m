function [Image_sys,Image_dia,para] = Recon_ungated(Image_ref,Image_PD,Data,para)

cardiac_signal = self_gating_image_space(crop_half_FOV(Image_ref));
respiration_signal = get_respiration_bins(crop_half_FOV(Image_ref));
nof = min(min(sum(cardiac_signal==1,1)),min(sum(cardiac_signal==2,1)));
Image_sys = zeros(size(Image_ref),'single');
Image_sys(:,:,nof+1:end,:,:) = [];
Image_dia = Image_sys;

for i=1:size(Image_ref,4)
    sys = cardiac_signal(:,i) == 1;
    respiration_signal_temp = respiration_signal(sys);
    Data_sys.filter = Data{i}.filter;
    Data_sys.kSpace = Data{i}.kSpace(:,:,sys,:,:,:,:);
    Data_sys.mask = Data{i}.mask(:,:,sys,:,:,:,:);
    Data_sys.sens_map = Data{i}.sens_map;
    Data_sys.first_est = Data{i}.first_est(:,:,sys,:,:);
    Data_sys.first_guess = Data{i}.first_guess(:,:,sys,:,:);

    if isfield(Data{i},'SMS')
        Data_sys.SMS = Data{i}.SMS;
    end
    
    para_sys = para{i};
    para_sys.Recon.bins(1,:) = respiration_signal_temp == 1;
    para_sys.Recon.bins(2,:) = respiration_signal_temp == 2;
    para_sys.Recon.bins(3,:) = respiration_signal_temp == 3;
    para_sys.Recon.bins(4,:) = respiration_signal_temp == 4;
    para_sys.Recon.bins(:,1) = true;
    para_sys.Recon.bins(:,end) = true;
    
    idx = sum(para_sys.Recon.bins,2)<=2;
    para_sys.Recon.bins(idx,:) = [];
        
    para{i}.respiration_signal = respiration_signal;
    para{i}.cardiac_signal = cardiac_signal;
    
    Data_sys.Motion = get_motion_SMS(Image_ref(:,:,sys,i,:),1,20);
    Data_sys.Motion_bins = get_motion_SMS_bins(Image_ref(:,:,sys,i,:),para_sys,1,20);
    
    para_sys.Recon.weight_tTV = para{i}.Recon.weight_tTV*3;
    para_sys.Recon.noi = 70;
    Image_sys_temp = STCR_pixel_tracking_bins(Data_sys,para_sys);
    Image_sys(:,:,:,i,:) = Image_sys_temp(:,:,1:nof,:,:);
end

for i=1:size(Image_ref,4)
    dia = cardiac_signal(:,i) == 2;
    respiration_signal_temp = respiration_signal(dia);
    Data_dia.filter = Data{i}.filter;
    Data_dia.kSpace = Data{i}.kSpace(:,:,dia,:,:,:,:);
    Data_dia.mask = Data{i}.mask(:,:,dia,:,:,:,:);
    Data_dia.sens_map = Data{i}.sens_map;
    Data_dia.first_est = Data{i}.first_est(:,:,dia,:,:);
    Data_dia.first_guess = Data{i}.first_guess(:,:,dia,:,:);

    if isfield(Data{i},'SMS')
        Data_dia.SMS = Data{i}.SMS;
    end
    
    para_dia = para{i};
    para_dia.Recon.bins(1,:) = respiration_signal_temp == 1;
    para_dia.Recon.bins(2,:) = respiration_signal_temp == 2;
    para_dia.Recon.bins(3,:) = respiration_signal_temp == 3;
    para_dia.Recon.bins(4,:) = respiration_signal_temp == 4;
    para_dia.Recon.bins(:,1) = true;
    para_dia.Recon.bins(:,end) = true;
    
    idx = sum(para_dia.Recon.bins,2)<=2;
    para_dia.Recon.bins(idx,:) = [];
        
    para{i}.respiration_signal = respiration_signal;
    para{i}.cardiac_signal = cardiac_signal;
    
    Data_dia.Motion = get_motion_SMS(Image_ref(:,:,dia,i,:),1,20);
    Data_dia.Motion_bins = get_motion_SMS_bins(Image_ref(:,:,dia,i,:),para_dia,1,20);
    
    para_dia.Recon.weight_tTV = para{i}.Recon.weight_tTV*3;
    para_dia.Recon.noi = 70;
    Image_dia_temp = STCR_pixel_tracking_bins(Data_dia,para_dia);
    Image_dia(:,:,:,i,:) = Image_dia_temp(:,:,1:nof,:,:);
end

disp('Saving image into Results...')
Image_PD = abs(crop_half_FOV(Image_PD));
Image_sys = abs(crop_half_FOV(Image_sys));
Image_dia = abs(crop_half_FOV(Image_dia));
save([para{1}.dir.save_recon_img_mat_dir,para{1}.dir.save_recon_img_name],'Image_sys','Image_dia','Image_PD','para','-v7.3');
disp('Reconstruction done');fprintf('\n')
