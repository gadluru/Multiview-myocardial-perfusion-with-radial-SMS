function para = load_default_para(option)

para.setting.PCA = 1;
para.NumberofPCAComponent = 8;

if strfind(option,'SMS')
    para.phase_mod = 1;
else
    para.phase_mod = 0;
end
if strfind(option,'golden')
    para.angle_mod = 1;
end
if strfind(option,'cSMS')
    para.cSMS = 1;
else
    para.cSMS = 0;
end
if contains(option,'GROG')
    para.Recon.interp_method = 'GROG';
    para.over_sampling = 1;
    para.core_size = [1,1];
end
if contains(option,'NUFFT')
    para.Recon.interp_method = 'NUFFT';
    para.over_sampling = 1.5;
    para.core_size = [6,6];
end
if contains(option,'GPU')
    para.setting.ifGPU = 1;
else
    para.setting.ifGPU = 0;
end

para.setting.debug = 1;
para.setting.plot = 1;

para.weight_sTV = 0.001;
para.weight_tTV = 0.02;

para.step_size = 2;
para.Recon.crop_half_FOV = 1;
para.image_orintation = 1;

para.dir.save_recon_img_mat_dir = strcat(pwd,'/ReconData/mat_files/');