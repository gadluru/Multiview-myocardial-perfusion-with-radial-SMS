addpath mfile/
clear para

Data_files = dir('RawData/*kSpace.mat');

para = load_default_para('2D SMS golden NUFFT');
para.dir.load_kSpace_name = Data_files.name;
para.dir.load_kSpace_dir = [Data_files.folder,'/'];

para.setting.plot =01;
para.setting.ifGPU = 1;
para.gated = 1;

para.Recon.noi = 30;
para.Recon.break = 0;
para.weight_sTV = 0.003;
para.weight_tTV = 0.02;
para.Recon.nor = 30;
para.image_orintation = [7,3,3];
para.kSpace_center = 145;

para.slice_pick = 1;

[image,para] = Reconstruction_multi_set_SMS_ungated(para);