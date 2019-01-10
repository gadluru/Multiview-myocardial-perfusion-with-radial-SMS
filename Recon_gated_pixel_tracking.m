addpath mfile/
clear para

% This will load the default reconstruction parameters for 2D radial SMS
% with GROG pre-iteration gridding. Another option is to use 'NUFFT'
% instead of 'GROG', which would use NUFFT during the iterations. The NUFFT
% codes were modified from Jeffery Fessler's oringinal codes at
% https://web.eecs.umich.edu/~fessler/
% to be compatible with GPU.
para = load_default_para('2D SMS golden GROG');
para.dir.load_kSpace_name = 'meas_MID00121_FID67424_UCAIR_SMS2_SA_LA_30_rays_10ml_stress_kSpace.mat';
para.dir.load_kSpace_dir = 'RawData/';

% This option, when = 1, will plot cost function and display updated image 
% at every iteration.
para.setting.plot = 0;

% If you have a MATLAB compatible GPU and enough GPU memory, set this to 1
% would run the iterative part of the code on GPU.
para.setting.ifGPU = 1;

para.gated = 1;

para.Recon.noi = 30;
para.Recon.break = 1;
para.weight_sTV = 0.003;
para.weight_tTV = 0.02;
para.Recon.nor = 30;
para.image_orintation = [7,3,3];
para.kSpace_center = 145;

% There are 3 SMS sets in the date. This option can be set to reconstructed
% any sets. 
para.slice_pick = 2;

% The final reconstructed image will be saved under
% /ReconData/mat_files/MID***/
% By loading in the reconstrtucted image and use the function
% 'show_yt_pause(reshape_MSMS(Image))' the images can be displayed.
Reconstruction_multi_set_SMS_ungated(para);
