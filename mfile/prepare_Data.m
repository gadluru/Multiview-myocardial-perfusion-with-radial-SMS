function [Data,para] = prepare_Data(kSpace_radial,para)
% [kSpace,phase_mod,para] = prepare_kSpace(kSpace_all,theta_all,phase_mod_all,RayPosition,para)

t1 = tic;
[sx,~,no_comp,~,ns] = size(kSpace_radial);
nSMS                = para.Recon.nSMS;
kCenter             = para.kSpace_center;
interp_method       = para.Recon.interp_method;
nor                 = para.Recon.nor;

theta     = get_angle_mod(para);% radial sampling angle
phase_mod = single(get_phase_mod(para));% if SMS
nof = length(theta)/nor;
theta = reshape(theta,[1,nor,nof]);
kSpace_radial = reshape(kSpace_radial,[sx,nor,nof,no_comp,1,ns]);
phase_mod = reshape(phase_mod,[1,nor,nof,1,nSMS]);

%%%%% pre-interpolation or NUFFT
disp('Pre-interpolate onto Cartesian space...')

[kx, ky] = get_k_coor(sx,theta,0,kCenter);

switch interp_method

    case 'GROG'

        [Data.G,Data.kSpace] = GROG.GROG_seperate_SMS(squeeze(kSpace_radial),kx,ky,phase_mod, para);

        if nSMS == 1            
            para.Recon.type = '2D';
        else
            para.Recon.type = 'seperate SMS';
        end

    case 'NUFFT'
        Data.kSpace = kSpace_radial;
        Data.N = NUFFT.init_new(kx,ky,para.over_sampling,para.core_size);
        im = NUFFT.NUFFT_adj_new(kSpace_radial.*conj(phase_mod),Data.N);
        if nSMS == 1
            Data.sens_map = get_sens_map(im,'2D');
        else
            Data.sens_map = get_sens_map(im,'SMS');
            Data.phase_mod = phase_mod;
        end
        Data.first_est = bsxfun(@times,im,conj(Data.sens_map));
        Data.first_est = sum(Data.first_est,4);

        para.Recon.type = 'NUFFT';

        scale_image = max(abs(Data.first_est(:)));
        para.Recon.weight_tTV = scale_image*para.weight_tTV;
        para.Recon.weight_sTV = scale_image*para.weight_sTV;
        toc;

        return
        
end

% orintation
if para.image_orintation == 0
    para.image_orintation = orintation_detection(abs(fftshift(ifft2(sum(sum(sum(Data.kSpace,3),4),7)))));
    Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
else
    Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
end

Data.kSpace = fftshift2(Data.kSpace);
Data.mask = logical(abs(Data.kSpace(:,:,:,1,:,:,:)));

switch para.Recon.type
    case 'seperate SMS'

        Data.SMS = exp(1i*(0:nSMS-1)*2*pi/nSMS.*(0:nSMS-1).');
        Data.SMS = permute(Data.SMS,[3,4,5,6,1,7,2]);
        kSpace_sms = sum(Data.kSpace.*conj(Data.SMS),7);

        Data.kSpace = ifft2(Data.kSpace);
        Data.kSpace = fftshift2(Data.kSpace);
        Data.kSpace = fft2(Data.kSpace);
        Data.kSpace = Data.kSpace .* Data.mask;

        im = ifft2(kSpace_sms);
        im = fftshift2(im);
        kSpace_sms = fft2(im).*logical(sum(Data.mask,7));
        sx_over = size(im,1);
        cut = (sx_over - para.Recon.image_size(1))/2;
        para.Recon.kSpace_size = [sx_over,sx_over];
        im([1:cut,end-cut+1:end],:,:,:,:) = [];
        im(:,[1:cut,end-cut+1:end],:,:,:) = [];
        Data.sens_map = get_sens_map(im,'SMS');
        Data.filter = ramp_filter_for_pre_interp(para);
        
        im = ifft2(kSpace_sms.*Data.filter);
        im([1:cut,end-cut+1:end],:,:,:,:) = [];
        im(:,[1:cut,end-cut+1:end],:,:,:) = [];
        Data.first_est = single(sum(bsxfun(@times, im, conj(Data.sens_map)),4));

    case '2D'
        im = ifft2(Data.kSpace);
        im = fftshift2(im);
        Data.kSpace = fft2(im).*Data.mask;
        sx_over = size(im,1);
        cut = (sx_over - para.Recon.image_size(1))/2;
        para.Recon.kSpace_size = [sx_over,sx_over];
        im([1:cut,end-cut+1:end],:,:,:,:) = [];
        im(:,[1:cut,end-cut+1:end],:,:,:) = [];
        Data.sens_map = get_sens_map(im,'2D');
        Data.filter = ramp_filter_for_pre_interp(para);
        
        im = ifft2(Data.kSpace.*Data.filter);
        im([1:cut,end-cut+1:end],:,:,:,:) = [];
        im(:,[1:cut,end-cut+1:end],:,:,:) = [];
        Data.first_est = single(sum(bsxfun(@times, im, conj(Data.sens_map)),4));
end

scale_image = max(abs(Data.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;

para.CPUtime.prepare_kSpace = toc(t1);toc(t1);fprintf('\n');

