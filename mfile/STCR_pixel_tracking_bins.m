function [Image,para] = STCR_pixel_tracking_bins(Data,para)
% [Image,para] = iteration_recon(first_est,kSpace,sens_map,phase_mod,para,varargin)
disp('Pixel-tracking STCR with binning reconstruction...');
t1 = tic;

ifplot     = para.setting.plot;
ifGPU      = para.setting.ifGPU;
weight_tTV = para.Recon.weight_tTV;
weight_sTV = para.Recon.weight_sTV;
beta_sqrd  = eps('single');


if isfield(Data,'first_guess')
    new_img_x = Data.first_guess;
    noi_start = length(para.step_size)+1;
    para.Recon.noi = para.Recon.noi + noi_start;
    para.Cost = para.Cost;
else
    new_img_x = Data.first_est;
    noi_start = 1;
    para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end

if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if ifGPU
    Data.kSpace        = gpuArray(Data.kSpace);
    new_img_x          = gpuArray(new_img_x);
    Data.sens_map      = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).kb_density_comp = gpuArray(Data.N(i).kb_density_comp);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end

    gpuInfo = gpuDevice;
    gpuSize = gpuInfo.AvailableMemory;
    imSize  = numel(new_img_x)*8;
    if imSize*para.Recon.no_comp > gpuSize*0.3
        para.Recon.type = [para.Recon.type,' less memory'];
    end
end

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
spatial  = @(im) compute_sTV_yt(im,weight_sTV,beta_sqrd);
ptv_all  = @(im) compute_tTV_pixel(im,weight_tTV*0.5,beta_sqrd,Data.Motion);
ptv_bin  = @(im) compute_tTV_pixel_bins(im,weight_tTV*0.5,beta_sqrd,Data.Motion_bins,para.Recon.bins);

for iter_no = noi_start:para.Recon.noi

%%%%% fidelity term/temporal/spatial TV

    [update_term,fidelity_norm] = fidelity(new_img_x);
    update_term = update_term + ptv_all(new_img_x);
    update_term = update_term + ptv_bin(new_img_x);
    update_term = update_term + spatial(new_img_x);
    
%%%%% conjugate gradient
    tic;
    if iter_no > noi_start
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    update_term_old = update_term; clear update_term
    
%%%%% line search    

    para.Cost = Cost_STCR_pixel_bins(fidelity_norm, new_img_x, Data, para, para.Cost); clear fidelity_update
    step_size = line_search_pixel_bins(new_img_x,update_term_old,Data,para);
    para.step_size(iter_no) = step_size;

    new_img_x = new_img_x + step_size .* update_term_old;
    para.CPUtime.update(iter_no) = toc;

%%%%% plot&save part 

    if ifplot == 1
        showImage(new_img_x,para.Cost)
    end
    
    % break when step size too small or cost not changing too much
    if iter_no > 1 && para.Recon.break
        if step_size<1e-4 || abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end

end

Image = squeeze(gather(new_img_x));
toc(t1);fprintf('\n')
