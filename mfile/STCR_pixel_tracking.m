function [Image,para] = STCR_pixel_tracking(Data,para)

disp('Pixel-Tracking STCR reconstruction...');
t1 = tic;

ifGPU      = para.setting.ifGPU;
ifplot     = para.setting.plot;
weight_sTV = para.Recon.weight_sTV;
weight_tTV = para.Recon.weight_tTV;
beta_sqrd  = eps('single');

if isfield(para.Recon,'PD_frames')
    PD_frames = para.Recon.PD_frames;
    if PD_frames(end) == 0
        PD_frames = 1:sum(PD_frames);
    end
    Data.first_est = Data.first_est(:,:,PD_frames,:,:);
    Data.kSpace = Data.kSpace(:,:,PD_frames,:,:,:,:);
    if isfield(Data,'mask')
        Data.mask = Data.mask(:,:,PD_frames,:,:,:,:);
    end
    if isfield(Data,'phase_mod')
        Data.phase_mod = Data.phase_mod(:,:,PD_frames,:,:,:,:);
    end
    if isfield(Data,'N')
        Data.N.S = Data.N.S(1:Data.N.sx_over.^2*PD_frames(end),1:Data.N.siz(1)*Data.N.siz(2)*PD_frames(end));
        Data.N.siz(3) = size(Data.first_est,3);
    end
end

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end
if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end

if isfield(Data,'first_guess')
    Image = Data.first_guess;
    noi_start = length(para.step_size)+1;
    para.Recon.noi = para.Recon.noi + noi_start;
%    para.Cost = para.PD_Cost;
else
    Image = Data.first_est;
    noi_start = 1;
    para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);
end

if ifGPU
    Image = gpuArray(Image);
    Data.first_est = gpuArray(Data.first_est);
    Data.sens_map = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    beta_sqrd = gpuArray(beta_sqrd);
    if isfield(Data,'N')
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).kb_density_comp = gpuArray(Data.N(i).kb_density_comp);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
end

%para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
sTV      = @(im) compute_sTV_yt(im,para.Recon.weight_sTV,beta_sqrd);
tTV      = @(im) compute_tTV_pixel(im,para.Recon.weight_tTV,beta_sqrd,Data.Motion);

for iter_no = noi_start:para.Recon.noi

    [update_term,fidelity_norm] = fidelity(Image);
    update_term = update_term + sTV(Image);
    update_term = update_term + tTV(Image);
    
    if iter_no > noi_start
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    
    update_term_old = update_term; clear update_term
    
    para.Cost = Cost_STCR_pixel(fidelity_norm, Image, Data.Motion, weight_sTV, weight_tTV, para.Cost); clear fidelity_update
    step_size = line_search_pixel(Image,update_term_old,Data,para);
    para.step_size(iter_no) = step_size;

    Image   = Image + step_size * update_term_old;
    
    if ifplot ==1
        showImage(Image,para.Cost)
    end

    if iter_no > 1 && para.Recon.break
        if step_size<1e-4 || abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end

end

Image = squeeze(gather(abs(Image)));
% para.PD_Cost = para.Cost;
% para.PD_step_size = para.step_size;
para = rmfield(para,'Cost');
para.CPUtime.PD = toc(t1);

disp('Reconstruction done.');toc(t1);fprintf('\n')