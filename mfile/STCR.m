function [Image,para] = STCR(Data,para)

disp('Iterative STCR reconstruction...');
t1 = tic;

ifGPU      = para.setting.ifGPU;
ifplot     = para.setting.plot;
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

scale_image = max(abs(Data.first_est(:)));
para.Recon.weight_sTV = para.weight_sTV*scale_image;
para.Recon.weight_tTV = para.weight_tTV*scale_image;

if isfield(Data,'phase_mod')
    Data.phase_mod_conj = conj(single(Data.phase_mod));
end
if isfield(Data,'sens_map')
    Data.sens_map_conj = conj(Data.sens_map);
end
if isfield(Data,'first_guess')
    Image = Data.first_guess;
else
    Image = single(Data.first_est);
end

if ifGPU
    Image = gpuArray(Image);
    Data.kSpace = gpuArray(Data.kSpace);
    Data.sens_map = gpuArray(Data.sens_map);
    Data.sens_map_conj = gpuArray(Data.sens_map_conj);
    beta_sqrd = gpuArray(beta_sqrd);
    if isfield(Data,'N')
        para.Recon.type = 'NUFFT less memory';
        for i=1:length(Data.N)
            Data.N(i).S = gpuArray(Data.N(i).S);
            Data.N(i).kb_density_comp = gpuArray(Data.N(i).kb_density_comp);
            Data.N(i).W = gpuArray(Data.N(i).W);
        end
    end
end

para.Cost = struct('fidelityNorm',[],'temporalNorm',[],'spatialNorm',[],'totalCost',[]);

fidelity = @(im) compute_fidelity_yt_new(im,Data,para);
sTV      = @(im) compute_sTV_yt(im,para.Recon.weight_sTV,beta_sqrd);
tTV      = @(im) compute_tTV_yt(im,para.Recon.weight_tTV,beta_sqrd);

for iter_no = 1:para.Recon.noi
    
    [update_term,fidelity_norm] = fidelity(Image);
    update_term = update_term + sTV(Image);
    update_term = update_term + tTV(Image);

    if iter_no > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
    end
    
    update_term_old = update_term; clear update_term

    para.Cost = Cost_STCR(fidelity_norm, Image, para, para.Cost); clear fidelity_update
    step_size = line_search(Image,update_term_old,Data,para);
    para.step_size(iter_no) = step_size;

    Image = Image + step_size * update_term_old;
    
    if ifplot ==1
        showImage(Image,para.Cost)
    end

    if iter_no > 1 && para.Recon.break
        if step_size<1e-4 || abs(para.Cost.totalCost(end) - para.Cost.totalCost(end-1))/para.Cost.totalCost(end-1) < 1e-4
            break
        end
    end

end

Image = squeeze(gather(Image));
if isfield(para.Recon,'PD_frames')
    para.PD_Cost = para.Cost;
    para = rmfield(para,'Cost');
end
disp('Reconstruction done.');toc(t1);fprintf('\n')