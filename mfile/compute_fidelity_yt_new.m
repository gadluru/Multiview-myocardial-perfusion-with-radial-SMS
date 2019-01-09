function [fidelity_update,fidelity_norm] = compute_fidelity_yt_new(image,Data,para)

switch para.Recon.type
    case {'2D','2D less memory'}

        fidelity_update = bsxfun(@times,image,Data.sens_map);
        fidelity_update = fft2(fidelity_update,para.Recon.kSpace_size(1),para.Recon.kSpace_size(2));
        if isfield(Data,'phase_mod')
            fidelity_update = bsxfun(@times,fidelity_update,Data.phase_mod);
            fidelity_update = sum(fidelity_update,5);
            fidelity_update = bsxfun(@times,fidelity_update,Data.phase_mod_conj);
        else
            fidelity_update = bsxfun(@times,fidelity_update,Data.mask);
        end

        fidelity_update = Data.kSpace - fidelity_update;
        
        fidelity_norm = fidelity_update/para.Recon.kSpace_size(1)/8;
                
        if isfield(Data,'filter')
            fidelity_update = fidelity_update.*Data.filter; % use filter to accelerate converge
        end
        
        fidelity_update = ifft2(fidelity_update);
        fidelity_update(para.Recon.image_size(1)+1:end,:,:,:,:) = [];
        fidelity_update(:,para.Recon.image_size(2)+1:end,:,:,:) = [];

        fidelity_update = bsxfun(@times,fidelity_update,Data.sens_map_conj);
        fidelity_update = sum(fidelity_update,4);

    case 'seperate SMS'

        fidelity_update = image.*Data.sens_map;
        fidelity_update = fft2(fidelity_update,para.Recon.kSpace_size(1),para.Recon.kSpace_size(2));
        
        fidelity_update = sum(fidelity_update.*Data.SMS,5);
        fidelity_update = fidelity_update.*Data.mask;
        
        fidelity_update = Data.kSpace - fidelity_update;
        % norm for line search
        fidelity_norm = sum(abs(fidelity_update(:)).^2/prod(para.Recon.kSpace_size)/64);
        fidelity_norm = sqrt(fidelity_norm);
        
        fidelity_update = sum(fidelity_update.*conj(Data.SMS),7);

        if isfield(Data,'filter')
            fidelity_update = fidelity_update.*Data.filter;
        end
        fidelity_update = ifft2(fidelity_update);
        fidelity_update = sum(fidelity_update.*Data.sens_map_conj,4);
        
    case 'seperate SMS less memory'
        fidelity_update_all = zeros(size(image),class(image));
        fidelity_norm = 0;

        for i=1:para.Recon.no_comp
            fidelity_update = bsxfun(@times,image,Data.sens_map(:,:,:,i,:,:));

            fidelity_update = fft2(fidelity_update);
            
            fidelity_update = sum(fidelity_update.*Data.SMS,5);
            fidelity_update = bsxfun(@times,fidelity_update,Data.mask);
            fidelity_update = Data.kSpace(:,:,:,i,:,:,:) - fidelity_update;
            
            fidelity_norm = fidelity_norm + sum(abs(fidelity_update(:)).^2/prod(para.Recon.kSpace_size)/64);
        
            fidelity_update = sum(fidelity_update.*conj(Data.SMS),7);
            
            if isfield(Data,'filter')
                fidelity_update = bsxfun(@times,fidelity_update,Data.filter);
            end
        
            fidelity_update = ifft2(fidelity_update);

            fidelity_update_all = fidelity_update_all + bsxfun(@times,fidelity_update,Data.sens_map_conj(:,:,:,i,:,:));
        end
        fidelity_update = fidelity_update_all;
        fidelity_norm = sqrt(fidelity_norm);
        return
        
    case 'NUFFT less memory'
        fidelity_update_all = zeros(size(image),class(image));
        fidelity_norm = 0;
        for i=1:para.Recon.no_comp
            fidelity_update = bsxfun(@times,image,Data.sens_map(:,:,:,i,:));
            fidelity_update = NUFFT.NUFFT_new(fidelity_update,Data.N);
            if para.Recon.nSMS ~= 1
                fidelity_update = bsxfun(@times,fidelity_update,Data.phase_mod);
                fidelity_update = sum(fidelity_update,5);
                fidelity_update = Data.kSpace(:,:,:,i) - fidelity_update;
                fidelity_norm = fidelity_norm + sum(abs(fidelity_update(:)).^2)/288/288/64;
                fidelity_update = bsxfun(@times,fidelity_update,Data.phase_mod_conj);
            else
                fidelity_update = Data.kSpace(:,:,:,i) - fidelity_update;
                fidelity_norm = fidelity_norm + sum(abs(fidelity_update(:)).^2)/288/288/64;
            end
            fidelity_update = NUFFT.NUFFT_adj_new(fidelity_update,Data.N);
            fidelity_update_all = fidelity_update_all + bsxfun(@times,fidelity_update,Data.sens_map_conj(:,:,:,i,:));
        end
        fidelity_norm = sqrt(fidelity_norm);
        fidelity_update = fidelity_update_all;
        return

    case 'NUFFT'
        fidelity_update = bsxfun(@times,image,Data.sens_map);
        fidelity_update = NUFFT.NUFFT_new(fidelity_update,Data.N);
        if para.Recon.nSMS ~= 1
            fidelity_update = bsxfun(@times,fidelity_update,Data.phase_mod);
            fidelity_update = sum(fidelity_update,5);
            fidelity_update = Data.kSpace - fidelity_update;
            fidelity_norm = sum(abs(fidelity_update(:)).^2)/288/288/64;
            fidelity_update = bsxfun(@times,fidelity_update,Data.phase_mod_conj);
        else
            fidelity_update = Data.kSpace - fidelity_update;%.*Data.kwic;
            fidelity_norm = sum(abs(fidelity_update(:)).^2)/288/288/64;
        end
        fidelity_norm = sqrt(fidelity_norm);
        fidelity_update = NUFFT.NUFFT_adj_new(fidelity_update,Data.N);
        fidelity_update = bsxfun(@times,fidelity_update,Data.sens_map_conj);
        fidelity_update = sum(fidelity_update,4);

end
