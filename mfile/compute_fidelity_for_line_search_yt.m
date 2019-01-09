function fidelity_norm = compute_fidelity_for_line_search_yt(image,Data,para)

switch para.Recon.type
    case 'seperate SMS less memory'
        fidelity_norm = 0;

        for i=1:para.Recon.no_comp
            fidelity_update = bsxfun(@times,image,Data.sens_map(:,:,:,i,:,:));

            fidelity_update = fft2(fidelity_update);
            
            fidelity_update = sum(fidelity_update.*Data.SMS,5);
            
            fidelity_update = bsxfun(@times,fidelity_update,Data.mask);
            fidelity_update = Data.kSpace(:,:,:,i,:,:,:) - fidelity_update;

            fidelity_norm = fidelity_norm + sum(abs(fidelity_update(:)/para.Recon.kSpace_size(1)/8).^2);

            clear fidelity_update_temp
        end
        fidelity_norm = sqrt(fidelity_norm);
    case 'seperate SMS'

        fidelity_update = image.*Data.sens_map;

        fidelity_update = fft2(fidelity_update);
            
        fidelity_update = sum(fidelity_update.*Data.SMS,5);
        
        fidelity_update = fidelity_update.*Data.mask;
        fidelity_update = Data.kSpace - fidelity_update;

        fidelity_norm = sum(abs(fidelity_update(:)/para.Recon.kSpace_size(1)/8).^2);
        fidelity_norm = sqrt(fidelity_norm);

    case {'2D','2D less memory'}
        
        fidelity_norm = bsxfun(@times,image,Data.sens_map);
        fidelity_norm(para.Recon.image_size(1)+1:para.Recon.kSpace_size(1),para.Recon.image_size(2)+1:para.Recon.kSpace_size(2),:,:,:,:,:,:,:) = 0;
        fidelity_norm = circshift(fidelity_norm,[(para.Recon.kSpace_size(1)-para.Recon.image_size(1))/2,(para.Recon.kSpace_size(2)-para.Recon.image_size(2))/2]);
        fidelity_norm = fft2(fidelity_norm);
        fidelity_norm = bsxfun(@times,fidelity_norm,Data.mask);
        fidelity_norm = (Data.kSpace - fidelity_norm)/para.Recon.kSpace_size(1)/8;
        
    case 'NUFFT'
        fidelity_update = bsxfun(@times,image,Data.sens_map);
        fidelity_update = NUFFT.NUFFT_new(fidelity_update,Data.N);
        if para.Recon.nSMS ~= 1
            fidelity_update = bsxfun(@times,fidelity_update,Data.phase_mod);
            fidelity_update = sum(fidelity_update,5);
            fidelity_update = Data.kSpace - fidelity_update;
            fidelity_norm = sum(abs(fidelity_update(:)).^2)/288/288/64;
        else
            fidelity_update = Data.kSpace - fidelity_update;
            fidelity_norm = sum(abs(fidelity_update(:)).^2)/288/288/64;
        end
        fidelity_norm = sqrt(fidelity_norm);
        
    case 'NUFFT less memory'
        fidelity_norm = 0;
        for i=1:para.Recon.no_comp
            fidelity_update = bsxfun(@times,image,Data.sens_map(:,:,:,i,:));
            fidelity_update = NUFFT.NUFFT_new(fidelity_update,Data.N);
            if para.Recon.nSMS ~= 1
                fidelity_update = bsxfun(@times,fidelity_update,Data.phase_mod);
                fidelity_update = sum(fidelity_update,5);
                fidelity_update = Data.kSpace(:,:,:,i) - fidelity_update;
                fidelity_norm = fidelity_norm + sum(abs(fidelity_update(:)).^2)/288/288/64;
            else
                fidelity_update = Data.kSpace(:,:,:,i) - fidelity_update;
                fidelity_norm = fidelity_norm + sum(abs(fidelity_update(:)).^2)/288/288/64;
            end
        end
        fidelity_norm = sqrt(fidelity_norm);
        return


        
end