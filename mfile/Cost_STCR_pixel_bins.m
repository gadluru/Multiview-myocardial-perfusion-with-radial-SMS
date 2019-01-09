function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR_pixel_bins(fUpdate, Image, Data, para, Cost_old)

%N = numel(Image);
sWeight = para.Recon.weight_sTV;
tWeight = para.Recon.weight_tTV;
bins = para.Recon.bins;

fNorm = sum(abs(fUpdate(:)).^2);

if tWeight ~= 0
    tNorm = abs(crop_half_FOV(Image(Data.Motion.idx_b) - Image(:,:,1:end-1,:,:)));
    tNorm = mean(tWeight) * sum(tNorm(:))*0.75;
    
    for i=1:size(bins,1)
        bin_temp = bins(i,:);
        Image_temp = Image(:,:,bin_temp,:,:);
        Motion_temp = Data.Motion_bins{i};
        tNorm_temp = abs(crop_half_FOV(Image_temp(Motion_temp.idx_b) - Image_temp(:,:,1:end-1,:,:)));
        %tNorm_temp = abs(crop_half_FOV(diff(Image_temp,1,3)));
        tNorm = tNorm + mean(tWeight) * sum(tNorm_temp(:))*0.75;
    end
    
else
    tNorm = 0;
end
Image = crop_half_FOV(Image);

if sWeight ~= 0
    sx_norm = abs(diff(Image,1,2));
    sx_norm(:,end+1,:,:,:)=0;
    sy_norm = abs(diff(Image,1,1));
    sy_norm(end+1,:,:,:,:)=0;
    sNorm = sqrt(abs(sx_norm).^2+abs(sy_norm).^2);
    sNorm = mean(sWeight) * sum(sNorm(:));
else
    sNorm = 0;
end

Cost = sNorm + tNorm + fNorm;

if nargin == 4
    Cost_new = Cost;
    return
end

Cost_new = Cost_old;

if isempty(Cost_old.fidelityNorm)==1
    Cost_new.fidelityNorm = gather(fNorm);
    Cost_new.temporalNorm = gather(tNorm);
    Cost_new.spatialNorm = gather(sNorm);
    Cost_new.totalCost = gather(Cost);
else    
    Cost_new.fidelityNorm(end+1) = gather(fNorm);
    Cost_new.temporalNorm(end+1) = gather(tNorm);
    Cost_new.spatialNorm(end+1) = gather(sNorm);
    Cost_new.totalCost(end+1) = gather(Cost);
end

end