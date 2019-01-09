function [bins,cardiac_signal] = self_gating_image_space(Image)
Image = abs(Image);
[sx,sy,nof,set,sms] = size(Image);
%signal = reshape(Image,[sx,sy,nof,set,sms]);
signal = Image./mean(mean(mean(Image,1),2),3);
max_t = max(signal,[],3);
signal_normalized = signal./max_t;
signal_sum = sum(signal_normalized(:,:,1:10,:,:),3);
signal_sum = reshape(signal_sum,[sx*sy,1,1,set,sms]);
signal_sum = sort(signal_sum,1);
threshold_mask1 = signal_sum(round(sx*sx*0.2),:,:,:,:);
mask1 = sum(signal_normalized(:,:,1:10,:,:),3) < threshold_mask1;
max_all = max(max(max_t));
mask2 = max_t > max_all*0.40;
mask_combined = mask1 & mask2;
mask_final = false(sx,sy,1,set,sms);
for i=1:set
    for j=1:sms
        mask_temp = false(sx,sy);
        CC = bwconncomp(mask_combined(:,:,:,i,j));
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = maxk(numPixels,3);
        mask_temp(CC.PixelIdxList{idx(1)}) = true;
        for ii=2:length(idx)
            if  numPixels(idx(ii)) / numPixels(idx(1)) > 0.3
                [x_temp,y_temp] = ind2sub([sx,sy],CC.PixelIdxList{idx(ii)});
                if sqrt((mean(x_temp)-sx/2)^2 + (mean(y_temp)-sy/2)^2) < 40
                    mask_temp(CC.PixelIdxList{idx(ii)}) = true;
                end
            end
        end
        mask_temp = bwdist(mask_temp);
        mask_final(:,:,1,i,j) = mask_temp<5;
    end
end
cardiac_signal = squeeze(sum(sum(sum(mask_final.*Image)),5));
close all

for i=1:set
    signal_all(:,i) = normalize_to_0_1(normalize_to_0_1(cardiac_signal(:,i)));
end

bins = zeros(size(signal_all),'single');
sys_threshold = 0.2;
dia_threshold = 0.8;

bins(signal_all < sys_threshold) = 1;
bins(signal_all > dia_threshold) = 2;

for i=1:set
    while sum(bins(:,i)==1) < nof/2*0.8
        sys_threshold = sys_threshold + 0.05;
        bins(signal_all < sys_threshold) = 1;
    end
    while sum(bins(:,i)==2) < nof/2*0.8
        dia_threshold = dia_threshold - 0.05;
        bins(signal_all > dia_threshold) = 2;
    end
end
bins(bins==0) = 3;
