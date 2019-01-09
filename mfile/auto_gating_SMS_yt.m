function para = auto_gating_SMS_yt(Image_ref,para)

if ~para{1}.gated
    cardiac_signal = self_gating_image_space(crop_half_FOV(Image_ref));
    respiration_signal = get_respiration_bins(crop_half_FOV(Image_ref));
    for i=1:size(Image_ref,4)
        para{i}.Recon.bins(1,:) = cardiac_signal(:,i)==1 & respiration_signal==1;
        para{i}.Recon.bins(2,:) = cardiac_signal(:,i)==2 & respiration_signal==1;
        para{i}.Recon.bins(3,:) = cardiac_signal(:,i)==3 & respiration_signal==1;
        
        para{i}.Recon.bins(4,:) = cardiac_signal(:,i)==1 & respiration_signal==2;
        para{i}.Recon.bins(5,:) = cardiac_signal(:,i)==2 & respiration_signal==2;
        para{i}.Recon.bins(6,:) = cardiac_signal(:,i)==3 & respiration_signal==2;
        
        para{i}.Recon.bins(7,:) = cardiac_signal(:,i)==1 & respiration_signal==3;
        para{i}.Recon.bins(8,:) = cardiac_signal(:,i)==2 & respiration_signal==3;
        para{i}.Recon.bins(9,:) = cardiac_signal(:,i)==3 & respiration_signal==3;
        
        para{i}.Recon.bins(10,:) = cardiac_signal(:,i)==1 & respiration_signal==4;
        para{i}.Recon.bins(11,:) = cardiac_signal(:,i)==2 & respiration_signal==4;
        para{i}.Recon.bins(12,:) = cardiac_signal(:,i)==3 & respiration_signal==4;

        idx = sum(para{i}.Recon.bins,2)<=2;
        para{i}.Recon.bins(idx,:) = [];
        
        para{i}.respiration_signal = respiration_signal;
        para{i}.cardiac_signal = cardiac_signal;
    end
else 
    respiration_signal = get_respiration_bins(crop_half_FOV(Image_ref));
    for i=1:size(Image_ref,4)
        para{i}.Recon.bins(1,:) = respiration_signal==1;
        para{i}.Recon.bins(2,:) = respiration_signal==2;
        para{i}.Recon.bins(3,:) = respiration_signal==3;
        para{i}.Recon.bins(4,:) = respiration_signal==4;
        
        para{i}.respiration_signal = respiration_signal;
    end
end