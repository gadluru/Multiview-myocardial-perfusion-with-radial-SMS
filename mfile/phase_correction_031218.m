function kSpace = phase_correction_031218(kSpace,phase_mod)
[nx,nor,nof,nc,nset] = size(kSpace);
center_idx = round(nx/2)+1;

coil_idx = squeeze(sum(sum(sum(abs(kSpace(center_idx,:,:,:,:)).^2,2),3),5));
[~,coil_idx] = max(coil_idx);

mod0 = phase_mod == 0;
mod0 = reshape(mod0,[1,nor,nof,1,size(mod0,2)]);
N = floor(nor/3)*3;
mod0(:,N+1:end,:,:,:) = false;

center = kSpace(center_idx,:,:,coil_idx,:).*mod0;
center_0 = center; % save for plot

phase_all = -pi:0.01:pi;
phase_all = exp(1i*phase_all).';

phase_shift = phase_all.*(center);
phase_shift_diff = angle(phase_shift);% - center_phase_frame_1;
phase_shift_diff_sos = sum(phase_shift_diff.^2,2);
[~,idx] = min(phase_shift_diff_sos);
phase_correction = phase_all(idx);
if size(phase_correction,1)~=1
    phase_correction = permute(phase_correction,[2 3 1]);
end
kSpace = phase_correction.*kSpace;
% center = kSpace(center_idx,:,:,coil_idx,:).*mod0;
% figure,plot(angle(center_0(:)))
% hold on
% plot(angle(center(:)))
% title 'Phase Correction Result'
% xlabel 'phase mod 0 measurements'
% ylabel 'phase'
% legend 'before' 'after'
% drawnow