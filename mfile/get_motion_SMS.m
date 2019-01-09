function Motion = get_motion_SMS(Image,step,motion_noi)

[sx,sy,nof,nSMS] = size(squeeze(Image));

for i=1:nSMS
    M = get_motion(Image(:,:,:,i),step,motion_noi);
    Motion.idx_b(:,:,:,:,i) = M.idx_b + sx*sy*nof*(i-1);
    Motion.idx_f(:,:,:,:,i) = M.idx_f + sx*sy*nof*(i-1);
    
    Motion.xb(:,:,:,:,i) = M.xb;
    Motion.xf(:,:,:,:,i) = M.xf;
    Motion.yb(:,:,:,:,i) = M.yb;
    Motion.yf(:,:,:,:,i) = M.yf;
end

