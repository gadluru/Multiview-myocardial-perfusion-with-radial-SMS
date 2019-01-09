function showImage(Image,Cost)
Image = crop_half_FOV(abs(Image));
[nx,ny,nof,~,nSMS,ns] = size(Image);
frame_num = floor(nof/4);
im = reshape(Image,[nx,ny,nof,nSMS*ns]);
if frame_num ~= 0
    im = im(:,:,[frame_num frame_num*2 frame_num*3],:);
    im = permute(im,[1 3 2 4]);
    im = reshape(im,[nx*3 ny*nSMS*ns]);
else
    im = im(:,:,1,:);
    im = squeeze(im);
    im = reshape(im,[nx ny*nSMS*ns]);
end

figure(100);
clf;
subplot(2,2,[1 3])
imagesc(im)
colormap gray
brighten(0.3)
axis image
axis off
        
subplot(2,2,[2 4])

hold on;
plot(Cost.totalCost);
plot(Cost.fidelityNorm,'c*-')
plot(Cost.temporalNorm,'kx-');
plot(Cost.spatialNorm,'k.-');
legend('Total Cost','Fidelity Norm','Temporal Norm','Spatial Norm')
%axis([1 lengthCost.totalCost)+1 0 max(Cost.totalCost)])
drawnow
end
