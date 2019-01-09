function show_yt_pause(img)
img = squeeze(img);
figure
for i=1:size(img,3)
    imagesc(abs(img(:,:,i)))
    colormap gray
    axis image
    brighten(0.4)
    title(i)
    %drawnow
    pause
end
end