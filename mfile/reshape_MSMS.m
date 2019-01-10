function Image = reshape_MSMS(Image)
% reshape multi-slice SMS image

[sx,sy,nof,ns,nsms] = size(Image);
Image = permute(Image,[1 4 2 5 3]);
Image = reshape(Image,[sx*ns,sy*nsms,nof]);