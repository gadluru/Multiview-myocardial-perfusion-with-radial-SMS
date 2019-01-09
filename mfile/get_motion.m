function Motion = get_motion(Image,step,noi)

Image = imfilter(Image,fspecial('gaussian',size(Image,1),2));

sx_large = size(Image,1);
sy_large = size(Image,2);
Image = abs(crop_half_FOV(Image));

[sx,sy,nof] = size(Image);
[X,Y] = meshgrid(1:sx_large,1:sy_large);
X = single(X);
Y = single(Y);

dx = (sx_large-sx)/2;
dy = (sy_large-sy)/2;

center_x = dx+1:dx+sx;
center_y = dy+1:dy+sy;

for i=1:nof-1
    Ib = Image(:,:,i);
    If = Image(:,:,i+1);
    [~,yb(:,:,i),xb(:,:,i),~] = Reg_GS_tv(If,Ib,step,noi);
    [~,yf(:,:,i),xf(:,:,i),~] = Reg_GS_tv(Ib,If,step,noi);
end

Motion.yb = yb;
Motion.xb = xb;
Motion.yf = yf;
Motion.xf = xf;

yb = yb+dy;
yf = yf+dy;
xb = xb+dx;
xf = xf+dx;

Xb = repmat(X,[1,1,nof-1]);
Yb = repmat(Y,[1,1,nof-1]);
Xf = Xb;
Yf = Yb;

Xb(center_x,center_y,:) = gather(xb);
Yb(center_x,center_y,:) = gather(yb);
Xf(center_x,center_y,:) = gather(xf);
Yf(center_x,center_y,:) = gather(yf);
%{
if GDFlag
    signal_sum = squeeze(sum(sum(Image,1),2));
    signal_diff = diff(signal_sum);
    [GD_frames(:,1),GD_frames(:,2)] = findpeaks(signal_diff);
    GD_frames = sortrows(GD_frames,'descend');
    GD_frames = GD_frames(1:2,2);
    %GD_frame = Image(:,:,GD_frames(1)+1);
    %heart_mask = bwdist(GD_frame==max(max(GD_frame)))<25 & GD_frame>max(max(GD_frame))*0.6;
    %heart_mask = refine_mask_edge(heart_mask);
    GD_frames = [GD_frames(1)-2:GD_frames(1)+2,GD_frames(2)-2:GD_frames(2)+2];
    NGD = length(GD_frames);
    Xb(:,:,GD_frames) = repmat(X,[1,1,NGD]);
    Xf(:,:,GD_frames) = repmat(X,[1,1,NGD]);
    Yb(:,:,GD_frames) = repmat(Y,[1,1,NGD]);
    Yf(:,:,GD_frames) = repmat(Y,[1,1,NGD]);
end
%}
Motion.yb = Yb;
Motion.xb = Xb;
Motion.yf = Yf;
Motion.xf = Xf;

Motion.idx_b = sub2ind([sx_large,sy_large],round(Yb),round(Xb));
Motion.idx_f = sub2ind([sx_large,sy_large],round(Yf),round(Xf));

add_idx = (0:nof-2)*sx_large*sy_large;
%add_idx(end) = add_idx(end-1);
add_idx = permute(add_idx,[3,1,2]);

Motion.idx_b = Motion.idx_b + add_idx + sx_large*sy_large;
Motion.idx_f = Motion.idx_f + add_idx;

%{
Xb = repmat(X,[1,1,nof]);
Yb = repmat(Y,[1,1,nof]);
Xf = Xb;
Yf = Yb;

Xb(center_x,center_y,1:end-1) = xb;
Yb(center_x,center_y,1:end-1) = yb;
Xf(center_x,center_y,2:end) = xf;
Yf(center_x,center_y,2:end) = yf;

add_idx = (0:nof)*sx_large*sy_large;
add_idx(end) = add_idx(end-1);
add_idx = permute(add_idx,[3,1,2]);



Motion.idx_b = sub2ind([sx_large,sy_large],round(Yb),round(Xb))+add_idx(:,:,2:end);
Motion.idx_f = sub2ind([sx_large,sy_large],round(Yf),round(Xf))+add_idx(:,:,1:end-1);
%}

