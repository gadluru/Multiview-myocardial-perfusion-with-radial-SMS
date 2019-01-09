function kSpace = coil_unstreaking_and_compression_SMS(kSpace,kSpace_info)
%in kSpace : [sx,nor,coil,nof,set]
idx = squeeze(sum(sum(sum(sum(kSpace)))) == 0);
kSpace(:,:,:,:,idx) = [];
kSpace = kSpace * 10^8;

kSpace = permute(kSpace,[1,2,4,3,5]);
%kSpace = phase_correction_022718(kSpace);
%kSpace = phase_correction_031218(kSpace,kSpace_info.phase_mod);
kSpace = permute(kSpace,[1,2,4,3,5]);
kSpace_raw = kSpace;
kSpace = permute(kSpace,[1,2,4,5,3]);
[sx,nor,nof,ns,nc] = size(kSpace);
kSpace = reshape(kSpace,[sx,nor*nof,ns,nc]);

if isfield(kSpace_info,'ray_mod')
    notwanted = kSpace_info.ray_mod(:,1)==0;
else
    PD = kSpace_info.ProtonDensityScans;
    notwanted = false(1,nor*nof);
    notwanted(1:PD*nor) = true;
end

kSpace(:,notwanted,:,:) = [];
phase_mod = kSpace_info.phase_mod(~notwanted,:);
angle_mod = kSpace_info.angle_mod(~notwanted,:);

%nor_ref = length(angle_mod);

nor_one = nor*3;
nor_ref = nor_one*10;
nor_ref = min(nor_ref,size(kSpace,2));
%nor_ref = 900;
threshold = 20;
%ratio = nor_ref/nor_one;

kSpace = kSpace(:,1:nor_ref,:,:);
kSpace = reshape(kSpace,[sx,nor_ref,1,ns,nc]);
phase_mod = phase_mod(1:nor_ref,:);
angle_mod = angle_mod(1:nor_ref,:);

phase(1,:,1,:,1,1) = exp(-1i*phase_mod*0);
phase(1,:,1,:,1,2) = exp(-1i*phase_mod*2*pi/3);
phase(1,:,1,:,1,3) = exp(-1i*phase_mod*4*pi/3);

angle_mod = permute(angle_mod,[3,1,2]);
[kx,ky] = get_k_coor(sx,angle_mod,0,round((sx+1)/2));

for i=1:3
    if size(angle_mod,3)>1
        for j=1:ns
            N{j} = NUFFT.init_new(kx(:,:,j),ky(:,:,j),1,[4,4]);
            N_one{j} = NUFFT.init_new(kx(:,1:nor_one,j),ky(:,1:nor_one,j),1,[4,4]);
            im_ref(:,:,j,:,i) = squeeze(NUFFT.NUFFT_adj_new(kSpace(:,:,:,j,:).*phase(:,:,:,j,:,i),N{j}));
            im_one(:,:,j,:,i) = squeeze(NUFFT.NUFFT_adj_new(kSpace(:,1:nor_one,:,j,:).*phase(:,1:nor_one,:,j,:,i),N_one{j}));
        end
    else
        N = NUFFT.init_new(kx,ky,1,[4,4]);
        N_one = NUFFT.init_new(kx(:,1:nor_one),ky(:,1:nor_one),1,[4,4]);
        for j=1:ns
            im_ref(:,:,j,:,i) = squeeze(NUFFT.NUFFT_adj_new(kSpace(:,:,:,j,:).*phase(:,:,:,1,:,i),N));
            im_one(:,:,j,:,i) = squeeze(NUFFT.NUFFT_adj_new(kSpace(:,1:nor_one,:,j,:).*phase(:,1:nor_one,:,1,:,i),N_one));
        end
    end
end
%im_sos = sos(im_ref,4);
%im_sos = permute(im_sos,[1,2,4,3,5]);

for i=1:ns
    for j=1:3
        im_one_temp = squeeze(im_one(:,:,i,:,j));
        im_ref_temp = squeeze(im_ref(:,:,i,:,j));
        %scale_temp = max(max(abs(im_ref_temp)))./max(max(abs(im_one_temp)));
        scale_temp = nor_ref/nor_one;
        n_temp = abs(im_one_temp.*scale_temp - im_ref_temp).^2;
        d_temp = abs(crop_half_FOV(im_ref_temp)).^2;
        score_temp(i,j,:) = sum(sum(n_temp))./sum(sum(d_temp));
    end
end

score = squeeze(sum(score_temp,2));
score = score/min(score(:));
score(score<threshold) = 1;
score(score>threshold) = 1./score(score>threshold)*threshold;
score = permute(score,[1,3,2]);
%{
figure,hold on
for i=1:size(score,1)
    plot(squeeze(score(i,:,:)))
end
title 'Coil Weight'
xlabel 'Coil Count'
ylabel 'Weight'
for i=1:size(score,1)
    Legend(i,:) = ['SMS set ', num2str(i)];
end
legend(Legend)
drawnow
%}
show_3D(im_one,score)
%N = round(nc/2);
%{
for i=1:3
    score_temp = score(i,:);
    score_temp = score_temp/min(score_temp);
    [~,order] = sort(score_temp,'descend');
    threshold = score_temp(order(9));
    %if threshold<1.3
        %threshold = 1.3;
    %end
    idx = score_temp>threshold;
    score(i,idx) = 1./score_temp(idx);
    score(i,idx) = score(i,idx)./max(score(i,idx))*0.7;
    score(i,~idx) = 1;
end
score = permute(score,[1,3,2]);
%}
clear kSpace
for i=1:ns
    kSpace_temp = kSpace_raw(:,:,:,:,i).*score(i,1,:);
    kSpace_temp = permute(kSpace_temp,[1,2,4,3]);
    kSpace_temp = reshape(kSpace_temp,[sx*nor*nof,nc]);
    coeff = pca(kSpace_temp);
    kSpace_temp = kSpace_temp*coeff(:,1:8);
    kSpace_temp = reshape(kSpace_temp,[sx,nor,nof,8]);
    kSpace(:,:,:,:,i) = kSpace_temp;
end
end

function show_3D(img,score)

    img = sum(abs(img),5);
    img = permute(img,[1,2,4,3]);
    [sx,~,ns,nf] = size(img);
    ns_sqrt = ceil(sqrt(ns));
    if ns_sqrt^2~=ns
        img(:,:,ns+1:ns_sqrt^2,:) = zeros(sx,sx,ns_sqrt^2-ns,nf);
        ns = ns_sqrt^2;
    end
    img = reshape(img,sx,sx*ns_sqrt,ns_sqrt,nf);
    img = permute(img,[1 3 2 4]);
    img = reshape(img,sx*ns_sqrt,sx*ns_sqrt,nf);

    for i=1:size(img,3)
        figure
        imagesc(abs(img(:,:,i)))
        colormap gray
        axis image
        axis off
        brighten(0.4)
        title(['SMS set ',num2str(i)])
        score_temp = squeeze(score(i,:,:));
        [x,y] = ind2sub([ns_sqrt,ns_sqrt],1:length(score_temp));
        for j=1:length(score_temp)
            text((x(j)-1)*sx+sx/2,(y(j)-1)*sx+sx/10,num2str(score_temp(j),'%.2f'),'Color',[0.85,0.325,0.098]);
        end
        drawnow
    end

end


