function tTV_update = compute_tTV_pixel(Image,weight,beta_square,Motion)

temp_f = Image(Motion.idx_b) - Image(:,:,1:end-1,:,:);
temp_b = Image(:,:,2:end,:,:) - Image(Motion.idx_f);
%{
    for i=1:size(temp_f,3)
        for j=1:size(temp_f,5)
            temp_f(:,:,i,1,j) = temp_f(:,:,i,1,j) - medfilt2(real(temp_f(:,:,i,1,j)),[3,3]) - 1i*medfilt2(imag(temp_f(:,:,i,1,j)),[3,3]);
            temp_b(:,:,i,1,j) = temp_b(:,:,i,1,j) - medfilt2(real(temp_b(:,:,i,1,j)),[3,3]) - 1i*medfilt2(imag(temp_b(:,:,i,1,j)),[3,3]);
        end
    end
%}
temp_f = temp_f./sqrt(abs(temp_f).^2+beta_square);
temp_b = temp_b./sqrt(abs(temp_b).^2+beta_square);
%temp_f = real(temp_f)./sqrt(real(temp_f).^2+beta_square) + 1i*imag(temp_f)./sqrt(imag(temp_f).^2+beta_square);
%temp_b = real(temp_b)./sqrt(real(temp_b).^2+beta_square) + 1i*imag(temp_b)./sqrt(imag(temp_b).^2+beta_square);

tTV_update = weight.*cat(3,temp_f(:,:,1,:,:),temp_f(:,:,2:end,:,:)-temp_b(:,:,1:end-1,:,:),-temp_b(:,:,end,:,:));

return

nof = size(Image,3);
for i=2:nof-1
    I0 = Image(:,:,i);
    If = Image(:,:,i+1);
    Ib = Image(:,:,i-1);

    %yb_temp = round(Motion.yb(:,:,i));
    %xb_temp = round(Motion.xb(:,:,i));
    
    %yf_temp = round(Motion.yf(:,:,i-1));
    %xf_temp = round(Motion.xf(:,:,i-1));
    
    %idx_b = sub2ind([288,288],yb_temp,xb_temp);
    %idx_f = sub2ind([288,288],yf_temp,xf_temp);
    
    %temp_f = If(idx_f) - I0;
    %temp_b = I0 - Ib(idx_b);
    
    %temp_f = If(Motion.idx_f(:,:,i-1)) - I0;
    %temp_b = I0 - Ib(Motion.idx_b(:,:,i));
    
    temp_f = If(Motion.idx_b(:,:,i)) - I0;
    temp_b = I0 - Ib(Motion.idx_f(:,:,i-1));
    
    temp_f = temp_f./sqrt(abs(temp_f).^2+beta_square);
    temp_b = temp_b./sqrt(abs(temp_b).^2+beta_square);
    
    tTV_update(:,:,i) = temp_f-temp_b;
end
tTV_update(:,:,end+1) = 0;
tTV_update = weight*tTV_update;
