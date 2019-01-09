function W = density_comp_area(kx,ky,theta_in)

[sx,nor,nof] = size(kx);

if nor==1
    W = 1:-2/sx:0;
    W = [W(1:sx/2),fliplr(W(1:sx/2))].';
    return
end

d = sqrt(kx.^2 + ky.^2);
ring_area = pi * ((d + 0.5).^2 - (d - 0.5).^2);
ring_area(d==0) = pi*0.5^2/nor;
theta = permute(theta_in,[2,1,3]);
theta(:,2,:) = repmat((1:nor)',1,1,nof);
for i=1:nof
    theta(:,:,i) = sortrows(theta(:,:,i),1);
end

for i=1:nor
    if i==1
        theta(i,3,:) = pi + theta(i+1,1,:) - theta(nor,1,:);
    elseif i==nor
        theta(i,3,:) = pi - (theta(i-1,1,:) - theta(1,1,:));
    else
        theta(i,3,:) = theta(i+1,1,:) - theta(i-1,1,:);
    end
end


theta(:,1,:) = [];
for i=1:nof
    theta(:,:,i) = sortrows(theta(:,:,i),1);
end
theta(:,1,:) = [];
theta = permute(theta,[2,1,3]);
W = bsxfun(@times,ring_area,theta/2/pi);
W = W / max(W(:));