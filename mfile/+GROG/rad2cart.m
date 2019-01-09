function kSpace_cart = rad2cart(kSpace_radial,G)

sx = G.siz(1);
nor = G.siz(2);
nof = G.siz(3);
nc = G.siz(4);
sx_over = G.sx_over;
core_size = G.core_size;

kSpace_radial = reshape(kSpace_radial,[sx*nor*nof,1,1,nc]);
if mod(core_size(1),2) == 0
    kSpace_cart = zeros([(sx_over+core_size(1))*(sx_over+core_size(2)) nof nc],'like',kSpace_radial);
else
    kSpace_cart = zeros([(sx_over+core_size(1)+1)*(sx_over+core_size(2)+1) nof nc],'like',kSpace_radial);
end

G_all = G.Dict_r2c(G.indx_r2c,:,:); 
G_all = reshape(G_all,[sx*nor*nof,prod(core_size),nc,nc]);

k_target = bsxfun(@times,G_all,kSpace_radial);
k_target = sum(k_target,4);

k_target = reshape(k_target,[sx*nor,nof,prod(core_size),nc]);

k_target = bsxfun(@times,k_target,G.weight1); % change weight here

k_target = permute(k_target,[1 3 4 2]);
k_target = reshape(k_target,[sx*nor*prod(core_size),nc,nof]);

for i=1:nof
    kSpace_cart(:,i,:) = single(G.rad2cart{i}*double(k_target(:,:,i)));
end

kSpace_cart = bsxfun(@times,kSpace_cart,G.weight2); 
if mod(core_size(1),2) == 0
    kSpace_cart = reshape(kSpace_cart,[sx_over+core_size(1) sx_over+core_size(2) nof nc]);
else
    kSpace_cart = reshape(kSpace_cart,[sx_over+core_size(1)+1 sx_over+core_size(2)+1 nof nc]);
end

