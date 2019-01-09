function G = init(kSpace_radial,kx,ky,over_sampling,core_size)
% G = init(kSpace_radial,kx,ky,over_sampling,core_size)

sx  = size(kSpace_radial,1);

nx = round(sx/8);
[Gx,Gy] = GROG.get_Gx_Gy(kSpace_radial(nx*3+1:5*nx,:,:,:), kx(nx*3+1:5*nx,:,:,:), ky(nx*3+1:5*nx,:,:,:));

skx = size(kx,1);
nor = size(kx,2);
nof = size(kx,3);
nc  = size(Gx,1);

sx_over = round(sx*over_sampling);

kx = kx*over_sampling;
ky = ky*over_sampling;
Gx = Gx^(1/over_sampling);
Gy = Gy^(1/over_sampling);

max_shift_x = core_size(1)/2;
max_shift_y = core_size(2)/2;

GxDict_size = max_shift_x*200 + 1;
GyDict_size = max_shift_y*200 + 1;
GxDict = single(zeros([1 nc nc GxDict_size]));
GyDict = single(zeros([nc nc 1 1 GyDict_size]));

dx = -max_shift_x:0.01:max_shift_x;
dy = -max_shift_y:0.01:max_shift_y;

for di=1:GxDict_size
    GxDict(1,:,:,di) = Gx^dx(di);
end
for di=1:GyDict_size
    GyDict(:,:,1,1,di) = Gy^dy(di);
end

G_r2c_Dict = bsxfun(@times,GxDict,GyDict);
G_r2c_Dict = squeeze(sum(G_r2c_Dict,2));
G_r2c_Dict = reshape(G_r2c_Dict,[nc nc GxDict_size*GyDict_size]);
G_r2c_Dict = permute(G_r2c_Dict,[3 1 2]);

GxDict = permute(GxDict,[2 3 1 4]);
GyDict = permute(GyDict,[3 1 2 4 5]);
G_c2r_Dict = bsxfun(@times,GxDict,GyDict);
G_c2r_Dict = squeeze(sum(G_c2r_Dict,2));
G_c2r_Dict = reshape(G_c2r_Dict,[nc nc GxDict_size*GyDict_size]);
G_c2r_Dict = permute(G_c2r_Dict,[3 1 2]);

x_cart = zeros(skx,nor,nof,core_size(1));
y_cart = zeros(skx,nor,nof,core_size(2));

if mod(core_size(1),2) == 0
    x_cart(:,:,:,1) = floor(kx) - core_size(1)/2 + 1;
else
    x_cart(:,:,:,1) = round(kx) - (core_size(1)-1)/2;
end
add_x(1,1,1,:) = 1:core_size(1)-1;
x_cart(:,:,:,2:end) = bsxfun(@plus,x_cart(:,:,:,1),add_x);

if mod(core_size(2),2) == 0
    y_cart(:,:,:,1) = floor(ky) - core_size(2)/2 + 1;
else
    y_cart(:,:,:,1) = round(ky) - (core_size(2)-1)/2;
end
add_y(1,1,1,:) = 1:core_size(2)-1;
y_cart(:,:,:,2:end) = bsxfun(@plus,y_cart(:,:,:,1),add_y);

x_cart = repmat(x_cart,1,1,1,core_size(2));
x_cart = reshape(x_cart,skx,nor,nof,prod(core_size));
y_cart = permute(y_cart,[1,2,3,5,4]);
y_cart = repmat(y_cart,1,1,1,core_size(1));
y_cart = reshape(y_cart,skx,nor,nof,prod(core_size));

dx = round((x_cart - kx)*100)/100;
dy = round((y_cart - ky)*100)/100;

xDict_r2c = round((dx+max_shift_x)*100)+1;
yDict_r2c = round((dy+max_shift_y)*100)+1;

index_Dict_r2c = sub2ind([GxDict_size GyDict_size],xDict_r2c,yDict_r2c);

xDict_c2r = round(-dx+max_shift_x)*100+1;
yDict_c2r = round(-dy+max_shift_y)*100+1;

index_Dict_c2r = sub2ind([GxDict_size GyDict_size],xDict_c2r,yDict_c2r);

if mod(core_size(1),2) == 0
    x_cart = x_cart+sx_over/2+core_size(1)/2;
else
    x_cart = x_cart+sx_over/2+core_size(1)/2+0.5;
end
if mod(core_size(2),2) == 0
    y_cart = y_cart+sx_over/2+core_size(1)/2;
else
    y_cart = y_cart+sx_over/2+core_size(1)/2+0.5;
end
%%%
x_cart = mod(x_cart,sx_over); x_cart(x_cart==0) = sx_over;
y_cart = mod(y_cart,sx_over); y_cart(y_cart==0) = sx_over;
%%%
if mod(core_size(1),2) == 0
    indx = sub2ind([(sx_over+core_size(1)),(sx_over+core_size(2)),nof],x_cart,y_cart);
else
    indx = sub2ind([(sx_over+core_size(1)+1),(sx_over+core_size(2)+1),nof],x_cart,y_cart);
end

rad2cart = cell(1,nof);
rad_num = (1:skx*nor*prod(core_size)).';

weight_scale = sqrt(max_shift_x^2+max_shift_y^2)*1.00001;
weight = weight_scale - sqrt(dx.^2 + dy.^2);

weight = permute(weight,[1,2,4,3]);
weight = reshape(weight,sx*nor*prod(core_size),nof);

if mod(core_size(1),2) == 0
    weight_cart = single(zeros((sx_over+core_size(1))*(sx_over+core_size(2)),nof));
else
    weight_cart = single(zeros((sx_over+core_size(1)+1)*(sx_over+core_size(2)+1),nof));
end

for i=1:nof
    indx_temp = indx(:,:,i,:);
    indx_temp = indx_temp(:);
    if mod(core_size(1),2) == 0
        rad2cart{i} = sparse(indx_temp,rad_num,1,(sx_over+core_size(1))*(sx_over+core_size(2)),skx*nor*prod(core_size));
    else
        rad2cart{i} = sparse(indx_temp,rad_num,1,(sx_over+core_size(1)+1)*(sx_over+core_size(2)+1),skx*nor*prod(core_size));
    end
    weight_cart(:,i) = rad2cart{i} * double(weight(:,i));
    
    weight_temp = full(sum(rad2cart{i},2));
    weight_temp(weight_temp ==0) = 1;
    weight3(:,i) = 1./weight_temp;
    clear indx_temp
end

G.rad2cart = rad2cart;
G.Dict_r2c = G_r2c_Dict;
G.Dict_c2r = G_c2r_Dict;
G.indx_r2c = index_Dict_r2c;
G.indx_c2r = index_Dict_c2r;
G.core_size = core_size;
G.sx_over = sx_over;
G.siz = [skx nor nof nc];
G.weight = weight3;
G.weight1 = reshape(weight,[skx*nor,prod(core_size),nof]);
G.weight1 = permute(G.weight1,[1 3 2]);
weight_cart(weight_cart==0)=1;
G.weight2 = 1./weight_cart;

G.W = density_comp_area(kx,ky,mod(atan(ky(1,:,:)./kx(1,:,:)),pi));
G.Gx = Gx;
G.Gy = Gy;
