function [I2,x,y,E] = Reg_GS_tv(I1,I0,step,noi)
[sx,sy] = size(I0);
Maxiter = noi;
a = 30;
s = 3;

kx = cos(2*pi*(0:sx/(sx-1):sx));
ky = cos(2*pi*(0:sy/(sy-1):sy));
W = 2*(kx+ky.'-2);
W = single((1-a*W).^-s);
 
[y,x] = meshgrid(1:sx,1:sy);

for iter=1:Maxiter
    
    I2 = interp2(I1,y,x);
    
    ddx = 0.5*(I2(3:end,:) - I2(1:end-2,:));
    ddy = 0.5*(I2(:,3:end) - I2(:,1:end-2));
    ddx = cat(1,I2(2,:) - I2(1,:),ddx,I2(end,:) - I2(end-1,:));
    ddy = cat(2,I2(:,2) - I2(:,1),ddy,I2(:,end) - I2(:,end-1));
    
    dI = I2-I0;
    dI = dI - medfilt2(dI,[3,3]);

    J = jacobian(y,x);

    dx = W.*fft2(ddx.*dI.*J);
    dx = -real(ifft2(dx));

    dy = W.*fft2(ddy.*dI.*J);
    dy = -real(ifft2(dy));

    d = sqrt(dx.^2+dy.^2);

    md = max(d(:));
 
    dx = dx/md;
    dy = dy/md;
    
    E(iter) = sum(d(:));
    if iter>1 && E(iter)>E(iter-1)
        step = step/2;
        %return
    end
    %{
    figure(2)
    plot(E)
    drawnow
    %}
    x = x + step*dx;
    y = y + step*dy;
    
    x(x<1) = 1;
    x(x>sx) = sx;
    y(y<1) = 1;
    y(y>sy) = sy;

%     figure(1)
%     subplot(2,1,1)
%     imagesc([I2,I0,dI])
%     colormap gray
%     axis image
%     brighten(0.4)
%     drawnow

end

function J = jacobian(x,y)
    %[Y,X] = meshgrid(1:sx,1:sy);
    [fxx,fxy] = gradient(x);
    [fyx,fyy] = gradient(y);
    J = fxx .* fyy - fxy .* fyx;
return