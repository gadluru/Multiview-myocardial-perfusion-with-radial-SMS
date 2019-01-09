function theta_all = get_angle_mod(para)

if length(para.angle_mod) > 1
    theta_all = para.angle_mod;
    return
end

sy = para.Recon.sy;
sz = para.Recon.sz;
data_angle  = para.dataAngle/180*pi;
golden_angle = ((sqrt(5)-1)/2)*pi;
theta_all = zeros(1,sy*sz);

switch para.angle_mod
    
    case 2
    beta = pi/24;
    for i=1:sy*sz
        if mod(i,3) == 1
            theta_all(i) = (i-1)/3*golden_angle;
        else
            theta_all(i) = theta_all(i-1) + beta;
        end
    end
    theta_all = mod(theta_all,data_angle);
    
    case 3
    theta_all = mod((0:sy*sz-1)*golden_angle/3,data_angle);
    
    case 4
    beta = pi/36;
    for i=1:sy*sz
        if mod(i,3) == 2
            theta_all(i) = (i-2)/3*golden_angle;
        elseif mod(i,3) == 1
            theta_all(i) = (i-1)/3*golden_angle - beta;
        else
            theta_all(i) = (i-3)/3*golden_angle + beta;
        end
    end
    theta_all = mod(theta_all,data_angle);
    
    case 5
    theta_all = mod((0:sy-1)*golden_angle,data_angle);
    theta_all = repmat(theta_all,[1 sz]);
    
    case 6
    beta = pi/55;
    for i=1:sy*sz
        if mod(i,3) == 2
            theta_all(i) = (i-2)/3*golden_angle;
        elseif mod(i,3) == 1
            theta_all(i) = (i-1)/3*golden_angle - beta;
        else
            theta_all(i) = (i-3)/3*golden_angle + beta;
        end
    end
    theta_all = mod(theta_all,data_angle);
    
    case 7
    beta = pi/24;
    for i=1:sy*sz
        if mod(i,3) == 2
            theta_all(i) = (i-2)/3*golden_angle;
        elseif mod(i,3) == 1
            theta_all(i) = (i-1)/3*golden_angle - beta;
        else
            theta_all(i) = (i-3)/3*golden_angle + beta;
        end
    end
    theta_all = mod(theta_all,data_angle);
    
    case 8 % not golden angle
    dtheta0 = data_angle/sy/4;
    dtheta = data_angle/sy;
    for time = 1:sz
        th0 = mod(time-1,4) * dtheta0;
        %th0 = mod(golden_angle*(time-1),dtheta);
        thetas = th0:dtheta:data_angle;
        for ray = 1:sy
            theta_all((time-1)*sy+ray) = thetas(ray);
        end
    end
    
    case 1
    theta_all = mod((0:sy*sz-1)*golden_angle,data_angle);
    
    case 9 %tiny golden angle
        tao = (1+sqrt(5))/2;
        N = 4;
        tiny_ga = pi/(tao+N-1);
        theta_all = mod((0:sy*sz-1)*tiny_ga,data_angle);
   
    case 10
        theta_all = mod((0:sy*sz-1)*golden_angle,data_angle);
        theta_all = reshape(theta_all,[sy,sz]);
        theta_all = sort(theta_all);
        theta_all = reshape(theta_all,[1,sy*sz]);
        
    

end


end