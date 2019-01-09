function phase_mod_all = get_phase_mod(para)

sy = para.Recon.sy;
sz = para.Recon.sz;

if length(para.phase_mod) == sy*sz && sum(para.phase_mod(:))~=0
    phase_mod_all(:,1) = ones(size(para.phase_mod));
    phase_mod_all(:,2) = exp(-1i*2/3*pi*para.phase_mod);
    phase_mod_all(:,3) = conj(phase_mod_all(:,2));
    idx = para.phase_mod == 3;
    phase_mod_all(idx,:) = 0;
    phase_mod_all = single(phase_mod_all);
    phase_mod_all = reshape(phase_mod_all,[1,sy*sz,1,3]);
    return
elseif length(para.phase_mod) == sy*sz && sum(para.phase_mod(:))==0
    phase_mod_all = ones([1 sy*sz]);
    return
end

switch para.phase_mod
    case 1
        phase_mod_all(1,:,1,1) = ones([1 sy*sz]);
        phase_mod_all(1,:,1,2) = repmat(exp(-1i*2/3*pi*(0:sy-1)),[1 sz]);
        phase_mod_all(1,:,1,3) = repmat(exp(-1i*4/3*pi*(0:sy-1)),[1 sz]);
        
    case 1.5
        phase_mod_all(1,:,1,1) = ones([1 sy*sz]);
        phase_mod_all(1,:,1,2) = exp(-1i*2/3*pi*(0:sy*sz-1));
        phase_mod_all(1,:,1,3) = exp(-1i*4/3*pi*(0:sy*sz-1));
        
    case 2
    phase_mod_all(1,:,1,1) = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    phase_mod_all(1,[1 4 8 11 15 18],1,2) = 1;
    phase_mod_all(1,[1 4 8 11 15 18],1,3) = 1;

    phase_mod_all(1,[3 7 10 13 14 17],1,2) = exp(-1i*2/3*pi);
    phase_mod_all(1,[3 7 10 13 14 17],1,3) = exp(-1i*4/3*pi);

    phase_mod_all(1,[2 5 6 9 12 16],1,2) = exp(-1i*4/3*pi);
    phase_mod_all(1,[2 5 6 9 12 16],1,3) = exp(-1i*2/3*pi);

    phase_mod_all(1,19:20,:,:) = zeros(1,2,1,3);
    phase_mod_all = repmat(phase_mod_all,[1 100 1 1]);
    case 3
    phase_mod_all(1,:,1,1) = [1 1 1 1 1 1 1 1 1];
    phase_mod_all(1,:,1,2) = [1 1 1 exp(-1i*2/3*pi) exp(-1i*2/3*pi) exp(-1i*2/3*pi) exp(-1i*4/3*pi) exp(-1i*4/3*pi) exp(-1i*4/3*pi)];
    phase_mod_all(1,:,1,3) = [1 1 1 exp(-1i*4/3*pi) exp(-1i*4/3*pi) exp(-1i*4/3*pi) exp(-1i*2/3*pi) exp(-1i*2/3*pi) exp(-1i*2/3*pi)];
    phase_mod_all = repmat(phase_mod_all,[1 ceil(sy*sz/9) 1 1]);
    phase_mod_all(:,sy*sz+1:end,:,:) = [];
    case 0
    phase_mod_all = ones([1 sy*sz]);
    case 4
    phase_mod_all(1,:,1,1) = [1 1 1 1 exp(-1i*2/3*pi) exp(-1i*4/3*pi) 1 exp(-1i*4/3*pi) exp(-1i*2/3*pi)];
    phase_mod_all(1,:,1,2) = [1 exp(-1i*2/3*pi) exp(-1i*4/3*pi) 1 exp(-1i*4/3*pi) exp(-1i*2/3*pi) 1 1 1];
    phase_mod_all(1,:,1,3) = [1 exp(-1i*4/3*pi) exp(-1i*2/3*pi) 1 1 1 1 exp(-1i*2/3*pi) exp(-1i*4/3*pi)];
    phase_mod_all = repmat(phase_mod_all,[1 ceil(sy*sz/9) 1 1]);
    phase_mod_all(:,sy*sz+1:end,:,:) = [];
    case 5
        phase_mod_all(1,:,1,1) = ones([1 sy]);
        phase_mod_all(1,:,1,2) = exp(-1i*2/3*pi*(0:sy-1));
        phase_mod_all(1,:,1,3) = exp(-1i*4/3*pi*(0:sy-1));
        
        temp1 = circshift(phase_mod_all,1,2);
        temp2 = circshift(temp1,1,2);
        phase_mod_all = cat(2,phase_mod_all,temp1,temp2);
        phase_mod_all = repmat(phase_mod_all,[1 ceil(sz/3) 1 1]);
        phase_mod_all(:,sz*sy+1:end,:,:) = [];
    case 6

        data_angle  = para.dataAngle/180*pi;
        golden_angle = ((sqrt(5)-1)/2)*pi;
        theta_all = mod((0:sy*sz-1)*golden_angle,data_angle);
        theta_all = reshape(theta_all,[sy,sz]);
        for i=1:sz
            theta_temp = theta_all(:,i);
            [~,I] = sort(theta_temp);
            phase_mod(:,1) = I;
            phase_mod(:,2) = ones([1 sy]);
            phase_mod(:,3) = exp(-1i*2/3*pi*(0:sy-1));
            phase_mod(:,4) = exp(-1i*4/3*pi*(0:sy-1));
            phase_mod_all(:,i,:) = phase_mod(:,2:end);
        end
        phase_mod_all = reshape(phase_mod_all,[1 sy*sz 1 3]);


end

end