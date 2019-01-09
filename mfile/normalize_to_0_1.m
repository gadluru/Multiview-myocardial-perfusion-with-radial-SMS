function normalized_signal = normalize_to_0_1(signal)
signal = signal(:);
N = length(signal);
signal = signal - min(signal);
smoothed_signal = smooth(signal);
[mins,locs] = findpeaks(-signal);
mins = -mins;
idx_notwanted = mins > smoothed_signal(locs);
mins(idx_notwanted) = [];
locs(idx_notwanted) = [];
cali = interp1(locs,mins,(1:N).');
cali(1:locs(1)) = min(signal(1:locs(1)),cali(locs(1)));
cali(locs(end):end) = min(signal(locs(end):end),cali(locs(end)));
cali_temp = signal - cali;
cali_temp(cali_temp>0) = 0;
[mins_0,locs_0] = findpeaks(-signal);
mins_0 = -mins_0;
cali_temp2 = mins_0-cali(locs_0);
j = 0;
while (sum(cali_temp)<0 || sum(cali_temp2)>0) && j<100
    [~,locs_add] = findpeaks(-cali_temp);
    locs = [locs;locs_add];
    mins = [mins;signal(locs_add)];
    [~,locs_add] = findpeaks(cali_temp2);
    locs_add = locs_0(locs_add);
    for i=1:length(locs_add)
        idx = locs==locs_add(i);
        if sum(idx)~=0
            locs(idx) = [];
            mins(idx) = [];
        end
    end
    locs = [locs;locs_add];
    mins = [mins;signal(locs_add)];
    [locs,order] = sort(locs);
    mins = mins(order);
    cali = interp1(locs,mins,(1:N).');
    cali(1:locs(1)) = min(signal(1:locs(1)),cali(locs(1)));
    cali(locs(end):end) = min(signal(locs(end):end),cali(locs(end)));
    cali_temp = signal - cali;
    cali_temp(cali_temp>0) = 0;
    cali_temp2 = mins_0-cali(locs_0);
    j = j+1;
end

cali_min = cali;
        
[maxs,locs] = findpeaks(signal);
idx_notwanted = maxs < smoothed_signal(locs);
maxs(idx_notwanted) = [];
locs(idx_notwanted) = [];
cali = interp1(locs,maxs,(1:N).');
cali(1:locs(1)) = max(signal(1:locs(1)),cali(locs(1)));
cali(locs(end):end) = max(signal(locs(end):end),cali(locs(end)));
cali_temp = cali - signal;
cali_temp(cali_temp>0) = 0;
[maxs_0,locs_0] = findpeaks(signal);
cali_temp2 = cali(locs_0)-maxs_0;
j = 0;
while (sum(cali_temp)<0 || sum(cali_temp2)>0) && j<100
    [~,locs_add] = findpeaks(-cali_temp);
    locs = [locs;locs_add];
    maxs = [maxs;signal(locs_add)];
    [~,locs_add] = findpeaks(cali_temp2);
    locs_add = locs_0(locs_add);
    for i=1:length(locs_add)
        idx = locs==locs_add(i);
        if sum(idx)~=0
            locs(idx) = [];
            maxs(idx) = [];
        end
    end
    locs = [locs;locs_add];
    maxs = [maxs;signal(locs_add)];
    [locs,order] = sort(locs);
    maxs = maxs(order);
    cali = interp1(locs,maxs,(1:N).');
    cali(1:locs(1)) = max(signal(1:locs(1)),cali(locs(1)));
    cali(locs(end):end) = max(signal(locs(end):end),cali(locs(end)));
    cali_temp = cali - signal;
    cali_temp(cali_temp>0) = 0;
    cali_temp2 = cali(locs_0)-maxs_0;
    j = j+1;
end

normalized_signal = (signal-cali_min)./(cali-cali_min);