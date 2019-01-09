function bins = get_resporation_bins(Image)
Image = abs(Image);
[sx,sy,nof,set,sms] = size(Image);
signal = reshape(Image,[sx*sy,nof,set,sms]);
%signal(:,1:10,:,:) = [];% delete proton density images
signal = permute(signal,[1,4,2,3]);
signal = reshape(signal,[sx*sy*sms,nof,set]);

signal_mean = mean(signal,2);
signal_local_mean = (signal + circshift(signal,[0,1,0]) + circshift(signal,[0,-1,0]))/3;
signal_local_mean(:,1,:) = (signal(:,1,:) + signal(:,2,:))/2;
signal_local_mean(:,end,:) = (signal(:,end,:) + signal(:,end-1,:))/2;

%signal_cardiac = abs(signal./signal_local_mean);

for i=1:set
    signal_mean_temp = signal_mean(:,:,i);
    signal_res_temp = signal_local_mean(:,:,i);
    [~,order] = sort(signal_mean_temp,'descend');
    N = length(order);
    N = round(N/2);
    signal_res_temp = signal_res_temp(order(1:N),:);
    max_temp = max(signal_res_temp,[],2);
    signal_res_temp = signal_res_temp./max_temp;
    mean_temp = mean(signal_res_temp,2);
    [~,order] = sort(mean_temp,'descend');
    N = round(N/2);
    signal_res_temp = signal_res_temp(order(1:N),:);
    signal_res_all(:,i,:) = signal_res_temp;
end

signal_res_all = reshape(signal_res_all,[N*set,nof]);
signal_res_all = signal_res_all.';
coeff = pca(signal_res_all);
signal_pca = signal_res_all*coeff(:,1:10);
signal_fft = abs(fft(signal_pca,[],1));
freq = round(nof/10):round(nof/5);
[~,idx] = max(max(signal_fft(freq,:),[],1));
signal = signal_pca(:,idx);
signal = signal - min(signal) + 1;

[mins,minlocs] = findpeaks(-signal);
mins = -mins;
signal_smoothed = smooth(signal);
notwanted = signal_smoothed(minlocs)<mins;
mins(notwanted) = [];
minlocs(notwanted) = [];
signal_min = interp1(minlocs,mins,(1:nof).');
%signal_min_nearest = interp1(minlocs,-mins,(1:nof).','nearest');
signal_min(1:minlocs(1)) = min(signal(1:minlocs(1)),signal_min(minlocs(1)));
signal_min(minlocs(end):end) = min(signal(minlocs(end):end),signal_min(minlocs(end)));
%signal_min(signal<signal_min) = signal_min_nearest(signal<signal_min);
flag = signal - signal_min;
flag(flag>0) = 0;

while sum(flag) < 0
    [~,loc_add] = findpeaks(-flag);
    minlocs = [minlocs;loc_add];
    mins = [mins;signal(loc_add)];
    [minlocs,order] = sort(minlocs);
    mins = mins(order);
    signal_min = interp1(minlocs,mins,(1:nof).');
    signal_min(1:minlocs(1)) = min(signal(1:minlocs(1)),signal_min(minlocs(1)));
    signal_min(minlocs(end):end) = min(signal(minlocs(end):end),signal_min(minlocs(end)));
    flag = signal - signal_min;
    flag(flag>0) = 0;
end

[pks,pklocs] = findpeaks(signal);
notwanted = signal_smoothed(pklocs)>pks;
pks(notwanted) = [];
pklocs(notwanted) = [];
signal_max = interp1(pklocs,pks,(1:nof).');
%signal_max_nearest = interp1(pklocs,pks,(1:nof).','nearest');
signal_max(1:pklocs(1)) = max(signal(1:pklocs(1)),signal_max(pklocs(1)));
signal_max(pklocs(end):end) = max(signal(pklocs(end):end),signal_max(pklocs(end)));
%signal_max(signal>signal_max) = signal_max_nearest(signal>signal_max);
flag = signal_max - signal;
flag(flag>0) = 0;
while sum(flag) < 0
    [~,loc_add] = findpeaks(-flag);
    pklocs = [pklocs;loc_add];
    pks = [pks;signal(loc_add)];
    [pklocs,order] = sort(pklocs);
    pks = pks(order);
    signal_max = interp1(pklocs,pks,(1:nof).');
    signal_max(1:pklocs(1)) = max(signal(1:pklocs(1)),signal_max(pklocs(1)));
    signal_max(pklocs(end):end) = max(signal(pklocs(end):end),signal_max(pklocs(end)));
    flag = signal_max - signal;
    flag(flag>0) = 0;
end
signal = (signal - signal_min)./(signal_max - signal_min);
%signal(signal>1) = 1;
%signal = signal-min(signal);
%signal = signal/max(signal);
N = round(nof/4);
[~,order] = sort(signal);
bins = zeros(nof,1);
bins(order(1:N)) = 1;
bins(order(N+1:2*N)) = 2;
bins(order(2*N+1:3*N)) = 3;
bins(order(3*N+1:end)) = 4;
%figure,plot(signal)


