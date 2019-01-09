function [kSpace_all,para] = PCA_kSpace_unstreaking_yt(para)

disp('Load k space data...');

PCA_dir = para.dir.PCA_dir;
ifPCA = ~isempty(dir(PCA_dir));

if ifPCA
    load(PCA_dir)
    [sx,sy,sz,no_comp,ns] = size(kSpace);
    disp('PCA result alreaty exists, skip PCA...');
else
    load(strcat(para.dir.load_kSpace_dir,para.dir.load_kSpace_name));
    if exist('AIFData','var')
        kSpace = AIFData; clear AIFData
        kSpace = permute(kSpace,[1 2 5 4 3]);
    end
    if(size(kSpace,4)==1)
        kSpace = permute(kSpace,[1 2 3 5 4]);
    end

    [sx,sy,nc,sz,ns] = size(kSpace);
    
    idx = sum(sum(sum(sum(kSpace,1),2),3),4)==0;

    kSpace(:,:,:,:,idx) = [];
    ns = size(kSpace,5);
    
%%%%% PCA on coils

    disp('Peforming PCA on coils...')
    if nc > para.NumberofPCAComponent
        kSpace = coil_unstreaking_and_compression_SMS(kSpace,kSpace_info);
        no_comp = para.NumberofPCAComponent;    
        save(PCA_dir,'kSpace','kSpace_info');
    else
        no_comp = nc;
        kSpace = permute(kSpace,[1 2 4 3]);
    end
    
end

if isfield(para,'slice_pick') && ns~=1
    kSpace = kSpace(:,:,:,:,para.slice_pick);
    if size(kSpace_info.phase_mod,2)~=1
        kSpace_info.phase_mod = kSpace_info.phase_mod(:,para.slice_pick);
        para.phase_mod = para.phase_mod(:,para.slice_pick);
        para.angle_mod = para.angle_mod(:,para.slice_pick);
    end
    ns = length(para.slice_pick);
end

kSpace = phase_correction_031218(kSpace,kSpace_info.phase_mod);
kSpace_all = reshape(kSpace,[sx,sy*sz,no_comp,1,ns]);

para.Recon.sy = sy;
para.Recon.sz = sz;
para.Recon.sx = sx;
para.Recon.no_comp = no_comp;
para.Recon.image_size = [sx,sx];


