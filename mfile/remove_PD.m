function [Data,para] = remove_PD(Data,para)

PD_frames = para.Recon.PD_frames;
if PD_frames(end) == 0
    PD_frames = 1:sum(PD_frames);
end
Data.first_est(:,:,PD_frames,:,:) = [];
Data.kSpace(:,:,PD_frames,:,:,:,:) = [];
if isfield(Data,'mask')
    Data.mask(:,:,PD_frames,:,:,:,:) = [];
end
if isfield(Data,'phase_mod')
    Data.phase_mod(:,:,PD_frames,:,:,:,:) = [];
end
if isfield(Data,'first_guess')
    Data = rmfield(Data,'first_guess');
end
if isfield(Data,'N')
    NPD = sum(PD_frames);
    Data.N.S = Data.N.S(Data.N.sx_over.^2*NPD+1:end,Data.N.siz(1)*Data.N.siz(2)*NPD+1:end);
    Data.N.siz(3) = size(Data.first_est,3);
end
para.Recon = rmfield(para.Recon,'PD_frames');
%para.Recon = rmfield(para.Recon,'RF_frames');