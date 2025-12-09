function [moptime,fidx,Ho,Lo,LoSwl,LoSea,Tp,TpSwl,TpSea]=GetMopStockdonParams(MopID,Rstart,Rend,MopQCflag)

% first get the Mop frequency spectra
[freq,bw,depth,moptime,fidx,a0] = GetMOPfreqSpectra(MopID,Rstart,Rend,MopQCflag);

% get unshoaling coefficients and deep water wavelength as a function
%  of the Mop frequencies and depth

% need group speed Cg at Mop depth
[Lf,Cf,Cgf]=LinearDispersion(freq,depth);
% need group speed Cgo and wavelengths Lo in deep water
[Lof,Cof,Cgof]=LinearDispersion(freq,2000);

% unshoal Mop freq spectra to deep water and convert to m^2 energies
a0deep=a0.*repmat((Cgf./Cgof),1,size(a0,2));

% get deep water significant wave height
Ho=4*sqrt(sum(a0deep.*repmat(bw,1,size(a0,2)),'omitnan'));

% find peak energy density frequency band for each spectrum
[Epeak,fpeak]=max(a0deep);
[EpeakSwl,fpeakSwl]=max(a0deep(1:10,:));
[EpeakSea,fpeakSea]=max(a0deep(11:end,:));fpeakSea=fpeakSea+10;

% use peak freq bands to define Tp and Lo for each spectrum
T=1./freq;
Tp=T(fpeak)';
TpSwl=T(fpeakSwl)';
TpSea=T(fpeakSea)';
Lo=Lof(fpeak)';
LoSwl=Lof(fpeakSwl)';
LoSea=Lof(fpeakSea)';

end