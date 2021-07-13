function [isurr,Isurr,Psurr]=surrogate_analysis(Data_surr,iY,iX1,iX2,nfft,fs)

dataS=Data_surr;
%%% PARAMETRIC ESTIMATION
maxIP=30;
[~,IP,~,~] = mos_idVAR(dataS',maxIP,0);
IP=8;% average value along the 91 subjects
%%% parameteric VAR identification
[eAmS,eSuS]=idVAR(dataS',IP,0); %series in rows
%%% cross spectral analysis
outS = fdVAR(eAmS,eSuS,nfft,fs);
fsurr=outS.f; %frequency axis
Psurr=outS.S; %spectral matrix

outSurr=iispectral(Psurr,iY,iX1,iX2);

isurr=outSurr.iy_x1_x2; %interaction information

%%% time domain measures (integral over the whole frequency axis)
Isurr=(2/fs)*trapz(fsurr',isurr);

end