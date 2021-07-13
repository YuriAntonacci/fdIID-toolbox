% Simulation study with 10 realizations of the VAR model (Eqs. 15, c=0.5)
clear; close all; clc
% Simulation parameters
N_Real=10; %number of realizations
M=4; %number of processes
%%% MVAR process parameters
iY=1; iX1=[2]; iX2=[3 4];% target and sources
%%% set poles and self-oscillations
par.poles{2}=([0.85 0.05]); 
par.poles{3}=([0.85 0.05]); 
par.poles{4}=([0.95 0.35]); 

%% simulation design
%%% common driver effect Z4->Z3 and Z4-> Z2 (c=0.5)
c=0.5;
C_21=c;
C_31=1-c;
par.coup=[2 1 1 C_21;3 1 2 C_31;4 2 1 0.5;4 3 3 0.5];
par.Su=ones(1,M); %variance of innovation processes

[Am,Su,Ak]=theoreticalVAR(M,par); %% Theoretical VAR parameters

Q=size(Am,2);%number of time series
p=size(Am,1)/Q;%model order
N=1000;%numer of time points (K=30)

for it=1:N_Real
    disp(it)
    U=randn(Q,N); %uncorrelated noise
    % Create Simulated Time-series with predefined structure
    [Y]=MVARfilter(Am',U);
   
    Proc.samp=Y;
    Proc.popt=p;
    Proc.model=Am';
    
    fs=100; % sampling frequency
    nfft=1001; %number of points on frequency axis (total)
    M=4;
    iY=1; iX1=[2]; iX2=[3 4];
    
    eAm=Proc.model;
    eSu=eye(M);
    outf = fdVAR(eAm,eSu,nfft,fs);
    f=outf.f; %frequency axis
    P_theo=outf.S; %spectral matrix
    
    %%% theoretical analysis
    outT=iispectral(P_theo,iY,iX1,iX2);
    %%% Interaction Information Decomposition
    fy_x1T(:,it)=outT.fy_x1;
    fy_x2T(:,it)=outT.fy_x2;
    fy_x1x2T(:,it)=outT.fy_x1x2;
    iy_x1_x2T(:,it)=outT.iy_x1_x2; 
    
    %%% time domain measures (integral over the whole frequency axis)
    Fy_x1T(it,1)=(2/fs)*trapz(f',fy_x1T(:,it));
    Fy_x2T(it,1)=(2/fs)*trapz(f',fy_x2T(:,it));
    Fy_x1x2T(it,1)=(2/fs)*trapz(f',fy_x1x2T(:,it));
    Iy_x1_x2T(it,1)=(2/fs)*trapz(f',iy_x1_x2T(:,it));
    
    %%% PARAMETRIC ESTIMATION
    maxIP=10;
    data=Y;
    % Testing for restricted form of weak-sense stationarity
    for ch=1:size(data,1)
        [~, mean_stat_flag(ch), var_stat_flag(ch),covtest] = isstationary(Y(ch,:)');
    end
    
    %%% MDL criterion for VAR model order (IP) estimation
    maxIP=30;
    [~,IP,~,~] = mos_idVAR(Y,maxIP,0);
    %%% parameteric VAR identification
    [eAm,eSu,~,Residuals]=idVAR(data,IP,0); %series in rows
    %%% Durbin-Watson test for whiteness (H0: no serial correlation) of VAR residuals
    [dw,pval] = whiteness(data,Residuals);
    %%% cross spectral analysis
    outf = fdVAR(eAm,eSu,nfft,fs);
    f=outf.f; %frequency axis
    P_AR1=outf.S; %spectral matrix
    
    %IID - Parametric approach (Section III.B)
    outAR=iispectral(P_AR1,iY,iX1,iX2);
    
    fy_x1P(:,it)=outAR.fy_x1;
    fy_x2P(:,it)=outAR.fy_x2;
    fy_x1x2P(:,it)=outAR.fy_x1x2;
    iy_x1_x2P(:,it)=outAR.iy_x1_x2; 
    
    Fy_x1P(it,1)=(2/fs)*trapz(f',fy_x1P(:,it));
    Fy_x2P(it,1)=(2/fs)*trapz(f',fy_x2P(:,it));
    Fy_x1x2P(it,1)=(2/fs)*trapz(f',fy_x1x2P(:,it));
    Iy_x1_x2P(it,1)=(2/fs)*trapz(f',iy_x1_x2P(:,it));
    
    
    % Non parametric approach (Blackman-Tuckey analysis)
    P=zeros(M,M,nfft);
    m=round((1.273*fs)/(2.5)); % number of truncation point for cross-correlation analysis (see reference 36)

    for ch=1:M
        x=data(ch,:)';
        x=x-mean(x);
        [rx, lags]=xcorr(x,m,'biased');
        nfft2=2*nfft;
        wi=parzenwin(2*m+1); % Parzen window 
        rxw=rx.*wi;
        px=fft(rxw,nfft2);% Fourier Transform of cross-correlation function
        px=(2/fs)*px(1:nfft); 
        f=(0:fs/(2*(nfft-1)):fs/2)'; %frequencies vector
        freq = f;
        for cc=1:M
            if cc~=ch
                y=data(cc,:)';
                y=y-mean(y);
                [rxy, lagxy]=xcorr(x,y,m,'biased');
                nfft2=2*nfft;
                wixy=parzenwin(2*m+1);
                rxyw=rxy.*wixy;
                pxy=fft(rxyw,nfft2);
                pxy=(2/fs)*pxy(1:nfft); 
                f=(0:fs/(2*(nfft-1)):fs/2)';
                freq = f;
                P(ch,cc,:)=pxy;
            end
            
        end
        P(ch,ch,:) = px;
    end
    
    outNP=iispectral(P,iY,iX1,iX2);
    %%% IID with non parametric approach
    fy_x1NP(:,it)=outNP.fy_x1;
    fy_x2NP(:,it)=outNP.fy_x2;
    fy_x1x2NP(:,it)=outNP.fy_x1x2;
    iy_x1_x2NP(:,it)=outNP.iy_x1_x2; %interaction information
     
    %%% time domain measures (integral over the whole frequency axis)
    Fy_x1NP(it,1)=(2/fs)*trapz(f',fy_x1NP(:,it));
    Fy_x2NP(it,1)=(2/fs)*trapz(f',fy_x2NP(:,it));
    Fy_x1x2NP(it,1)=(2/fs)*trapz(f',fy_x1x2NP(:,it));
    Iy_x1_x2NP(it,1)=(2/fs)*trapz(f',iy_x1_x2NP(:,it));
    
    
end

%% plot Figure 3 (panels a-d)
%%% median and percentiles along the iterations of the simulation study

lo=5; hi=95; % percentiles for showing distributions
%%% Parametric approach
fy_x1x2_mP=squeeze(median(fy_x1x2P,2));
fy_x1x2_loP=squeeze(prctile(fy_x1x2P,lo,2));
fy_x1x2_hiP=squeeze(prctile(fy_x1x2P,hi,2));
fy_x1x2_theo=squeeze(median(fy_x1x2T,2));

fy_x1_mP=squeeze(median(fy_x1P,2));
fy_x1_loP=squeeze(prctile(fy_x1P,lo,2));
fy_x1_hiP=squeeze(prctile(fy_x1P,hi,2));
fy_x1_theo=squeeze(median(fy_x1T,2));

fy_x2_mP=squeeze(median(fy_x2P,2));
fy_x2_loP=squeeze(prctile(fy_x2P,lo,2));
fy_x2_hiP=squeeze(prctile(fy_x2P,hi,2));
fy_x2_theo=squeeze(median(fy_x2T,2));

iy_x1_x2_mP=squeeze(median(iy_x1_x2P,2));
iy_x1_x2_loP=squeeze(prctile(iy_x1_x2P,lo,2));
iy_x1_x2_hiP=squeeze(prctile(iy_x1_x2P,hi,2));
iy_x1_x2_theo=squeeze(median(iy_x1_x2T,2));

%% Non-parametric approach
fy_x1x2_mNP=squeeze(median(fy_x1x2NP,2));
fy_x1x2_loNP=squeeze(prctile(fy_x1x2NP,lo,2));
fy_x1x2_hiNP=squeeze(prctile(fy_x1x2NP,hi,2));

fy_x1_mNP=squeeze(median(fy_x1NP,2));
fy_x1_loNP=squeeze(prctile(fy_x1NP,lo,2));
fy_x1_hiNP=squeeze(prctile(fy_x1NP,hi,2));

fy_x2_mNP=squeeze(median(fy_x2NP,2));
fy_x2_loNP=squeeze(prctile(fy_x2NP,lo,2));
fy_x2_hiNP=squeeze(prctile(fy_x2NP,hi,2));

iy_x1_x2_mNP=squeeze(median(iy_x1_x2NP,2));
iy_x1_x2_loNP=squeeze(prctile(iy_x1_x2NP,lo,2));
iy_x1_x2_hiNP=squeeze(prctile(iy_x1_x2NP,hi,2));

freq=f;

%% Plot of fy_x1x2
HCF=figure('units','inches','position',[0 0 11.7 8.3]);
orient(HCF,'landscape')
[ha, pos] = tight_subplot(1,4, 0.05,0.05,0.05);
hold on
set(0,'defaultAxesFontSize', 8);
set(0,'defaultTextFontSize', 8);

colth=[128 0 0]/255; % color of theoretical value
colm=[86 86 158]/255; % color of median value
colm_np=[0.8500, 0.3250, 0.0980];
colsh=[189 185 219]/255; %color of shades
colsh_np=[0.9290, 0.6940, 0.1250];
axes(ha(1));
set(ha(1:4),'XTickLabel',{'0','10','20','30','40','50'})
set(ha(1),'YTickLabel',{'0','1','2','3','4','5'})

h1=area(freq,[fy_x1x2_loP fy_x1x2_mP-fy_x1x2_loP fy_x1x2_hiP-fy_x1x2_mP]); hold on
set(h1(1),'FaceColor','w','FaceAlpha',0); set(h1(2),'FaceColor',colsh,'FaceAlpha',0); set(h1(3),'FaceColor',colsh,'FaceAlpha',0);
set(h1(1),'EdgeColor','k'); set(h1(2),'EdgeColor',colsh); set(h1(3),'EdgeColor',colsh); 
set(h1(2),'FaceAlpha',0.5);
set(h1(3),'FaceAlpha',0.5);
hold on
h2=area(freq,[fy_x1x2_loNP fy_x1x2_mNP-fy_x1x2_loNP fy_x1x2_hiNP-fy_x1x2_mNP]); hold on
set(h2(1),'FaceColor','w','FaceAlpha',0); set(h2(2),'FaceColor',colsh_np,'FaceAlpha',0); set(h2(3),'FaceColor',colsh_np,'FaceAlpha',0);
set(h2(1),'EdgeColor','k'); set(h2(2),'EdgeColor',colsh_np); set(h2(3),'EdgeColor',colsh_np);% set(h2(3),'LineWidth',1.5);
p3=plot(freq,fy_x1x2_mNP,'Color',colm_np,'LineWidth',1.5);
set(h2(2),'FaceAlpha',0.5);
set(h2(3),'FaceAlpha',0.5);
hold on
p1=plot(freq,fy_x1x2_mP,'b','LineWidth',1.5);
hold on
p2=plot(freq,fy_x1x2_theo,'r','LineWidth',1.5);
title('f_{x_1x_2;y}')

%% Plot of fx1_y
axes(ha(2));
set(ha(2),'YTickLabel',{'0','0.5','1','1.5','2','2.5'})

h3=area(freq,[fy_x1_loP fy_x1_mP-fy_x1_loP fy_x1_hiP-fy_x1_mP]); hold on
set(h3(1),'FaceColor','w','FaceAlpha',0); set(h3(2),'FaceColor',colsh,'FaceAlpha',0); set(h3(3),'FaceColor',colsh,'FaceAlpha',0);
set(h3(1),'EdgeColor','k'); set(h3(2),'EdgeColor',colsh); set(h3(3),'EdgeColor',colsh); 
set(h3(2),'FaceAlpha',0.5);
set(h3(3),'FaceAlpha',0.5);
hold on
h4=area(freq,[fy_x1_loNP fy_x1_mNP-fy_x1_loNP fy_x1_hiNP-fy_x1_mNP]); hold on
set(h4(1),'FaceColor','w','FaceAlpha',0); set(h4(2),'FaceColor',colsh_np,'FaceAlpha',0); set(h4(3),'FaceColor',colsh_np,'FaceAlpha',0);
set(h4(1),'EdgeColor','k'); set(h4(2),'EdgeColor',colsh_np); set(h4(3),'EdgeColor',colsh_np); 
p3=plot(freq,fy_x1_mNP,'Color',colm_np,'LineWidth',1.5);
set(h4(2),'FaceAlpha',0.5);
set(h4(3),'FaceAlpha',0.5);
hold on
p1=plot(freq,fy_x1_mP,'b','LineWidth',1.5);
hold on
p2=plot(freq,fy_x1_theo,'r','LineWidth',1.5);
title('f_{x_1;y}')

%% Plot of fx2_y
axes(ha(3));
set(ha(3),'YTickLabel',{'0','0.5','1','1.5','2','2.5'})

h5=area(freq,[fy_x2_loP fy_x2_mP-fy_x2_loP fy_x2_hiP-fy_x2_mP]); hold on
set(h5(1),'FaceColor','w','FaceAlpha',0); set(h5(2),'FaceColor',colsh,'FaceAlpha',0); set(h5(3),'FaceColor',colsh,'FaceAlpha',0);
set(h5(1),'EdgeColor','k'); set(h5(2),'EdgeColor',colsh); set(h5(3),'EdgeColor',colsh);
set(h5(2),'FaceAlpha',0.5);
set(h5(3),'FaceAlpha',0.5);
hold on
h6=area(freq,[fy_x1_loNP fy_x1_mNP-fy_x1_loNP fy_x1_hiNP-fy_x1_mNP]); hold on
set(h6(1),'FaceColor','w','FaceAlpha',0); set(h6(2),'FaceColor',colsh_np,'FaceAlpha',0); set(h6(3),'FaceColor',colsh_np,'FaceAlpha',0);
set(h6(1),'EdgeColor','k'); set(h6(2),'EdgeColor',colsh_np); set(h6(3),'EdgeColor',colsh_np);% set(h6(3),'LineWidth',1.5);
p3=plot(freq,fy_x1_mNP,'Color',colm_np,'LineWidth',1.5);
set(h6(2),'FaceAlpha',0.5);
set(h6(3),'FaceAlpha',0.5);
hold on
p1=plot(freq,fy_x2_mP,'b','LineWidth',1.5);
hold on
p2=plot(freq,fy_x2_theo,'r','LineWidth',1.5);
title('f_{x_2;y}')
%% Plot of i_x1_x2_y

axes(ha(4));
set(ha(4),'YLim',[-4 3])
set(ha(4),'YTickLabel',{'-4','-3','-2','-1','0','1','2','3'})
h7=area(freq,[iy_x1_x2_loP iy_x1_x2_mP-iy_x1_x2_loP iy_x1_x2_hiP-iy_x1_x2_mP]); hold on
set(h7(1),'FaceColor','w','FaceAlpha',0); set(h7(2),'FaceColor',colsh,'FaceAlpha',0); set(h7(3),'FaceColor',colsh,'FaceAlpha',0);
set(h7(1),'EdgeColor','k'); set(h7(2),'EdgeColor',colsh); set(h7(3),'EdgeColor',colsh);
set(h7(2),'FaceAlpha',0.5);
set(h7(3),'FaceAlpha',0.5);
hold on
h8=area(freq,[iy_x1_x2_loNP iy_x1_x2_mNP-iy_x1_x2_loNP iy_x1_x2_hiNP-iy_x1_x2_mNP]); hold on
set(h8(1),'FaceColor','w','FaceAlpha',0); set(h8(2),'FaceColor',colsh_np,'FaceAlpha',0); set(h8(3),'FaceColor',colsh_np,'FaceAlpha',0);
set(h8(1),'EdgeColor','k'); set(h8(2),'EdgeColor',colsh_np); set(h8(3),'EdgeColor',colsh_np); %set(h8(3),'LineWidth',1.5);
p9=plot(freq,iy_x1_x2_mNP,'Color',colm_np,'LineWidth',1.5);
set(h8(2),'FaceAlpha',0.5);
set(h8(3),'FaceAlpha',0.5);
hold on
p10=plot(freq,iy_x1_x2_mP(:,1),'b','LineWidth',1.5);
hold on
p11=plot(freq,iy_x1_x2_theo,'r','LineWidth',1.5);
title('i_{x_1;x_2;y}')

