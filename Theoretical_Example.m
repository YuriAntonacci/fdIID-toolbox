%% theoretical simulation with 4 processes, showing coexistence of synergy and redundancy
clear; close all; clc
%%% parameters
fs=100; % sampling frequency
nfft=1001; %number of points on frequency axis (total)

cpar=[0:0.02:1]; %% coupling parameter

%% Simulation - Theoretical VAR process - Equation 15
M=4; %%% number of processes
iY=1; iX1=[2]; iX2=[3 4]; %% target and sources
%%% set poles and self-oscillations
par.poles{2}=([0.85 0.05]); 
par.poles{3}=([0.85 0.05]); 
par.poles{4}=([0.95 0.35]); 

for cnt=1:length(cpar)
    %%% set the interactions as in Figure 1
    c=cpar(cnt);
    C_21=c;
    C_31=1-c;
    par.coup=[2 1 1 C_21;3 1 2 C_31;4 2 1 0.5;4 3 3 0.5];
    par.Su=ones(1,M); %variance of innovation processes
 
    [Am,Su,Ak]=theoreticalVAR(M,par); %% VAR parameters

    % cross-spectral analysis: compute spectral matrices from VAR parameters
    outfd = fdVAR(Am',Su',nfft,fs);
    P=outfd.S; % spectral matrix
    f=outfd.f; % frequency vector

    % frequency domain dynamic source interaction measures 
    out_fd=iispectral(P,iY,iX1,iX2);
    
    iy_x1_x2(:,cnt)=out_fd.iy_x1_x2;
    fy_x1(:,cnt)=out_fd.fy_x1; 
    fy_x2(:,cnt)=out_fd.fy_x2; 
    fy_x1x2(:,cnt)=out_fd.fy_x1x2;
    Iy_x1_x2(cnt)=out_fd.Iy_x1_x2;
    Fy_x1(cnt)=out_fd.Fy_x1; 
    Fy_x2(cnt)=out_fd.Fy_x2; 
    Fy_x1x2(cnt)=out_fd.Fy_x1x2;
   
end


%% plot IID - Figure 2
set(0,'defaultAxesFontSize', 8);
set(0,'defaultTextFontSize', 8);
set(0,'defaultTextFontName', 'Times');
set(0,'defaultAxesFontName', 'Times');

colrange=[];
lim1=floor(length(cpar)/2);
for i=1:lim1, colrange=[colrange; 0 (i-1)/lim1 1-(i-1)/lim1]; end
for i=1:length(cpar)-lim1, colrange=[colrange; (i-1)/lim1 1-(i-1)/lim1 0]; end

ymax1=2; ymax2=6;

HCF=figure('units','inches','position',[0 0 11.7 8.3]);
orient(HCF,'landscape')
[ha, pos] = tight_subplot(1,5, 0.05,0.05,0.05);

for cnt=1:length(cpar)
    axes(ha(2));
    pfig(cnt)=plot(f,fy_x1x2(:,cnt),'color',colrange(cnt,:),'Linewidth',1.1); hold on
    xlim([0 fs/2]); ylim([min(min(fy_x1x2)) max(max(fy_x1x2))]);  xlabel('f') ; title('f_{x_1x_2;y}')

    axes(ha(3));
    pfig(cnt)=plot(f,fy_x1(:,cnt),'color',colrange(cnt,:),'Linewidth',1.1); hold on
    xlim([0 fs/2]); ylim([min(min(fy_x1)) max(max(fy_x1))]);  xlabel('f') ;title('f_{x_1;y}')

    axes(ha(4));
    pfig(cnt)=plot(f,fy_x2(:,cnt),'color',colrange(cnt,:),'Linewidth',1.1); hold on
    xlim([0 fs/2]); ylim([min(min(fy_x2)) max(max(fy_x2))]);   xlabel('f') ;title('f_{x_2;y}')

    axes(ha(5));
    pfig(cnt)=plot(f,iy_x1_x2(:,cnt),'color',colrange(cnt,:),'Linewidth',1.1); hold on
    xlim([0 fs/2]); ylim([min(min(iy_x1_x2)) max(max(iy_x1_x2))]);  xlabel('f'); title('i_{x_1;x_2;y}')
   
end

%%% Plot of integrated measures in the time domain (panel a Figure 2)
axes(ha(1));
legend([pfig(1) pfig(lim1+1) pfig(length(cpar))],['c=' num2str(cpar(1))],['c=' num2str(cpar(lim1+1))],['c=' num2str(cpar(length(cpar)))])
axes(ha(1));
plot(cpar,Fy_x1,'-xr','Linewidth',1.1,'MarkerSize',4);
hold on;
plot(cpar,Fy_x2,'-xb','Linewidth',1.1,'MarkerSize',4);
hold on;
plot(cpar,Fy_x1x2,'-xg','Linewidth',1.1,'MarkerSize',4);
hold on;
plot(cpar,Iy_x1_x2,'-xk','Linewidth',1.1,'MarkerSize',4);
xlabel('c')
hold on
legend('F_{x_1;1}','F_{x_2;y}','F_{x_1x_2;y}','I_{x_1;x_2;y}','Location','southwest');
ylim([-1.5 2])
hold on






