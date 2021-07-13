%%% Information decomposition of time series recorded from a ring of 32 
%%% chaotic oscillators (Section IV.B -main document)

clear all
close all
clc

load('top_row.mat');% load data

data=double(ts);
fs=100000; % 100 Khz sampling frequency

%%% decimation of time series
for tt=1:size(data,1)
    data1(tt,:) = decimate(data(tt,:)',10,'fir');
end
fs_eq=fs/10; % equivalent sampling frequency

%%%Definition of the Target and the Source 1 
T=circularBuffer([1:32]');% introducing T as circular buffer
S=[3:32,1,2]';

data_d=data1; 
maxIP=30;
Nsurr=10; % number of surrogates (100 in the main document)
%%% with this settings are necessary 400 seconds for each target and ten
%%% surrogates

for tar=1%:length(T) selection of the first target (to reduce computational time)
    disp(T(tar))
    tic
    W=T; % W it is equal to Target vector
    W([tar-3:tar+3])=[]; % delete the elements with a distance <=4 from T(tar)
    for ll=1:length(W)
        disp(ll)
        DATA=data_d([T(tar) S(tar) W(ll)],:);
        % Testing for restricted form of weak-sense stationarity
        for ch=1:size(DATA,1)
            [~, mean_stat_flag(ch), var_stat_flag(ch)] = isstationary(squeeze(DATA(ch,:)'));
        end
        
        [~,IP,~,~] = mos_idVAR(DATA,maxIP,0);
        %         IP=17; %% average value along the entire ring
        [eAm,eSu,~,Residuals]=idVAR(DATA,IP,0);
        %%% Durbin-Watson test for whiteness (H0: no serial correlation) of VAR residuals
        [dw,pval] = whiteness(DATA,Residuals);
        
        COUPLES(ll,:,tar)=[T(tar) S(tar) W(ll)];
        %         DISTANCES(ll,tar)=abs(T(tar)-W);
        p=nextpow2(size(data_d,2));
        nfft=2^p;
        % Parameteric estimation of Power Spectral Density
        outvar_est = fdVAR(eAm,eSu,nfft,fs_eq);
        eP1=outvar_est.S; %spectral matrix
        f=outvar_est.f; % vector of frequencies
        %%% set target and sources series
        iY=1;% T(i)
        iX1=[2];% S(l)
        iX2=[3]; % W(m)
        
        % frequency and time domain source interaction measures
        out_t=iispectral(eP1,iY,iX1,iX2); %
        iy_x1_x2(:,ll,tar)=out_t.iy_x1_x2;
        fy_x1(:,ll,tar)=out_t.fy_x1;
        fy_x2(:,ll,tar)=out_t.fy_x2;
        fy_x1x2(:,ll,tar)=out_t.fy_x1x2;
        Iy_x1_x2(ll,tar)=out_t.Iy_x1_x2;
        Fy_x1(ll,tar)=out_t.Fy_x1;
        Fy_x2(ll,tar)=out_t.Fy_x2;
        Fy_x1x2(ll,tar)=out_t.Fy_x1x2;
        %%% surrogate analysis
        DATAs=DATA';
        for ns=1:Nsurr
            [X1surr,X2surr]=surriaafft2(DATAs,iY,iX1,iX2);
            surr=[DATAs(:,[iY]),X1surr,X2surr];
            Data_surr(:,:,ns)=surr;
        end
        
        
        for NN=1:Nsurr
            dataS=Data_surr(:,:,NN);
            [~,IP,~,~] = mos_idVAR(dataS',maxIP,0);
            %%% parameteric VAR identification
            [eAmS,eSuS]=idVAR(dataS',IP,0); %series in rows
            %%% cross spectral analysis
            outS = fdVAR(eAmS,eSuS,nfft,fs_eq);
            fsurr=outS.f; %frequency axis
            Psurr(:,:,:,NN)=outS.S; %spectral matrix
            
            outSurr=iispectral(Psurr(:,:,:,NN),iY,iX1,iX2);
            iy_x1_x2S(:,NN,ll,tar)=outSurr.iy_x1_x2; %interaction information
            %%% time domain measures (integral over the whole frequency axis)
            Iy_x1_x2S(NN,ll,tar)=(2/fs)*trapz(f',iy_x1_x2S(:,NN));
            
        end
       
    end
    toc
    
    
end

%%% evaluate the distance of the source X2=[4:29] from the Target 1
DISTANCE=[4:16,15:-1:4]';

%% Plot of the results (as in Figure 7 - main document)
dist_sel=[1,5,9,13]; % selection of different distances (4,8,12,16)
colors=[0 0 1;0 0 0;0 1 0;1 0 0];
% average value along the target for time and frequency domain measures
Fy_x1x2=mean(Fy_x1x2,2); 
Fy_x1=mean(Fy_x1,2);
Fy_x2=mean(Fy_x2,2);
Iy_x1_x2=mean(Iy_x1_x2,2);
Iy_x1_x2_surr=mean(Iy_x1_x2S,3);

fy_x1x2=mean(fy_x1x2,3);
fy_x1=mean(fy_x1,3);
fy_x2=mean(fy_x2,3);
iy_x1_x2=mean(iy_x1_x2,3);
iy_x1_x2_surr=mean(iy_x1_x2S,4); 
i_mean=iy_x1_x2;

%%% Figure of Frequnecy domain Analysis
Measures={'f_{x_1;y}','f_{x_2;y}','f_{x_1x_2;y}','i_{x_1;x_2;y}'};
HCF=figure('units','inches','position',[0 0 11.7 8.3]);
orient(HCF,'landscape')
[ha, pos] = tight_subplot(3,2, 0.05,0.05,0.05);
hold on
%%% IID frequency domain
set(0,'defaultAxesFontSize', 8);
set(0,'defaultTextFontSize', 8);
axes(ha(1));
for tt=1:4
    plot(f,fy_x1(:,[dist_sel(tt)]),'Color',colors(tt,:),'LineWidth',1);
    hold on
end
title(Measures{1});
ylim([-1 10])
axes(ha(2));
for tt=1:4
    plot(f,fy_x2(:,[dist_sel(tt)]),'Color',colors(tt,:)','LineWidth',1);
    hold on
end
title(Measures{2});
ylim([-0.6 8])
axes(ha(3));
for tt=1:4
    plot(f,fy_x1x2(:,[dist_sel(tt)]),'Color',colors(tt,:),'LineWidth',1);
    hold on
end
title(Measures{3});
ylim([-1 10])
axes(ha(4));
for tt=1:4
    plot(f,iy_x1_x2(:,[dist_sel(tt)]),'Color',colors(tt,:),'LineWidth',1);
    hold on
end
title(Measures{4});
ylim([-3 8])
set(ha(1:5),'XTickLabel',{'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5'})

%%% Interaction information measure with surrogates analysis
axes(ha(5));
DIST=1; %% only distance of X2=4 from the target (Figure 7e)
lo=2.5; hi=97.5; % percentiles for showing distributions

colsh=[0.9290, 0.6940, 0.1250];
colm=[0.8500, 0.3250, 0.0980];
iy_x1_x2_surr_m=mean(iy_x1_x2_surr(:,:,DIST),2);
iy_x1_x2_lo_surr=squeeze(prctile(iy_x1_x2_surr(:,:,[DIST]),lo,2));
iy_x1_x2_hi_surr=squeeze(prctile(iy_x1_x2_surr(:,:,[DIST]),hi,2));
iy_x1_x2_m=squeeze(iy_x1_x2(:,[DIST]));

h7=area(f,[iy_x1_x2_lo_surr iy_x1_x2_surr_m-iy_x1_x2_lo_surr iy_x1_x2_hi_surr-iy_x1_x2_surr_m]); hold on
set(h7(1),'FaceColor','w','FaceAlpha',0); set(h7(2),'FaceColor',colsh,'FaceAlpha',0); set(h7(3),'FaceColor',colsh,'FaceAlpha',0);
set(h7(1),'EdgeColor','k'); set(h7(2),'EdgeColor',colsh); set(h7(3),'EdgeColor',colsh);% set(h7(3),'LineWidth',1.5);
plot(f,iy_x1_x2_surr_m,'Color',colm,'LineWidth',1);
set(h7(2),'FaceAlpha',0.5);
set(h7(3),'FaceAlpha',0.5);
hold on
plot(f,iy_x1_x2_m,'b','LineWidth',1);
ylim([-3 7])
%% plot of Integrated measures (Time Domain - Figure 7f) 
axes(ha(6))
distance=[1:25];
plot(distance,Fy_x1,'-xr','Linewidth',1,'MarkerSize',3);
hold on;
plot(distance,Fy_x2,'-xb','Linewidth',1,'MarkerSize',3);
hold on;
plot(distance,Fy_x1x2,'-xg','Linewidth',1,'MarkerSize',3);
hold on;
plot(distance,Iy_x1_x2,'-xk','Linewidth',1,'MarkerSize',3);
hold on
ylim([0, 2])
set(ha(6),'XTick',distance);
set(ha(6),'XTickLabels',{'4','5','6','7','8','9','10','11','12','13','14','15','16','15','14','13','12','11','10','9','8','7','6','5','4'});
TIME={'F_{x_1;y}','F_{x_2;y}','F_{x_1x_2;y}','I_{x_1;x_2;y}'};
legend(TIME)

