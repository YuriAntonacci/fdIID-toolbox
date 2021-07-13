clear all
close all
clc

% Parameters
fs = 1; % Sampling frequency (1 month)
nfft = 501; % Spectrum samples
freq = 0:pi/(nfft-1):pi; % Frequency vector (from 0 to pi)
w_real = freq/(2*pi)*fs; % Frequency vector 
bw = 0.05; % Bandwidth
surrogate_num = 100; % Number of surrogate time series
tau = round((1.273*fs)/bw); % Autocorrelation lag 

% Load the data
load('ClimateData.mat');
x = data{:,:};

% Select the channels
iY = 1; % Target
iX1 = [2]; % Driver block 1
iX2 = [3]; % Driver block 2
i_chan = [1,2,3]; % All channels to use
X = x(:,i_chan);
[N,M]=size(X);

% Testing for restricted form of weak-sense stationarity
for ch=1:size(X,2)
    [~, mean_stat_flag(ch), var_stat_flag(ch)] = isstationary(squeeze(X(:,ch)));
end

% Create surrogate time series
surrogates = zeros(N,M,surrogate_num);

for ns = 1:surrogate_num
    [X1surr,X2surr] = surriaafft2(X,10,iX1,iX2);
    surr= [X1surr,X2surr,X(:,[iY])];
    surrogates(:,:,ns)=surr;
end

% Result variables
i_WC = zeros(nfft,1);
i_WC_surr = zeros(nfft,surrogate_num);
P = zeros(M,M,nfft);

% Evaluate the measure for the data and the surrogate time series
for n = 1:surrogate_num+1
    
    if n == 1   % Original time series
        x = X;
    else        % Surrogates
        x = squeeze(surrogates(:,:,n-1));
    end

    % Valuto lo spettro in maniera non parametrica    
    
    % Ciclo tra tutte le serie
    for i = 1:M

        % Evaluate autocorrelation
        [rx, lags] = xcorr(x(:,i),tau,'biased');
        
        nfft2=2*nfft;
        wi=parzenwin(2*tau+1);
        rxw=rx.*wi;
        px=fft(rxw,nfft2);
        px=(2/fs)*px(1:nfft);
        
        % Diagonal elements of the PSD
        P(i,i,:) = px;
        
        % Evaluate cross-correlation
        for j = 1:M
            if i ~= j
                
                [rxy, lagxy]=xcorr(x(:,i),x(:,j),tau,'biased');
                
                nfft2=2*nfft;
                wixy=parzenwin(2*tau+1);
                rxyw=rxy.*wixy;
                pxy=fft(rxyw,nfft2);
                pxy=(2/fs)*pxy(1:nfft); 
                
                % Off-diagonal elements of the PSD
                P(i,j,:)=pxy;
            end
        end 
    end
    
    % Evaluate the spectral interaction information
    outNP = iispectral(P,iY,iX1,iX2);
    iy_x1_x2_WC = outNP.iy_x1_x2;
    
    % Store the results
    if n == 1
        i_WC = iy_x1_x2_WC;
        fy_x1x2 = outNP.fy_x1x2;
        fy_x1 = outNP.fy_x1;
        fy_x2 = outNP.fy_x2;
    else
        i_WC_surr(:,n-1) = iy_x1_x2_WC;       
    end
end

% i_WC contains the spectral decomposition of the interaction information
% i_WC_sull contains the same measure for each surrogate
% Use w_real as the X-axis vector
%% Plot of the results (Figure 8)
Measures={'f_{x_1;y}','f_{x_2;y}','f_{x_1x_2;y}'};
HCF=figure('units','inches','position',[0 0 11.7 8.3]);
orient(HCF,'landscape')
[ha, pos] = tight_subplot(1,2, 0.05,0.05,0.05);
hold on
colors=[0 0 1;0 0 0;0 1 0;1 0 0];
%%% IID frequency domain
set(0,'defaultAxesFontSize', 8);
set(0,'defaultTextFontSize', 8);
axes(ha(1));
IND=301;
plot(w_real(1:IND),fy_x1(1:IND),'Color',colors(1,:),'LineWidth',1.5);
hold on
plot(w_real(1:IND),fy_x2(1:IND),'Color',colors(2,:),'LineWidth',1.5);
hold on
plot(w_real(1:IND),fy_x1x2(1:IND),'Color',colors(3,:),'LineWidth',1.5);
ylim([-0.1 1])
legend(Measures)
%%% Interaction information measure with surrogates analysis
axes(ha(2));
lo=2.5; hi=97.5; % percentiles for showing distributions
colsh=[0.9290, 0.6940, 0.1250];
colm=[0.8500, 0.3250, 0.0980];
iy_x1_x2_surr_m=mean(i_WC_surr((1:IND),:),2);
iy_x1_x2_lo_surr=squeeze(prctile(i_WC_surr((1:IND),:),lo,2));
iy_x1_x2_hi_surr=squeeze(prctile(i_WC_surr((1:IND),:),hi,2));
iy_x1_x2_m=squeeze(i_WC((1:IND)));

h7=area(w_real(1:IND),[iy_x1_x2_lo_surr iy_x1_x2_surr_m-iy_x1_x2_lo_surr iy_x1_x2_hi_surr-iy_x1_x2_surr_m]); hold on
set(h7(1),'FaceColor','w','FaceAlpha',0); set(h7(2),'FaceColor',colsh,'FaceAlpha',0); set(h7(3),'FaceColor',colsh,'FaceAlpha',0);
set(h7(1),'EdgeColor','k'); set(h7(2),'EdgeColor',colsh); set(h7(3),'EdgeColor',colsh);
plot(w_real(1:IND),iy_x1_x2_surr_m,'Color',colm,'LineWidth',1.5);
set(h7(2),'FaceAlpha',0.5);
set(h7(3),'FaceAlpha',0.5);
hold on
plot(w_real(1:IND),iy_x1_x2_m,'b','LineWidth',1.5);
ylim([-0.1 1])
set(ha(2),'XTickLabel',{'0','0.05','0.1','0.15','0.2','0.25','0.3'});
legend('i_{x_1;x_2;y}')