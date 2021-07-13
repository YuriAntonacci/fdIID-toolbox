%%% Analysis of real EEG data recorded during Motor task execution (Section IV.A -main document)
%%% Information about electrodes placement and expertimental design can be
%%% found at: https://physionet.org/content/eegmmidb/1.0.0/
clear all
close all
clc
%%% set as current directory: fd_IID/DataSets/EEG_Analysis

datadir=cd;
subjects=1;
listing=dir(datadir);
sessions={'03','07','11'}; %only REST and close Right fist

duration=160*4; %% duration of each trial (4s)
SUB=1:subjects; 
Fs=160; % sampling frequency

%%% load labels of the EEG channels
load('label.mat')

%%% Define X1, X2 and the Target Y
x1=[5,6,7];% source1=[FC2, FC4, FC6]
x2=[12,13,14];% source 2=[C2, C4, C6]
y=[8 9 10]; % target= [C1, C3, C5]
iter_surr=10; % in the main document it is equal to 100
for ss=1:length(subjects)
    tic
    fold=sprintf('%s',listing(SUB(ss)+2).name);
    disp(fold)
    datadir1=fullfile(datadir,fold);
    cd(datadir1);
    for rr=1:length(sessions)
        namefile=sprintf('%sR%s.edf',fold,sessions{rr});
        % read edf file from data folder
        [eeg,header]=lab_read_edf(namefile);
        fs=header.samplingrate;
        if (fs~=160)
            error('wrong sampling rate!');
        end
        % desing of butterworth filter (2-35 Hz) used as band-pass filter
        [b,a]=butter(2,[7,30]/(fs/2));
        % design of butterworth filter (59-61 Hz) used as notch filter
        [bn, an] = butter(2, [59 61]./(fs/2), 'stop');
        % detrend of EEG signals
        eeg=eeg-repmat(mean(eeg),64,1);
        % applying the filters
        for j=1:64
            eeg(j,:)=filter(b,a,eeg(j,:)-mean(eeg(j,:)));
            eeg(j,:)=filter(bn,an,eeg(j,:)-mean(eeg(j,:)));
        end
        %%% Segmentation of EEG signals
        for seg=1:size(header.events.TYP,2)
            if header.events.TYP{seg}=='T0' % Condition REST
                pos=header.events.POS(seg);
                REST(:,:,seg,rr)=eeg(:,pos+1:pos+duration);
                RR(seg,rr)=seg;
            elseif header.events.TYP{seg}=='T2' % Condition RIGHT
                pos=header.events.POS(seg);
                RIGHT(:,:,seg,rr)=eeg(:,pos+1:pos+duration);
                RI(seg,rr)=seg;
            end
        end
    end
    [ind_rest]=find(RR~=0);
    REST=reshape(REST,size(REST,1),size(REST,2),[]);
    REST=REST(:,:,[ind_rest]);
    
    [ind_right]=find(RI~=0);
    RIGHT=reshape(RIGHT,size(RIGHT,1),size(RIGHT,2),[]);
    RIGHT=RIGHT(:,:,[ind_right]);
    
    TOT=[x1,x2,y];
    rest=REST([TOT],:,:);
    right=RIGHT([TOT],:,:);
    DATA{1}=rest(:,:,:); % selection of the first 22 trials 
    DATA{2}=right;
    %%% REST and RIGHT
    
    for cc=1:length(DATA) % experimental condition
        data=DATA{cc};
        for uu=1:size(data,3) % number of trials
            % Testing for restricted form of weak-sense stationarity
            for ch=1:size(data,1)
                [~, mean_stat_flag(ch), var_stat_flag(ch)] = isstationary(squeeze(data(ch,:,uu)'));
            end
            %%% Estimating optimal order of VAR model
            maxIP=30;
            [~,IP,~,~] = mos_idVAR(data(:,:,uu),maxIP,0);
            IP=12;% average value along the 91 subjects
            %%% VAR model identification procedure
            p=nextpow2(size(data,2));
            nfft=2^p;
            [eAm,eSu,~,Residuals]=idVAR(data(:,:,uu),IP,0);
            %%% Durbin-Watson test for whiteness (H0: no serial correlation) of VAR residuals
            [dw,pval] = whiteness(data(:,:,uu),Residuals);
            
            outvar_est = fdVAR(eAm,eSu,nfft,Fs);
            eP1=outvar_est.S; %spectral matrix
            f=outvar_est.f;
            iY=[7,8,9]; %% set signals for the target
            iX1=[1,2,3]; %% set signals for the sources
            iX2=[4,5,6];
            out_t=iispectral(eP1,iY,iX1,iX2);
            
            %%% Surrogate analysis
            DATAsurr=data(:,:,uu)';
            for ns=1:iter_surr
                [X1surr,X2surr]=surriaafft2(DATAsurr,1,iX1,iX2);
                surr=[X1surr,X2surr,DATAsurr(:,[iY])];
                Data_surr(:,:,ns)=surr;
            end
            %%% distribution of Interaction Information measure along 100
            %%% surrogates
            parfor NN=1:iter_surr
                [isurr(:,NN),Isurr(NN)]=surrogate_analysis(Data_surr(:,:,NN),iY,iX1,iX2,nfft,fs);
            end
            
            clear DATAsurr Data_surr dataS
            
            
            iy_x1_x2SS(:,:,uu)=isurr;
            Iy_x1_x2SS(:,uu)=Isurr;
            iy_x1_x2f(:,uu)=out_t.iy_x1_x2;
            fy_x1f(:,uu)=out_t.fy_x1;
            fy_x2f(:,uu)=out_t.fy_x2;
            fy_x1x2f(:,uu)=out_t.fy_x1x2;
            Iy_x1_x2f(uu)=out_t.Iy_x1_x2;
            Fy_x1f(uu)=out_t.Fy_x1;
            Fy_x2f(uu)=out_t.Fy_x2;
            Fy_x1x2f(uu)=out_t.Fy_x1x2;
            
        end
        % average value along trials
        iyx1x2S(:,:,cc)=mean(iy_x1_x2SS,3);
        Iyx1x2S(:,cc)=mean(Iy_x1_x2SS,2);
        iyx1x2(:,cc)=mean(iy_x1_x2f,2);
        fyx1(:,cc)=mean(fy_x1f,2);
        fyx2(:,cc)=mean(fy_x2f,2);
        fyx1x2(:,cc)=mean(fy_x1x2f,2);
        Iyx1x2(cc)=mean(Iy_x1_x2f);
        Fyx1(cc)=mean(Fy_x1f);
        Fyx2(cc)=mean(Fy_x2f);
        Fyx1x2(cc)=mean(Fy_x1x2f);
        
        clear iy_x1_x2f fy_x1f fy_x2f fy_x1x2f Iy_x1_x2f Fy_x1f Fy_x2f Fy_x1x2f iy_x1_x2SS Iy_x1_x2SS
        
    end
    
    %%% Plot of the obtained results (Figure 5)
    freq=f(2:end/2+1);
    colR=[0.105882352941176 0.619607843137255 0.466666666666667]; % color REST
    colT=[0.850980392156863 0.372549019607843 0.00784313725490196]; % color TASK
    COL=[colR;colT];
    HCF=figure('units','inches','position',[0 0 11.7 8.3]);
    orient(HCF,'landscape')
    [ha, pos] = tight_subplot(2,3, 0.05,0.05,0.05);
    hold on
    %%% IID frequency domain analysis
    set(0,'defaultAxesFontSize', 8);
    set(0,'defaultTextFontSize', 8);
    
    F=cat(3,fyx1,fyx2,fyx1x2,iyx1x2);
    Fsurr=iyx1x2S;
    
    em_F=F; % freq x cond x measures 
    Measures={'f_{x_1;y}','f_{x_2;y}','f_{x_1x_2;y}','i_{x_1;x_2;y}'};
    for tt=1:4
        axes(ha(tt));
        for cc=1:2
            plot(freq,em_F(2:end/2+1,cc,tt),'Color',COL(cc,:),'LineWidth',1.5);
            
            hold on
        end
        title(Measures{tt})
    end

    %%%% surrogate analsysis for rest and for right conditions
    col_m=[0 0 0]; % color surrogate mean
    col_h=[17/255 17/255 17/255]; % color percentiles
    axes(ha(5))
    lo=2.5; hi=97.5;
    %%% mean surrogates along the subject
    iy_x1_x2_surr_m_tot=squeeze(mean(Fsurr,2));
    iy_x1_x2_lo_surr=squeeze(prctile(Fsurr,lo,2));
    iy_x1_x2_hi_surr=squeeze(prctile(Fsurr,hi,2));
    iy_x1_x2_m=squeeze(em_F(:,:,4));% take interaction information
    IND=[2:length(freq)+1];

    h7=area(freq,[iy_x1_x2_lo_surr([IND],1) iy_x1_x2_surr_m_tot([IND],1)-iy_x1_x2_lo_surr([IND],1) iy_x1_x2_hi_surr([IND],1)-iy_x1_x2_surr_m_tot([IND],1)]); hold on
    set(h7(1),'FaceColor','w','FaceAlpha',0); set(h7(2),'FaceColor',col_h,'FaceAlpha',0); set(h7(3),'FaceColor',col_h,'FaceAlpha',0);
    set(h7(1),'EdgeColor','k'); set(h7(2),'EdgeColor',col_h); set(h7(3),'EdgeColor',col_h);
    plot(freq,iy_x1_x2_surr_m_tot([IND],1),'Color',col_m,'LineWidth',1.5);
    set(h7(2),'FaceAlpha',0.5);
    set(h7(3),'FaceAlpha',0.5);
    hold on
    plot(freq,iy_x1_x2_m(2:end/2+1,1),'Color',COL(1,:),'LineWidth',1.5);
    
    clear h7
    axes(ha(6))
    h7=area(freq,[iy_x1_x2_lo_surr([IND],2) iy_x1_x2_surr_m_tot([IND],2)-iy_x1_x2_lo_surr([IND],2) iy_x1_x2_hi_surr([IND],2)-iy_x1_x2_surr_m_tot([IND],2)]); hold on
    set(h7(1),'FaceColor','w','FaceAlpha',0); set(h7(2),'FaceColor',col_h,'FaceAlpha',0); set(h7(3),'FaceColor',col_h,'FaceAlpha',0);
    set(h7(1),'EdgeColor','k'); set(h7(2),'EdgeColor',col_h); set(h7(3),'EdgeColor',col_h);
    plot(freq,iy_x1_x2_surr_m_tot([IND],2),'Color',col_m,'LineWidth',1.5);
    set(h7(2),'FaceAlpha',0.5);
    set(h7(3),'FaceAlpha',0.5);
    hold on
    plot(freq,iy_x1_x2_m(2:end/2+1,2),'Color',COL(2,:),'LineWidth',1.5);
    
end
