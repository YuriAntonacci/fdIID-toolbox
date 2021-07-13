% genera iterative amplitude adjusted fourier tranform surrogates 
% algoritmo di Schreiber e Schmitz - Physical Review Letters 1996
% lo faccio per due serie e provo a mantenere il cross-spettro

% DATA: serie da surrogare
% iX1=[vettore sorgente1 - anche multivariata]
% iX2=[vettore sorgente2 - anche multivariata]
% nit: numero di iterazioni volute (default 7)
% stop: se metto 'spe' esce con lo spettro conservato, se metto 'dis' esce con la distribuzione conservata

function [xs,ys]=surriaafft2(DATA,nit,iX1,iX2,stop)

error(nargchk(1,4,nargin));%min e max di input arguments
if nargin < 5, stop='spe'; end %default matcha lo spettro
if nargin < 3, nit=7; end %default 7 iterazioni

% clear;close all;
% percorso='D:\johnny\lavoro\integrate_nlpred\elaborati_loo_si\';% percorso dei dati da analizzare
% nomefile='b-ca.prn';
% rs=load([percorso nomefile]);
% y=rs(:,1);
% y=(y-mean(y))/std(y);
x=DATA(:,[iX1]);
y=DATA(:,[iX2]);

%% inizializzazione
for xx=1:size(x,2)
    [xsorted(:,xx),xindice]=sort(x(:,xx));
    mx(:,xx)=abs(fft(x(:,xx)));
end

for yy=1:size(y,2)
    [ysorted(:,yy),yindice]=sort(y(:,yy));
    my(:,yy)=abs(fft(y(:,yy)));
end

xys=surrfft2(x,y);
xs=xys(:,1:size(x,2));
ys=xys(:,size(x,2)+1:end);


%% ciclo

for i=1:nit
    
    % step 1: impone lo spettro
    fasexs=angle(fft(xs));
    fxs=mx.*(cos(fasexs)+j*sin(fasexs));
    xs=ifft(fxs);
    xs=real(xs);
    xs=xs-mean(xs,1);
    
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+j*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys,1);

    % step 2: impone la distribuzione
    [xssorted,xsindice]=sort(xs);
    xpermuted=zeros(length(x),size(x,2));
    for i=1:length(x)
        for ll=1:size(xpermuted,2)
            xpermuted(xsindice(i,ll),ll)=xsorted(i,ll);
        end
    end
    xs=xpermuted;

    
    [yssorted,ysindice]=sort(ys);
    ypermuted=zeros(length(y),size(y,2));
    for i=1:length(y)
        for ll=1:size(ypermuted,2)
            ypermuted(ysindice(i,ll),ll)=ysorted(i,ll);
        end
    end
    ys=ypermuted;

end

%se volevo conservare lo spettro, faccio 1 altro mezzo giro dove impongo solo quello
if stop=='spe'
    fasexs=angle(fft(xs));
    fxs=mx.*(cos(fasexs)+j*sin(fasexs));
    xs=ifft(fxs);xs=real(xs);
    xs=xs-mean(xs,1);
    
    faseys=angle(fft(ys));
    fys=my.*(cos(faseys)+j*sin(faseys));
    ys=ifft(fys);ys=real(ys);
    ys=ys-mean(ys,1);
end




