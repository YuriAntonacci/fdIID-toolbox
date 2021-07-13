%genera x y surrogati con la phase randomisation mantenedo phase lags - mantiene power spectrum e cross-spectrum
%NB: il size della serie pi� lunga deve essere pari

function xys=surrfft2(x,y)

szx=size(x);szx=szx(1);
szy=size(y);szy=szy(1);
if szx>szy
   pad=zeros(szx-szy,1);
   y=[y;pad];
   szy=szx;
end
if szy>szx
   pad=zeros(szy-szx,1);
   x=[x;pad];
   szx=szy;
end

fx=fft(x);
fy=fft(y);
mx=abs(fx);
my=abs(fy);
fasex=angle(fx);
fasey=angle(fy);

fasesurra=2*pi*rand(szx/2-1,1)-pi;
fasesurrb=-flipdim(fasesurra,1);

fasex(2:szx/2,:)=fasex(2:szx/2,:)+repmat(fasesurra,1,size(fx,2));
fasex(szx/2+2:szx,:)=fasex(szx/2+2:szx,:)+repmat(fasesurrb,1,size(fx,2));
fasey(2:szx/2,:)=fasey(2:szx/2,:)+repmat(fasesurra,1,size(fy,2));
fasey(szx/2+2:szx,:)=fasey(szx/2+2:szx,:)+repmat(fasesurrb,1,size(fy,2));

bin=2*pi*rand(1)-pi;
fasex(1,:)=fasex(1,:)+bin;fasex(szx/2+1,:)=fasex(1,:);
fasey(1,:)=fasey(1,:)+bin;fasey(szx/2+1,:)=fasey(1,:);

fxs=mx.*(cos(fasex)+j*sin(fasex));
fys=my.*(cos(fasey)+j*sin(fasey));
xs=ifft(fxs);xs=real(xs);xs=xs-mean(xs,1);
ys=ifft(fys);ys=real(ys);ys=ys-mean(ys,1);
xys=[xs ys];