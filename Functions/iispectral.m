%% frequency domain source interaction measures
% input: spectral matrix and dimensions of blocks
% output frequency domain measures of linear association and source interaction

function out=iispectral(P,iY,iX1,iX2)

M=size(P,1);
My=length(iY);
M1=length(iX1);
M2=length(iX2);
assert(M1+M2+My==M)

%% partition spectral matrix in three blocks
Py=P(iY,iY,:);
Px1=P(iX1,iX1,:);
Px2=P(iX2,iX2,:);
Pyx1=P([iY iX1],[iY iX1],:);
Pyx2=P([iY iX2],[iY iX2],:);
Px1x2=P([iX1 iX2],[iX1 iX2],:);

%% block coherences
Nf=size(P,3);
dyn=nan*ones(Nf,1); dx1n=dyn; dx2n=dyn; dyx1n=dyn; dyx2n=dyn; dx1x2n=dyn; dyx1x2n=dyn;
for n=1:Nf
    dy=abs(det(Py(:,:,n)));  dyn(n)=dy;
    dx1=abs(det(Px1(:,:,n))); dx1n(n)=dx1;
    dx2=abs(det(Px2(:,:,n))); dx2n(n)=dx2;
    dyx1=abs(det(Pyx1(:,:,n))); dyx1n(n)=dyx1;
    dyx2=abs(det(Pyx2(:,:,n))); dyx2n(n)=dyx2;
    dx1x2=abs(det(Px1x2(:,:,n))); dx1x2n(n)=dx1x2;
    dyx1x2=abs(det(P(:,:,n))); dyx1x2n(n)=dyx1x2;
    
 
%     Cy_x1(n)=abs(1-dyx1/(dy*dx1));
%     Cy_x2(n)=abs(1-dyx2/(dy*dx2));
%     Cy_x1x2(n)=abs(1-dyx1x2/(dy*dx1x2));
end

% epsi=1/10000;
% dyx1n(dyx1n<epsi)=0;
% dyx2n(dyx2n<epsi)=0;
% dyx1x2n(dyx1x2n<epsi)=0;

%% non-logarithmic measures of association and source interaction
Cy_x1=(1-dyx1n./(dyn.*dx1n));
Cy_x2=(1-dyx2n./(dyn.*dx2n));
Cy_x1x2=(1-dyx1x2n./(dyn.*dx1x2n));

Cy_x1_x2=Cy_x1+Cy_x2-Cy_x1x2;

%% logarithmic measures of association and source interaction
fy_x1=-log(1-Cy_x1);
fy_x2=-log(1-Cy_x2);
fy_x1x2=-log(1-Cy_x1x2);

iy_x1_x2=fy_x1+fy_x2-fy_x1x2;

%% time domain measures
% ( total EEG power: var(x)=(2/fs)*(integral of Px from 0 to fs/2) -> (2/fs) * (sum(Px) * (fs/2)/Nf ) = sum(Px)/Nf )
Fy_x1=sum(fy_x1)/Nf;
Fy_x2=sum(fy_x2)/Nf;
Fy_x1x2=sum(fy_x1x2)/Nf;
Iy_x1_x2=Fy_x1+Fy_x2-Fy_x1x2;


%% out
out.iy_x1_x2=iy_x1_x2;
out.fy_x1x2=fy_x1x2;
out.fy_x1=fy_x1;
out.fy_x2=fy_x2;
out.Iy_x1_x2=Iy_x1_x2;
out.Fy_x1x2=Fy_x1x2;
out.Fy_x1=Fy_x1;
out.Fy_x2=Fy_x2;
out.Cy_x1=Cy_x1;
out.Cy_x2=Cy_x2;
out.Cy_x1x2=Cy_x1x2;
out.Cy_x1_x2=Cy_x1_x2;

end