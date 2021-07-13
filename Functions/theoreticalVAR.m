%% THEORETICAL COEFFICIENTS FOR SIMULATED VAR PROCESSES

function [Am,Su,Ak,z]=theoreticalVAR(M,par)
% clear; close all;clc;
% M=2;
% par.poles{1}=([0.8 0.1]); % Autonomous oscillations in process 1
% % par.poles{2}=([0.8 0.1; 0.9 0.3]); % Autonomous oscillations in process 2
% par.poles{2}=([0.9 0.3]); % Autonomous oscillations in process 2
% par.coup=[1 2 1 0.5; 2 1 2 -0.2]; % in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
% par.Su=[1 1]; % variance of innovation processes

Su=eye(M); %innovation covariance matrix (diagonal)
for k=1:M
    Su(k,k)=par.Su(k);
end

for m=1:M
    npoles=size(par.poles{m},1);
    z{m}=[]; % list of poles
    for n=1:npoles
        r=par.poles{m}(n,1);
        f=par.poles{m}(n,2);
        if f==0 % add real pole
            z{m}=[z{m}; r];
        else %add complex conjugate pole
            za=r*(cos(2*pi*f)+1i*sin(2*pi*f)); %pole
            zb=r*(cos(2*pi*f)-1i*sin(2*pi*f)); %complex conjugate pole
            z{m}=[z{m}; za; zb];
        end
    end
    Apol=poly(z{m});
    Amd{m}=-Apol(2:length(Apol));
    pd(m)=length(Amd{m}); %max delay for autonomous oscillations in process m
end

if isempty(par.coup)
    p=max(pd);
else
    p=max(max(pd),max(par.coup(:,3))); %maximum delay
end

Ak=zeros(M,M,p); %blocks of coefficients
for m=1:M
    for k=1:length(Amd{m})
        Ak(m,m,k)=Amd{m}(k);
    end
end
for k=1:size(par.coup,1)
    Ak(par.coup(k,1),par.coup(k,2),par.coup(k,3))=par.coup(k,4);
end

Am=[]; %group coefficient blocks Ak in a single matrix Am
for kk=1:p
    Am=[Am; Ak(:,:,kk)];
end

% stability check
Am1=Am';
E=eye(M*p);AA=[Am1;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
if lambdamax>=1
    error('The simulated VAR process is not stable');
end


end


