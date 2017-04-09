clear all; close all;

L=20;
Fc=100;
dt=1/Fc;
f0=15;
t=(0:L-1)*dt;
t(floor(length(t)/2))=t(floor(length(t)/2))+dt/2; %errore campionamento
xc=@(t,f0) sin(2*pi*f0*t).*(t); %segnale analogico
xd=xc(t,f0); %segnale discreto

N=10;
FcN=Fc*N; %nuovo campionamento
LN=N*L;
tN=(0:LN-1)/FcN;
xdN=xc(tN,f0);

%resampling

ydN=HermiteSplineInterpolation(xd,t,tN);
figure(1); plot(tN,ydN,'-',tN,xdN,'--',t,xd,'o');

t2=sort(rand(1,L)*L/Fc);
xd2=xc(t2,f0);
ydN2=HermiteSplineInterpolation(xd2,t2,tN);
figure(2); plot(tN,ydN2,'-',tN,xdN,'--',t2,xd2,'o');
