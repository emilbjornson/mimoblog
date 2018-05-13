% Demonstration of out-of-band radiation in Massive MIMO
% using real amplifier in the lab, through WebLab
%
% For theory, see https://arxiv.org/pdf/1802.02475.pdf
% For WebLab, see http://dpdcompetition.com/rfweblab/
%
% Erik G. Larsson, 2018
% erik.g.larsson@liu.se

clear all
close all

N=20000;  % number of samples

% generate first sinusoid
f1=pi/10;
z1=exp(j*[0:N-1]*f1);
%y1=weblab(z1(:)*0.01);
%close all
%plot(180/pi*unwrap(angle(z1)))

% generate second sinusoid
f2=2*pi/10;
z2=exp(j*[0:N-1]*f2);
%y2=weblab(z2(:)*0.01);
%hold on
%plot(180/pi*unwrap(angle(z2)))
 
% amplify signals for array with line-of-sight
M=50;  % number of antennas
Y=zeros(N,M);
phi1=-0.17;   % angle-of-arrival for first sinusoid
phi2=0.6;     % angle-of-arrival for second sinusoid
for m=1:M
    % compute signal transmitted by antenna m
    z=exp(j*pi*sin(phi1)*(m-1))*z1 + exp(j*pi*sin(phi2)*(m-1))*z2;
    % run through the amplifier via weblab
    Y(:,m)=weblab(z(:)*0.03);

    % scaling 0.03 seems to be a reasonable operating point
    % if pushed up too high, amplifier might burn... (but weblab
    % gives error message before this happens)
    
    % uncomment the line below for cross-validation with polynomial model
    % Y(:,m)=z+0.2*z.*(abs(z).^2);  
end

save weblabout.mat

% Perform rough spectral analysis of the result
L=500;
S=linspace(-pi/2,pi/2,L);
p1=zeros(L,1);
p2=zeros(L,1);
Ptot=zeros(L,1);
for l=1:L
  a = exp(j*[0:M-1]*pi*sin(S(l)));
  a = a(:);
  y=a'*Y.';
  x=abs(fft(y)).^2;
  p1(l)=sum(x(round(N/2*f1/pi)-5:round(N/2*f1/pi)+5));
  p2(l)=sum(x(round(N/2*f2/pi)-5:round(N/2*f2/pi)+5));
  Ptot(l)=sum(x);
end

set(0, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');

figure(1)
clf
semilogy(180/pi*S,p1,'k')
hold on
semilogy(180/pi*S,p2,'r')
semilogy(180/pi*S,Ptot,'b')
line([phi1 phi1]*180/pi,[1e15 1e5])
line([phi2 phi2]*180/pi,[1e15 1e5])
xlabel('angle relative to array boresight [degrees]')
ylabel('relative radiated power [dB]')
legend('at $f_1$','at $f_2$','total power over the entire band')

boldify
