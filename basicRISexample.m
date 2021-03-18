close all;
clear;


%% Set simulation parameters

%Set transmit power (in mW)
P = 10;

%Bandwidth
B = 20e6;


%Noise figure (in dB)
noiseFiguredB = 13;


%Compute the noise power in dBm
sigma2dBm = -174 + 10*log10(B) + noiseFiguredB;
sigma2 = db2pow(sigma2dBm);


%Set the amplitude reflection coefficient
alpha = 1;


%Define the range of different number of reflection elements in the RIS
Nrange = 0:1000;


%Pathloss from transmitter to RIS
beta_Tx_RIS = 10^(-7.5);

%Pathloss from RIS to receiver
beta_RIS_Rx = 10^(-7.5);


%Prepare to save transmit powers
rate_RIS = zeros(length(Nrange),2);
rate_SISO = zeros(length(Nrange),2);


%Go through all number of RIS elements
for k = 1:length(Nrange)
    
    %Extract number of RIS elements
    N = Nrange(k);
    
    for c = 1:2
        
        %Two cases for the direct path
        if c == 1
            beta_direct = 10^(-10);
        elseif c == 2
            beta_direct = 10^(-7.5);
        end
        
        %Rate without RIS
        rate_SISO(k,c) = log2(1+ P*beta_direct/sigma2);
        
        %Rate with RIS
        rate_RIS(k,c) = log2(1 + P*(sqrt(beta_direct) + N*alpha*sqrt(beta_Tx_RIS*beta_RIS_Rx)).^2/sigma2);
        
    end
    
end


%Plot simulation results
figure;
hold on; box on;

plot(Nrange,rate_RIS(:,1),'r-','LineWidth',2);
plot(Nrange,rate_SISO(:,1),'k--','LineWidth',2);
plot(Nrange,rate_RIS(:,2),'r-','LineWidth',2);
plot(Nrange,rate_SISO(:,2),'k--','LineWidth',2);
xlabel('Number of elements','Interpreter','Latex');
ylabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
legend({'With RIS','Baseline'},'Interpreter','latex','Location','Best');
set(gca,'fontsize',18);
