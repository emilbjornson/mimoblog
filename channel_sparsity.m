% Illustration of a 64-antenna array's ability to resolve the channel paths
% from three distinct angular directions. The script uses the function
% functionSpatialSignature3DLoS.m from the book "Massive MIMO networks".
%
% It was used in the blog post "Channel Sparsity in Massive MIMO"
% https://ma-mimo.ellintech.se/2019/09/22/channel-sparsity-in-massive-mimo/
%
% Emil Björnson, 2019
% emil.bjornson@liu.se

close all;
clear;


%Define the number of antennas in the array in three cases
M_H_range = [1 8 64]; %Number of antennas per row
M_V_range = [64 8 1]; %Number of antennas per column

%Define the range of azimuth angles to be considered when computing the
%received signal gain
varphiRangeDeg = (-90:0.1:90);
varphiRange = varphiRangeDeg*(pi/180);

%Define the channel properties
varphiDeg = [30 -20 40]; %Azimuth angles of the true channel paths
varphi = varphiDeg*(pi/180);
theta = 0; %Common elevation angle of the true channel paths


%Wavelength (normalized)
lambda = 1;

%Go through the three cases
for k = 1:length(M_H_range)
    
    
    M_H = M_H_range(k); %Exctract number of antennas per row
    M_V = M_V_range(k); %Exctract number of antennas per column
    d_H = 0.5*lambda; %Horizontal antenna spacing
    d_V = 0.5*lambda; %Vertical antenna spacing
    
    %Define the antenna geometry
    M = M_H*M_V; %Total number of antennas
    U = zeros(3,M); %Matrix containing the position of the antennas
    
    i = @(m) mod(m-1,M_H); %Horizontal index
    j = @(m) floor((m-1)/M_H); %Vertical index
    
    for m = 1:M
        U(:,m) = [0; i(m)*d_H; j(m)*d_V]; %Position of the mth element
    end
    
    
    %Compute the array response for various directions (using a function
    %developed in the book Massive MIMO networks)
    hRange = zeros(M,length(varphiRange));
    for n = 1:length(varphiRange)
        hRange(:,n) = functionSpatialSignature3DLoS(U,varphiRange(n),theta,lambda);
    end
    
    %Compute the channel by adding up the array responses for the three
    %different paths. The last two paths have only half the amplitude, to
    %simply model that not all paths must be equally strong
    h = functionSpatialSignature3DLoS(U,varphi(1),theta,lambda);
    for n = 2:length(varphi)
        h = h + 0.5*functionSpatialSignature3DLoS(U,varphi(n),theta,lambda);
    end
    
    %Compute how the received signal power is distributed over different
    %azimuth directions. This is the channel gain that will be plotted
    gains = abs(h'*hRange).^2;
    
    %Plot the results
    figure(1);
    hold on; box on;
    plot(varphiRangeDeg,10*log10(gains),'LineWidth',1);
    xlabel('Angle of arrival');
    ylabel('Signal gain [dB]');
    
end

legend('$1 \times 64$','$8 \times 8$','$64 \times 1$','Interpreter','latex','Location','SouthWest');

for n = 1:length(varphiDeg)
    plot(varphiDeg(n)*[1 1],[-20 50],'k');
end

xlim([-90 90]);
ylim([-20 50]);
set(gca,'fontsize',14);
