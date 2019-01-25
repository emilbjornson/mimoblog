% Illustration of the effect of beamforming when using a distributed array
%
% Emil Björnson, 2019
% emil.bjornson@liu.se

close all;
clear

%Set the wavelength
wavelength = 10;

%Size of the simulation area
squareSide = 5*wavelength;

%Generate the grid of points where the normalized channel gain is computed
x = 2:1:squareSide-2;
y = 2:1:squareSide-2;


%Put out the antennas
pointsperdim = 20; %Number of antennas per side in the square
locations = squareSide*linspace(0.025,0.975,pointsperdim); %Equally space locations of the antennas
antenna_locations= [locations 1i*locations locations+1i*squareSide squareSide+1i*locations];


%Select the target location of the beamforming
targetlocation = round(squareSide/4)+1i*round(squareSide/4);


%Prepare to save numbers
phaseshifts = zeros(length(x),length(y),length(antenna_locations));
distances = zeros(length(x),length(y),length(antenna_locations));
angletarget = zeros(length(x),length(y),length(antenna_locations));


%% Go through all antennas
for n = 1:length(antenna_locations)
    
    %Go through all spatial sample points
    for k = 1:length(x)
        
        for j = 1:length(y)
            
            %Compute the distance from the antenna to the sample point
            distances(k,j,n) = abs(x(k)+1i*y(j) - antenna_locations(n));
            
            %Compute the phase-shift from the antenna to the sample point
            phaseshifts(k,j,n) = 2*pi*distances(k,j,n)/wavelength;
            
        end
        
    end
    
    %Compute the phase-shift from the antenna to the target point
    angletarget(:,:,n) = 2*pi*abs(targetlocation - antenna_locations(n))/wavelength;
    
end

%Compute the distances from all antennas to the target point
distances_target = abs(antenna_locations-targetlocation);


%Compute the combined channel gain from all antennas to the different
%sample points. We are considering a free-space propagation pathloss
channel_gain = abs(sum(exp(1i*(phaseshifts-angletarget))./sqrt(4*pi*distances.^2),3)).^2;

%Compute the combined channel gain from all antennas to the target point
channel_gain_target = abs(sum(exp(1i*0)./sqrt(4*pi*distances_target.^2))).^2;

%Set a normalization factor for the simulations
normalization = channel_gain_target;


%% Plot simulation figure


%First figure
figure;

surf(y/wavelength,x/wavelength,channel_gain/normalization);
colormap(flipud(parula));
colorbar;
shading interp;
hold on;

for n = 1:length(antenna_locations)
    plot3(imag(antenna_locations(n))/wavelength,real(antenna_locations(n))/wavelength,0,'k*','MarkerSize',5);
end

plot3(imag(targetlocation)/wavelength,real(targetlocation)/wavelength,(channel_gain_target)/normalization,'ko','MarkerSize',5,'MarkerFaceColor','k');
hold off;

xlabel('Distance (wavelengths)');
ylabel('Distance (wavelengths)');
zlabel('Normalized channel gain');


%Second figure
figure;

surf(y/wavelength,x/wavelength,channel_gain/normalization);
colormap(flipud(parula));
colorbar;
shading interp;
hold on;

for n = 1:length(antenna_locations)
    plot3(imag(antenna_locations(n))/wavelength,real(antenna_locations(n))/wavelength,0,'k*','MarkerSize',5);
end

plot3(imag(targetlocation)/wavelength,real(targetlocation)/wavelength,(channel_gain_target)/normalization,'ko','MarkerSize',5,'MarkerFaceColor','k');
hold off;

xlabel('Distance (wavelengths)');
ylabel('Distance (wavelengths)');
zlabel('Normalized channel gain');
view(2);
