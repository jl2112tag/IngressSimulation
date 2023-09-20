%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liquid ingress in porous media simulation tool
% Jongmin Lee, 2022
% Terahertz applications group, 
% University of Cambridge, UK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% define variables
capillaryNum = 50; % number of capillaries
layerNum = 100; % number of capillary layers
capillaryHeight = 10; % capillary height for each layer in micrometers
interlayerConnectivity = 0.1; % pore connection probablity ranged from 0 to 1
dynamicViscosity = 6.922e-4; % water dynamic viscosity at 37 degress, Pa*s
surfaceTension = 70.41e-3; % water in contact with air at 35 degress, Nm^-1
contactAngle = 30; % degree
numOfBins = 10;

% capillary penetration time in a 3D structure
% case 1: normal distribution with a mean (meanRadius) and standard deviation (sigmaRadius) in micrometers
meanRadius = 5; % capillary radius in micrometer
sigmaRadius = 1; % capillary radius standard deviation
x=linspace(0,meanRadius*2,capillaryNum);
radiusPDF = normpdf(x,meanRadius,sigmaRadius); % capillary radius distribution function, x-axis: radius in micrometers, y-axis: probability

xq = linspace(0,max(cumtrapz(radiusPDF)),capillaryNum+1);
xq = (xq(1:end-1)+xq(2:end))/2;
radiusVector = interp1(cumtrapz(radiusPDF),x,xq);

pTimeCoef = 2*dynamicViscosity*(capillaryHeight*10^-6)^2/(surfaceTension*cos(deg2rad(contactAngle)));
pTimeVector = pTimeCoef./(radiusVector*10^-6); % penetration time vector in seconds

pTimeMat = repmat(pTimeVector,capillaryNum,1);
pTimeMatX = transpose(pTimeMat);
pTimeMatUnit(:,:,1) = pTimeMat;
pTimeMatUnit(:,:,2) = pTimeMatX;
pTimeVol = repmat(pTimeMatUnit,1,1,ceil(layerNum/2));


% 3D capillary connectivity map
mat_C = ones(capillaryNum,capillaryNum,layerNum);
mat_P = rand(size(mat_C));
mat_C(mat_P < interlayerConnectivity) = 0; % 1 for connection nodes
