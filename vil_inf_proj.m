clear;
clc;
close all;

%% Section 1
% 5 Hz filtered data is imported with 
% Villarrica_infrasound_data.mat

load('Villarrica_infrasound_data.mat')

sps = 100; %Hz

% Assignes filtered data to channel variables
ch1f = data_filt(:,1);
ch2f = data_filt(:,2);
ch3f = data_filt(:,3);

% Channel 1 vs. 2
% Cross correlates filtered data between ch1 and ch2
oneVtwoF = xcorr(ch1f, ch2f, 'coeff');
N = length(oneVtwoF);
midpoint = ((N-1)/2)+1;
[Yaf, Iaf] = max(oneVtwoF);
lagaF = (midpoint - Iaf);

% Channel 2 vs. 3
twoVthreeF = xcorr(ch2f, ch3f, 'coeff');
[Ybf, Ibf] = max(twoVthreeF);
lagbF = (midpoint - Ibf);

% Channel 3 vs. 1
threeVoneF = xcorr(ch3f, ch1f, 'coeff');
[Ycf, Icf] = max(threeVoneF);
lagcF = (midpoint - Icf);

maxValsFiltered = [Yaf Ybf Ycf]
% Outputs sound wave lag time between channels: [1&2  2&3  3&1]
% From the lag times, we can determine approximate direction.
lagTimesFiltered = [lagaF lagbF lagcF]
fprintf('Waves are propogating to the northwest \n \n');

%% Section 2

% Position vectors relative to channel 2
% which is located the origin:
pos12 = [17 18];
pos32 = [-7 15];

sound = -340; % approximate speed of sound in m/s
              % with negative polarity pointing to where 
              % waves are going, not where they're from.

i = 0;
for theta = 1:360
    i = i + 1;
    S(i,[1 2]) = [sind(theta) cosd(theta)];
    d12(i) = dot(S(i,:),pos12);
    d32(i) = dot(S(i,:),pos32);
    t12 = (d12/sound)*sps;
    t32 = (d32/sound)*sps;
end

thetaSpace = linspace(0,359,360);
eps = (t12-lagaF).^2 + (t32-lagbF).^2;

[errY, errX] = min(eps);

% adjusts angle starting point to due N
fprintf('at an angle of %d degrees from north \n \n', errX-90);


%%% Plots %%%
figure(3)
plot(thetaSpace,eps);
title('Error Function');
% the plotted 'Angle' is computed with zero degrees at due east
% (default MatLab polarplot) so the non azimuthal angle is
% plotted here and adjusted in line 69 to show this with
% respect to due north.
xlabel('Angle');
ylabel('Epsilon');
xlim([0 360]);


figure(4)
errXrad = errX *(pi/180);
phi = [0 errXrad];
rho = [0 25];
polarplot(phi,rho, 'b>-');
hold on;
polarplot(0,0,'gd');
polarplot(-1.1342+pi,16.5529,'gd');
polarplot(0.8140,24.7588,'gd');
title('Compass Direction');
text(0,0,'ch2');
text(-1.1342+pi,16.5529,'ch3');
text(0.8140,24.7588,'ch1');
pax = gca;
angles = 0:45:360;
pax.ThetaTick = angles;
labels = {'E','NE','N','NW','W','SW','S','SE'};
pax.ThetaTickLabel = labels;
hold off;
