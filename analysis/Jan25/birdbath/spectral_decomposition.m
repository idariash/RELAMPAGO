% Ivan Arias
% 2021/01/19

% Spectral decomposition of birdbath scans

clear all
addpath(genpath('/net/denali/storage2/radar2/people/idariash/home/radartoolbox'))
addpath('/net/nasstore/students/GRAD/ECE/idariash/home/CSU/CHIVO/time_series/scr')

filename = '/net/denali/storage2/radar2/people/idariash/home/CSU/tookits/RVP9_timeSeries/res/20190125/col-radar.20190125.213000.786.BIRDBATH.1.H+V.50.mat';

Azimuth = 135;
delta_azimuth = 10;
azimuth_offset = 0;

Elevation = 85;
radarConstant = 66.660; %from netCDF file variable r_calib_radar_constant_h
PRT1 = 0.0011;
PRT2 = 0.0017;
PRT = 0.0008333;
lambda = 0.05;
chivito = chivo_timeSeries(filename);
chivito.azimuth = chivito.azimuth - azimuth_offset;
ZdrSpectrum = chivito.computeSpectrumZdr_birdbath(Azimuth, Elevation, PRT, delta_azimuth);
PowerSpectrum = chivito.computeSpectrum_Z_birdbath(Azimuth, Elevation, PRT, radarConstant, delta_azimuth);
%ZdrSpectrum = 1/4*conv2(ZdrSpectrum, [1 1; 1 1], 'same');
%PowerSpectrum = 1/4*conv2(PowerSpectrum, [1 1; 1 1], 'same');
[r, c] = size(PowerSpectrum);

Ts = PRT;
M = r;
range = 50*(1:c)/1e3; % range in km
vel=-velocityaxis(lambda,Ts,M)-0.6;

%%

ZdrSpectrum_toPlot = 10*log10(abs(ZdrSpectrum)); %10*log10(Ts/M*abs(spectrum));
ZdrSpectrum_toPlot = ZdrSpectrum_toPlot - 0.65;
PowerSpectrum_toPlot =  PowerSpectrum;

ZdrSpectrum_toPlot(PowerSpectrum_toPlot < -10) = nan;
%spectrum_toPlot = 10*log10(Ts/M*abs(spectrum));
%spectrum_toPlot(abs(spectrum_toPlot) > 10) = nan;
height = range*sind(Elevation);
figure
pcolor(vel, range, ZdrSpectrum_toPlot')
shading flat
colormap jet
colorbar
xlabel v(m/s)
ylabel Range(km)
ylim([0 20])
caxis([-2 6])
grid on
title(['Azi: ' num2str(Azimuth) ' | Elv: ' num2str(Elevation)])
hc=colorbar;
ylabel(hc,'Diff. Reflectivity (dB)');

PowerSpectrum_toPlot(PowerSpectrum_toPlot < -10) = nan;
figure
pcolor(vel, range, PowerSpectrum_toPlot')
shading flat
colormap jet
colorbar
xlabel v(m/s)
ylabel('Height (km)')
ylim([0 10])
grid on
%title(['Azi: ' num2str(Azimuth) ' | Elv: ' num2str(Elevation)])
hc=colorbar;
ylabel(hc,'Reflectivity (dBZ)');


%%
h_0 = 2.0;
delta = 0.05;
smooth_filter = ones(1, 10);
smooth_filter = smooth_filter/sum(smooth_filter);
spectrum_gate = PowerSpectrum(:, abs(range - h_0) < delta);
spectrum_gate = nanmean(spectrum_gate');
spectrum_gate = conv(spectrum_gate, smooth_filter, 'same');
figure(3)
plot(vel, spectrum_gate)
grid on
xlabel('velocity (m/s)')
ylabel('Reflectivity (dBZ)')
title(['Range: ' num2str(h_0) ' km'])
xlim([-15,15])
ylim([-10,70])


spectrum_gate_zdr = ZdrSpectrum_toPlot(:, abs(range - h_0) < delta);

spectrum_gate_zdr = nanmean(spectrum_gate_zdr');
spectrum_gate_zdr = conv(spectrum_gate_zdr, smooth_filter, 'same');
%spectrum_gate_zdr(spectrum_gate < -1) = nan;
Index = (abs(vel - (-10)) < 4) | (abs(vel - 5) < 5);
spectrum_gate_zdr(~Index) = nan;
figure(4)
plot(vel, spectrum_gate_zdr)
grid on
xlabel('velocity (m/s)')
ylabel('Diff. Reflectivity (dB)')
title(['Range: ' num2str(h_0) ' km'])
xlim([-15,15])
ylim([-2,8])