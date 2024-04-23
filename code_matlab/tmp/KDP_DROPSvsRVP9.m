% Ivan Arias
% KPD comparison
% 2019/02/21

close all
filename_rvp9 = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/NetCDF/20181126/cfrad.20181126_071002.505_to_20181126_071536.425_col-radar_PPICOWH360_SUR.nc';
filename_drops = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/DROPS/20181126_07to08UTC/cfrad.20181126_071002.505_to_20181126_071536.425_col-radar_PPICOWH360_SUR.nc'; 

KDP_rvp9 = ncread(filename_rvp9, 'KDP');
KDP_drops = ncread(filename_drops, 'KDP');
range = ncread(filename_drops, 'range');

RHOHV_rvp9 = ncread(filename_rvp9, 'RHOHV');
RHOHV_drops = ncread(filename_drops, 'RHOHV');

SNR_rvp9 = ncread(filename_rvp9, 'SNR');

DBZ_drops = ncread(filename_drops, 'DBZ');
DBZ_rvp9 = ncread(filename_rvp9, 'DBZ');


KDP_rvp9(KDP_rvp9 < -10) = nan;
KDP_drops(KDP_drops < -10) = nan;

KDP_rvp9(RHOHV_rvp9 < 0.95) = nan;
KDP_drops(RHOHV_drops < 0.95) = nan;

KDP_rvp9(SNR_rvp9 < 10) = nan;
KDP_drops(SNR_rvp9 < 10) = nan;


%% First sweep or second, depends

KDP_rvp9_1st_sweep = KDP_rvp9(360*length(range)+1:720*length(range));
% KDP_rvp9_1st_sweep = vec2mat(KDP_rvp9_1st_sweep, length(range));
% KDP_rvp9_1st_sweep = KDP_rvp9_1st_sweep';

KDP_drops_1st_sweep = KDP_drops(360*length(range)+1:720*length(range));
% KDP_drops_1st_sweep = vec2mat(KDP_drops_1st_sweep, length(range));
% KDP_drops_1st_sweep = KDP_drops_1st_sweep';

DBZ_rvp9_1st_sweep = DBZ_rvp9(360*length(range)+1:720*length(range));
% DBZ_rvp9_1st_sweep = vec2mat(DBZ_rvp9_1st_sweep, length(range));
% DBZ_rvp9_1st_sweep = DBZ_rvp9_1st_sweep';

DBZ_drops_1st_sweep = DBZ_drops(360*length(range)+1:720*length(range));
% DBZ_drops_1st_sweep = vec2mat(DBZ_drops_1st_sweep, length(range));
% DBZ_drops_1st_sweep = KDP_drops_1st_sweep';


%% Indexing KDP behaivor

%KDP_rvp9(KDP_rvp9 < 0.009 | KDP_rvp9 > 0.021) = nan;

KDP_indexing = (KDP_rvp9_1st_sweep > 0.009 & KDP_rvp9_1st_sweep < 0.021);

DBZ_drops_KDP = DBZ_drops(KDP_indexing);
DBZ_drops_KDP(DBZ_drops_KDP < -1000) = nan;

DBZ_rvp9_KDP = DBZ_rvp9(KDP_indexing);
DBZ_rvp9_KDP(DBZ_rvp9_KDP < -1000) = nan;

%KDP_drops_1st_sweep(KDP_indexing) = nan;

%Ploting


scatter(KDP_drops_1st_sweep, KDP_rvp9_1st_sweep);
x = -2:4;
hold on 
plot(x,x)
figure
hist(DBZ_drops_KDP,16);
figure
hist(DBZ_rvp9_KDP,16);




