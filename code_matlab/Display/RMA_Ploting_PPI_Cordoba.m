% Ivan Arias (C)
% Plot RMA1, Cordoba Radar
% May-22/2018


addpath('/net/denali/storage2/radar2/tmp/Ivan/Utiles/PPI_Plotting');
addpath('/net/denali/storage2/radar2/tmp/Ivan/Utiles/Matlab');
DataPath = '/home/idariash/Desktop/Link to Ivan/PuertoRico/Data/NetCDF/20170429/';
figure_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Radar_Cordoba/Figures/';
Directory = dir([DataPath, '*.nc']);
rhoHV_threshold = 0.9;
for N = 1:4 %length(Directory)
    close all
    filename = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Data/20171129/cfrad.20171129_192316.000_to_20171129_192432.000_1_SUR.nc'; %[DataPath, Directory(i).name];
    RadarType = 'NexRad';
    Campaign = 'Cordoba';
    close all

    % get the fields
    Zh = ncread(filename, 'TH');
    Zv = ncread(filename, 'TV');
    ZDR = ncread(filename, 'TDR');
    VEL = ncread(filename, 'VRAD');
    RHOHV = ncread(filename, 'RHOHV');
    PHIDP = ncread(filename, 'PHIDP');
    SW = ncread(filename, 'WRAD');
    KDP = ncread(filename, 'KDP');

    %Zh(RHOHV < 0.5) = nan;
    %ZDR(RHOHV < 0.9) = nan;

    %Note: need to change the range in km not in m
    range = ncread(filename, 'range')/1e3;

    % get the elevation information
    elv = ncread(filename, 'elevation');
    azi = ncread(filename, 'azimuth');


    %configuration to show Nth sweep 
    elv = elv( (N-1)*360 + 1 : N*360 );
    azi = azi( (N-1)*360 + 1 : N*360 );
    
    Zh_Nth_sweep = Zh( : ,  (N-1)*360 + 1 : N*360 );
    Zv_Nth_sweep = Zh( : ,  (N-1)*360 + 1 : N*360 );
    ZDR_Nth_sweep = ZDR( : ,  (N-1)*360 + 1 : N*360 );
    VEL_Nth_sweep = VEL( : ,  (N-1)*360 + 1 : N*360 );
    RHOHV_Nth_sweep = RHOHV( : ,  (N-1)*360 + 1 : N*360 );
    PHIDP_Nth_sweep = PHIDP( : ,  (N-1)*360 + 1 : N*360 );
    KDP_Nth_sweep = KDP( : ,  (N-1)*360 + 1 : N*360 );
    SW_Nth_sweep = SW( : ,  (N-1)*360 + 1 : N*360 );

    

    % Plot
    figure(1); %DBZh
    Zh_Nth_sweep( Zh_Nth_sweep<-20) = nan;
    fig_01 = plotPPI(range', azi', Zh_Nth_sweep, 'dBZ', [-30 70], .1, 25,'nwsz',1);
    %[date, hour] = get_time(filename, RadarType, Campaign);
    str = ['Horizontal Reflectivity, Elevation:  ' num2str(elv(4)) '°'];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(fig_01 ,fullfile([figure_path, '/DBZh/'], ['DBZh_Elv', num2str(elv(4)), '.png']),'png');
    
    figure(2); %DBZv
    Zv_Nth_sweep(Zv_Nth_sweep<-20) = nan;
    fig_02 = plotPPI(range', azi', Zv_Nth_sweep, 'dBZ', [-30 70], .1, 25,'nwsz',1);
    str = ['Vertical Reflectivity, Elevation:  ' num2str(elv(4)) '°'];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(fig_02 ,fullfile([figure_path, '/DBZv/'], ['DBZv_Elv', num2str(elv(4)), '.png']),'png');
    
    figure(3); %ZDR 
    ZDR_Nth_sweep(ZDR_Nth_sweep<-10) = nan;
    fig_03=plotPPI(range', azi', ZDR_Nth_sweep, 'dB', [-12 12], 1, 25,'nwsz',1);
    str = ['Differential Reflectivity, Elevation:  ' num2str(elv(4)) '°'];
    title(str, 'Color', 'k'); 
    %grid on; 
    %colormap('parula')
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(fig_03 ,fullfile([figure_path, '/ZDR/'], ['ZDR_Elv', num2str(elv(4)), '.png']),'png');
    
    figure(4); %VEL
    VEL_Nth_sweep(VEL_Nth_sweep<-150) = nan;
    fig_04 = plotPPI(range', azi', VEL_Nth_sweep, 'm/s', [-50 50], .01, 25,'nwsv',1);
    str = ['Radial Velocity, Elevation:  ' num2str(elv(4)) '°'];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(fig_04 ,fullfile([figure_path, '/VEL/'], ['VEL_Elv', num2str(elv(4)), '.png']),'png');
    
    figure(5);  % RHOHV
    RHOHV_Nth_sweep(RHOHV_Nth_sweep<-2) = nan;
    fig_05 = plotPPI(range', azi', RHOHV_Nth_sweep, '\rho_{hv}', [-2 2], .01, 25,'jet',1);
    str = ['RHOHV, Elevation:  ' num2str(elv(4)) '°'];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(fig_05 ,fullfile([figure_path, '/RHOHV/'], ['RHOHV_Elv', num2str(elv(4)), '.png']),'png');

    figure(6); % PHIDP
    PHIDP_Nth_sweep(PHIDP_Nth_sweep<-180) = nan;
    fig_06=plotPPI(range', azi', PHIDP_Nth_sweep, 'deg', [-10 180], .1, 25,'jet',1);
    str = ['PHIDP, Elevation:  ' num2str(elv(4)) '°'];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(fig_06 ,fullfile([figure_path, '/PHIDP/'], ['PHIDP_Elv', num2str(elv(4)), '.png']),'png');

    figure(7); %KDP
    KDP_Nth_sweep(KDP_Nth_sweep<-10) = nan;
    fig_07 = plotPPI(range', azi', KDP_Nth_sweep, 'deg/km', [-3 3], .01, 25,'winter',1);
    str = ['KDP, Elevation:  ' num2str(elv(4)) '°'];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    colormap winter
    saveas(fig_07 ,fullfile([figure_path, '/KDP/'], ['KDP_Elv', num2str(elv(4)), '.png']),'png');
 
    figure(8) % Scatter ZvsZDR
    Z_vector = Zv_Nth_sweep;
    ZDR_vector = ZDR_Nth_sweep;
    ZDR_vector(ZDR_Nth_sweep<-3) = NaN;
    %ZDR_vector(ZDR_Nth_sweep > 8) = NaN;
    ZDR_vector(RHOHV_Nth_sweep < 0.9) = NaN;
    Z_vector(RHOHV_Nth_sweep < 0.9) = NaN;
    Z_vector = reshape(Z_vector, 1, 360*484);
    ZDR_vector = reshape(ZDR_vector, 1, 360*484);
    ZvsZDR = scatter(Z_vector, ZDR_vector);
    xlabel('Reflectivity (dBZ)')
    ylabel('Differential Reflectivity (dBZ)')
    grid on
    str = [hour, ' UTC, ZvsZDR, ', ' SQI: ', num2str(SQI_threshold), ', rhoHV: ', num2str(rhoHV_threshold)];
    title(str, 'Color', 'k');
    saveas(ZvsZDR,fullfile([figure_path, '/Scatter_ZvsZDR/'], ['KDP_Elv', num2str(elv(4)), '.png']),'png');
    
    % fig_name= sprintf('Zm');
    % path_f='N:\tmp\Ivan\D3R\IGARSS_Paper\Plots\';
    % saveas(aa,fullfile(path_f,[campaign, '_', fig_name]),'png');
    % 
    % figure(2); 
    % bb=plotPPI(range, azi', double(data.cZ), 'dBZ', [-30 75], 2.5, 10,'nwsz',1);
    % %title('Intrinsic reflectivity', 'Color', 'k'); 
    % grid on
    % set(gcf, 'Position', get(0,'Screensize'));
    % xlim([-40 40]);
    % ylim([-40 40]);
    % fig_name = sprintf('Z');
    % saveas(bb,fullfile(path_f,[campaign, '_', fig_name]),'png');
    % % 
    % % %
    % figure(3); 
    % data.RainRate = ncread(filename, 'RainRate');
    % cc=plotPPI(range, azi', double(data.RainRate), 'mm/h', [0 25], 2.5, 10, 1);
    % %title('Rain Rate', 'Color', 'k'); grid on
    % set(gcf, 'Position', get(0,'Screensize'));
    % xlim([-40 40]);
    % ylim([-40 40]);
    % fig_name= sprintf('RainRate');
    % saveas(cc,fullfile(path_f,[campaign, '_', fig_name]),'png');

    % %
    % figure; 
    % dd=plotPPI(range, azi', double(data.SW), '', [-1 4], 0.25, 5,1);
    % title('Spectral Width', 'Color', 'k'); grid on
    % set(gcf, 'Position', get(0,'Screensize'));
    % xlim([-40 40]);
    % ylim([-40 40]);
    % fig_name= sprintf('SpectralWidth');
    % saveas(dd,fullfile(path_f,fig_name),'png');
    % 
    % %
    % figure; 
    % ee=plotPPI(range, azi', double(data.RH), ' ', [0.92 1.00], 0.005, 5,'idl',1);
    % title('CoPol Correlation', 'Color', 'k'); grid on
    % set(gcf, 'Position', get(0,'Screensize'));
    % xlim([-40 40]);
    % ylim([-40 40]);
    % fig_name= sprintf('CoPolCorrelation');
    % saveas(ee,fullfile(path_f,fig_name),'png');
    % 
    % %
    % figure; 
    % ff=plotPPI(range, azi', double(data.VEL), 'm/sec', [-28 28], 5, 5,'idl',1);
    % title('Velocity', 'Color', 'k'); grid on
    % set(gcf, 'Position', get(0,'Screensize'));
    % xlim([-40 40]);
    % ylim([-40 40]);
    % fig_name= sprintf('Velocity2');
    % saveas(ff,fullfile(path_f,fig_name),'png');
    % 
    % figure; 
    % gg=plotPPI(range, azi', double(data.SNR), '', [0 15], 1, 5,1);
    % title('SNR', 'Color', 'k'); grid on
    % set(gcf, 'Position', get(0,'Screensize'));
    % xlim([-40 40]);
    % ylim([-40 40]);
    % fig_name= sprintf('SNR');
    % saveas(gg,fullfile(path_f,fig_name),'png');
    % 
    % close all
    % 
end