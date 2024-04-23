% 2019/11/07
% Ivan Arias
% Kdp vs Z dispersion to find funy regions in the cloud boundary

addpath /net/denali/storage2/radar2/tmp/Ivan/Utiles/Matlab
datapath = '/net/denali/storage/radar/RELAMPAGO/DROPS/2018/12/14/';

directory = dir([datapah, '*.nc']);

Z_total = [];
Kdp_total = [];
Height_total = [];

for i = 1:length(directory)
    filename = [datapath, directory(i).name];
    Z = ncread(filename, 'corrected_reflectivity');
    rhohv = ncread(filename, 'corrected_cross_correlation_ratio');
    Kdp = ncread(filename, 'corrected_specific_differential_phase');
    range = double(ncread(filename, 'range'))'/1e3;
    azimuth = ncread(filename, 'azimuth');
    elevation = ncread(filename, 'elevation');
    ray_n_gates = ncread(filename, 'ray_n_gates');
    
    Range = [];
    Elevation = [];
    Azimuth = [];
    for I = 1:length(ray_n_gates)
        Range = [Range range(1:ray_n_gates(I))];

    end

    for I = 1:length(elevation)
        elevation_n_gates = elevation(I)*ones(1,ray_n_gates(I));
        azimuth_n_gates = azimuth(I)*ones(1,ray_n_gates(I));
        Elevation = [Elevation elevation_n_gates];
        %Azimuth = [Azimuth azimuth_n_gates];
    end
    
    index_toFilter = rhohv > 0.95 & Z > -10 & Kdp > -5 & Range' < 100;
    
    Height = Range.*sind(Elevation);
    
    Z = Z(index_toFilter);
    Kdp = Kdp(index_toFilter);
    Height = Height(index_toFilter);
    
    
    Z_total = [Z_total Z'];
    Kdp_total = [Kdp_total Kdp']; 
    Height_total = [Height_total Height];
    disp(i)
    
end
return 
%% 

Height_high_Kdp = Height_total(Kdp_total > 4);
Ref_high_Kdp = Z_total(Kdp_total > 4);
figure
histogram(Height_high_Kdp, 16)

%%
[A,CZplot,CZDRplot]=p_icp(Kdp_total, Z_total, 50);
mean_A = nanmean(A,2);
figure; pcolor(CZplot, CZDRplot, log10(A)); shading flat; grid on; hold on
caxis([0 4]);
% hold on 
% plot(x,x, 'r')
hc=colorbar;
xlabel(hc,'log_{10}(N_{obs})');
xlim([-10 30])
ylim([-10 80])
