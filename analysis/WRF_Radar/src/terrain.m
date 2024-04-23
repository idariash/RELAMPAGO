
filename = '/net/k2/storage/people/idariash/home/CSU/RELAMPAGO/analysis/Mountain/terr.nc';

terr = ncread(filename, 'topo')';
longitude = ncread(filename, 'X');
latitude = ncread(filename, 'Y');

[longitudes, latitudes] = meshgrid(longitude, latitude);

chivo_lat = -31.6342;
chivo_lon = -64.1686;

mask_SdC = abs(latitudes - chivo_lat) < 1 & abs(longitudes - chivo_lon) < 4;
terr(~mask_SdC) = nan;

smoth_filter = ones(10,10);
smoth_filter = 1/sum(sum(smoth_filter))*smoth_filter;

terr = conv2(terr, smoth_filter, 'same');

terr_max = max(terr)/1000;
longitude_SdC = longitude(~isnan(terr_max));
terr_SdC = terr_max(~isnan(terr_max));
distance = 111.320*cosd(chivo_lat)*(longitude_SdC - chivo_lon);


hold on
yyaxis right
plot(distance, terr_SdC,'r', 'LineWidth', 2)
xlim([-100, 100])
colormap(cmocean('rain'))