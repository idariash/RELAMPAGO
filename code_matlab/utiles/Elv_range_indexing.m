%Ivan Arias
% 2018/02/18
% Indexing Elevation and range for PPI RELAMPAGO

Range = [];
Elevation = [];
Azimuth = [];
for I = 1:length(ray_n_gates)
    Range = [Range range(1:ray_n_gates(I))];
    
end

for I = 1:length(elevation)
    elevation_n_gates = elevation(I)*ones(1,ray_n_gates(I));
    Elevation = [Elevation elevation_n_gates];
end

for I = 1:length(azimuth)
    azimuth_n_gates = azimuth(I)*ones(1,ray_n_gates(I));
    Azimuth = [Azimuth azimuth_n_gates];
end