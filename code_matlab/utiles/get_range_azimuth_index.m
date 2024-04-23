% Ivan Arias
% 2020/01/14
% get the Range and Azimuth and Elevation Indexation

function [Range, Azimuth, Elevation] = get_range_azimuth_index(range, azimuth, elevation, ray_n_gates)

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
    Azimuth = [Azimuth azimuth_n_gates];
end