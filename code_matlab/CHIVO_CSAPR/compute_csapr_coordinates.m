% Ivan Arias
% 2019/12/20
% Compute polar cordinates for CSAPR from CHIVO coordinates

function [r_csapr, azi_csapr, elv_csapr] = compute_csapr_coordinates(r_chivo, azi_chivo, elv_chivo)

r_csapr = sqrt((r_chivo*cosd(90 - azi_chivo) + 54)^2 + ...
    (r_chivo*sind(90 - azi_chivo) + 54)^2);

azi_csapr = 90 - atand((r_chivo*sind(90 - azi_chivo) + 54)/ ... 
    (r_chivo*cosd(90 - azi_chivo) + 54));

elv_csapr = atand((0.5 + r_chivo*tand(elv_chivo))/r_csapr);

