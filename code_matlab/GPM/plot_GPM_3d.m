% lat_swath = lat_ku(swath_number,:);
% lon_swath = lon_ku(swath_number,:);
% Lat_lon_index = -32.5 < lat_ku & lat_ku < -31.5 & -66 < lon_ku & lon_ku < -64;
% ZcKu_RHI = ZcKu(65:176,  Lat_lon_index);
% [r,~,c] = size(ZcKu_RHI);
% ZcKu_RHI = reshape(ZcKu_RHI, r, c);
% lat_rhi = lat_swath(Lat_lon_index);
% lon_rhi = lon_swath(Lat_lon_index);
% z = (r:-1:1)*0.125;
% figure


%%
ZcKu(ZcKu<15) = nan;
ZcKu_volume = ZcKu(60:176,:,1170:1210);
ZcKu_volume = permute(ZcKu_volume, [2 3 1]);
figure
cross_track = 1:49;
along_track = 1170:1210;
height = 0.125*(117:-1:1);
[x, y, z] = meshgrid(along_track, cross_track, height);

xslice = [1177 1198 1207];
yslice = [23 37];
zslice = [1];%[10 20 30 40 50 60 70 80 90 110 110];

slice(x,y,z, ZcKu_volume ,xslice,yslice,zslice)
xlabel('Along track')
ylabel('Cross track')
zlabel('Height (km)')
shading flat
colormap jet
set(gca, 'YDir','reverse')
%set(gca, 'XDir','reverse')
return 
%%
figure
x = -2:.2:2;
y = -2:.25:2;
z = -2:.16:2;

[x,y,z] = meshgrid(x,y,z);
v = x.*exp(-x.^2-y.^2-z.^2);

xslice = [-1.2,.8,2];    % location of y-z planes
yslice = 2;              % location of x-z plane
zslice = [-2,0];         % location of x-y planes

slice(x,y,z,v,xslice,yslice,zslice)
xlabel('x')
ylabel('y')
zlabel('z')
