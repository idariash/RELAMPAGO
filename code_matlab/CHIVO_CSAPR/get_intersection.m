% Compute intersection of RHI

Azimuth = input('Azimuth CHIVO ');
phi = Azimuth*pi/180;
x1 = -62.1378;
y1 = -54.6392;

%x = (35 - 35*tan(pi/3))/(tan(pi/3) - tan(pi/2 - phi))
%y = tan(pi/2 - phi)*x

x = (x1*tan(pi/3) -y1)/(tan(pi/3) - tan(pi/2 - phi))
y = tan(pi/2 - phi)*x

CHIVO = sqrt(x^2 +y^2)
CSAPR = sqrt((x - x1)^2 + (y - y1)^2)