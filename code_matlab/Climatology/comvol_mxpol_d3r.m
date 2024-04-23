
figure
%D3R
lineLength              = 150;
angle1                  = 1.4;                                            % change elevation
angle2                  = angle1 + 1;                                       % change beamwidth
x(1)                    = 0;                                                %change x and y accordingly
y(1)                    = 0.480;
x(2)                    = x(1) + lineLength * cosd(angle1);
y(2)                    = y(1) + lineLength * sind(angle1);
a(1)                    = 0;
b(1)                    = 0.480;                                            %change x and y accordingly
a(2)                    = a(1) + lineLength * cosd(angle2);
b(2)                    = b(1) + lineLength * sind(angle2);
hold on;
plot(x, y,'b','LineWidth',2);
plot(a, b,'b','LineWidth',2);
xlim([-50 80]);
ylim([0 5]);

%JCO
lineLength              = 150;
angle1                  = 180 - 0.48;                                       % change elevation
angle2                  = angle1 -1;                                  % change beam width
x(1)                    = 79;                                            %change x and y accordingly
y(1)                    = 1.1;
x(2)                    = x(1) + lineLength * cosd(angle1);
y(2)                    = y(1) + lineLength * sind(angle1);
a(1)                    = 79;                                             %change a and b accordingly
b(1)                    = 1.1;
a(2)                    = a(1) + lineLength * cosd(angle2);
b(2)                    = b(1) + lineLength * sind(angle2);
plot(x, y,'r','LineWidth',2);
plot(a, b,'r','LineWidth',2);
xlim([0 80]);
ylim([0 5]);


xlabel('Distance from CHIVO (km)','FontWeight', 'normal','FontSize', 14);
ylabel('Height (km)','FontWeight', 'normal','FontSize', 14);
title('Common observation volume','FontWeight', 'normal','FontSize', 14);

start_x=32;
end_x=48;
y_range=[-25 25];
hold on;
plot([start_x,start_x],y_range,'k','LineWidth',1);
plot([end_x,end_x],y_range,'k','LineWidth',1);