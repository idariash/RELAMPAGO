% clear all
% Scattering_Light_Rain_Drops
%density plot
close all
%addpath /net/denali/storage2/radar2/tmp/Ivan/Utiles/Matlab
[A,CZplot,CZDRplot]=p_icp(Z, ZDR,400);
%A(A<-10)=NaN;
mean_A = nanmean(A,2);
figure; pcolor(CZplot, CZDRplot, log10(A)); shading flat; grid on; hold on
caxis([1 3]);
%%
xlim([10 35])
ylim([-2 2])
%%
my_handle=colorbar('ytick',2:0.2:5.4,'FontWeight','bold');
set(get(my_handle,'Title'),'string','log_{10}(N_{obs})','FontWeight','bold');
axis([10 40 -4 4]);
set(gcf, 'Position', [100, 100, 800, 800])
set(gca,'XTick',[0:1:50],'FontSize',20,'FontWeight','bold');
set(gca,'YTick',[-2:0.2:2],'FontSize',20,'FontWeight','bold');
xlabel(' Z (dBZ)','FontSize',20,'FontWeight','bold');
ylabel(' ZDR (dB)','FontSize',20,'FontWeight','bold');
%%

%errorbar
 X = Z;
 Y = ZDR;
inter_Zh=min(X):1:max(X);
%inter_Zh=10:2:34;
for i=2:length(inter_Zh)
    ss=find(X>inter_Zh(i-1) & X<inter_Zh(i) );
    numberOfElement = sum(ss);
    x_mean(i)=(inter_Zh(i-1)+inter_Zh(i))/2;
    int_n(i)=length(ss);
    y_mean(i)=nanmean(Y(ss));
    y_std(i)=nanstd(Y(ss))/sqrt(numberOfElement);
    disp(y_std(i))
end
hold on;errorbar(x_mean(2:end),y_mean(2:end),y_std(2:end),'-kx','linewidth',1.5); grid on;
% figure
% histogram(Z)
% title('Histogram Z (dBZ)')
% figure
% hist(ZDR,64)
% title('Histogram ZDR (dB)')
% figure
% errorbar(x_mean(2:end),y_mean(2:end),y_std(2:end),'-kx','linewidth',1.5); grid on;
y_mean(x_mean<0) = nan;
y_mean(x_mean>10) = nan;
disp(nanmean(y_mean))