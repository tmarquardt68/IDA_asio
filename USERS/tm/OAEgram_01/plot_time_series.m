function plot_time_series(modDPlinesSeries,modSFlinesSeries)

a=1;
clf
% Plot modDP level
h_ax(a)=subplot('position',[0.0446    0.7700    0.9429    0.21]);  a=a+1;
plot(20*log10(abs(modDPlinesSeries(1,:))),'bo:'), hold on
plot(20*log10(abs(modDPlinesSeries(2,:))),'bx-')
plot(20*log10(abs(modDPlinesSeries(3,:))),'bx-','LineWidth',2),
plot(20*log10(abs(modDPlinesSeries(4,:))),'bx-')
plot(20*log10(abs(modDPlinesSeries(5,:))),'bo:')
plot(20*log10(abs(modDPlinesSeries(6,:))),'ro:')
plot(20*log10(abs(modDPlinesSeries(7,:))),'rx-')
plot(20*log10(abs(modDPlinesSeries(8,:))),'ro-','LineWidth',2),
plot(20*log10(abs(modDPlinesSeries(9,:))),'rx-')
plot(20*log10(abs(modDPlinesSeries(10,:))),'ro:')
grid on,title modDP,zoom on

% Plot modDP CM
h_ax(a)=subplot('position',[0.0446    0.5100    0.9429    0.21]); a=a+1;
plot(20*log10(abs(modDPlinesSeries(11,:))),'bo:'), hold on
plot(20*log10(abs(modDPlinesSeries(12,:))),'bx-')
plot(20*log10(abs(modDPlinesSeries(13,:))),'bx-','LineWidth',2),
plot(20*log10(abs(modDPlinesSeries(14,:))),'bx-')
plot(20*log10(abs(modDPlinesSeries(15,:))),'bo:')
plot(20*log10(abs(modDPlinesSeries(16,:))),'ro:')
plot(20*log10(abs(modDPlinesSeries(17,:))),'rx-')
plot(20*log10(abs(modDPlinesSeries(18,:))),'ro-','LineWidth',2),
plot(20*log10(abs(modDPlinesSeries(19,:))),'rx-')
plot(20*log10(abs(modDPlinesSeries(20,:))),'ro:')
grid on,title modDPCM,zoom on

% Plot modSF level
h_ax(a)=subplot('position',[0.0446    0.27    0.9429    0.200]);  a=a+1;
plot(20*log10(abs(modSFlinesSeries(1,:))),'ko:'), hold on
plot(20*log10(abs(modSFlinesSeries(2,:))),'kx-')
plot(20*log10(abs(modSFlinesSeries(3,:))),'kx-','LineWidth',2),
plot(20*log10(abs(modSFlinesSeries(4,:))),'kx-')
plot(20*log10(abs(modSFlinesSeries(5,:))),'ko:')
plot(20*log10(abs(modSFlinesSeries(6,:))),'k+-','LineWidth',2),
grid on,title modSF,zoom on

% Plot modCM_f2 level
h_ax(a)=subplot('position',[0.0446   0.03   0.9429    0.200]); 
plot(20*log10(abs(modSFlinesSeries(7,:))),'ko:'), hold on
plot(20*log10(abs(modSFlinesSeries(8,:))),'kx-')
plot(20*log10(abs(modSFlinesSeries(9,:))),'kx-','LineWidth',2),
plot(20*log10(abs(modSFlinesSeries(10,:))),'kx-')
plot(20*log10(abs(modSFlinesSeries(11,:))),'ko:')
grid on,title modCMf2,zoom on

xlabel('no of 6s-presentions'),zoom on
linkaxes(h_ax, 'x')
