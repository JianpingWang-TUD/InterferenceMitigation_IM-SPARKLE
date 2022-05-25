%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix completion                                    %
% Factorization nuclear norm to sum of Frobenious norm %
% Simulation: Three stationary targets                 %
% Edit by J.Wang, May 13, 2020                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear; clc;
addpath('./lowRaS_MC/')

path4figure = './figs/';
flag_plot=0;
%%  load data
load('./data/Data4Demo.mat')

t_us = t*1e6;   % change the unit of time to micro-second.
%% ======================Interference Mitigation =========================
%% =============IM-SPARKLE==========================
disp('IM-SPARKLE...')
beta_1 = 0.1;           % 0.5
mu     = 0.02;          % 0.05
tau    = 0.02;          % 0.02
k_beta = 1.6;
k_mu   = 1.2; 
R = 10;

tic;
[x_lowRaS, i_lowRaS, rerr] = lowRaS_Hankel(sig_full_trc, R, beta_1, mu, tau, k_beta, k_mu);
toc;
hf = figure
semilogy(1:length(rerr), rerr,'r-','linewidth',1);
xlabel('Iteration','FontSize', 11)
ylabel('Relative error','fontsize',11)
grid on
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])
if flag_plot==1
print(hf,'-dpng', '-r300', [path4figure 'simu_SNR_' num2str(SNR)  '_converg.png']);
print(hf,'-dpdf', '-r300', [path4figure 'simu_SNR_' num2str(SNR) '_converg.pdf'],'-opengl');
saveas(hf, [path4figure 'simu_SNR_' num2str(SNR)  '_converg.fig'])
end
%% Evalutation metric

SINR_0 = 20*log10(norm(sig_Rx_trc)/norm(sig_Rx_trc- sig_full_trc) )

SINR_lowRaS = 20*log10(norm(sig_Rx_trc)/norm(sig_Rx_trc.'- x_lowRaS) )
corr_lowRaS = (x_lowRaS)'*(sig_Rx_trc.')/(norm(x_lowRaS) * norm(sig_Rx_trc))


%% Image display
hfig_1 =figure;
Ftsz = 11;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])
plot(t_us, real(sig_full_trc), 'b-')%t_us, real(sig_Rx_trc), 'b-',...
hold on
%x_rect1 = [70 104 104 70 70];       y_rect1 = [-7 -7 7 7 -7];
x_rect1 = [101 152 152 101 101];       y_rect1 = [-7 -7 7 7 -7];
x_rect2 = [185 245 245 185 185];    y_rect2 = [-7.5 -7.5 7.5 7.5 -7.5];
x_rect3 = [280 320 320 280 280];    y_rect3 = [-6 -6 6 6 -6];
plot(x_rect1, y_rect1,'r--',...
     x_rect2, y_rect2,'r--',...
     x_rect3, y_rect3,'r--')
text(x_rect2(3),y_rect2(3)-0.2, '\leftarrow Interference','Color','red','Fontsize', Ftsz)
grid on
axis tight
ylim([-8 8])
xlabel('Time [\mus]', 'FontSize', Ftsz)
ylabel('Amplitude', 'FontSize', Ftsz)
title('Real part', 'FontSize', Ftsz)
%legend('reference', 'Contaminated Sig')
if flag_plot==1
print(hfig_1,'-dpng', '-r300', [path4figure 'simu_sig_full_SNR_' num2str(SNR)  '.png']);
print(hfig_1,'-dpdf', '-r300', [path4figure 'simu_sig_full_SNR_' num2str(SNR) '.pdf'],'-opengl');
saveas(hfig_1, [path4figure 'simu_sig_full_SNR_' num2str(SNR)  '.fig'])
end

hfig_2 =figure;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 4],'PaperSize',[6 4])
subplot(211)
plot(t_us, real(sig_Rx_trc),'-r',...
     t_us, real(x_lowRaS),'g-.')
% % hold on
% % x_rect_inset = [100 110 110 100 100]; y_rect_inset = [-2.7 -2.7 2.7 2.7 -2.7];
% % plot(x_rect_inset,y_rect_inset,'r--')
grid on
% axis tight
ylim([-1.8, 1.8])
legend('ref sig','Proposed'); %,'Location','south','Orientation','horizontal'
xlabel('Time [\mus]', 'FontSize', Ftsz)
ylabel('Amplitude', 'FontSize', Ftsz)
title('Recovered signal', 'FontSize', Ftsz)

%axes('Position',[0.45 0.78 0.45 0.12])
subplot(212)
II = t_us >=200 & t_us<=206;
% II = t_us >=250 & t_us<=255;
plot(t_us(II), real(sig_Rx_trc(II)),'r-',...
     t_us(II), real(x_lowRaS(II)), 'g-.')
grid on
axis tight
xlabel('Time [\mus]', 'Fontsize', Ftsz)
ylabel('Amplitude', 'Fontsize', Ftsz)
title('Close-up of the recovered signal')
if flag_plot==1
print(hfig_2,'-dpng', '-r300', [path4figure 'simu_extrac_usable_sig_SNR_' num2str(SNR) '.png']);
print(hfig_2,'-dpdf', '-r300', [path4figure 'simu_extrac_usable_sig_SNR_' num2str(SNR) '.pdf'],'-opengl');
saveas(hfig_2, [path4figure, 'simu_extrac_usable_sig_SNR_' num2str(SNR) '.fig'])
end
%%
hfig_3 =figure;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])
plot(t_us,real(x_lowRaS),'b--')
grid on
xlabel('Time [\mus]', 'FontSize', Ftsz)
ylabel('Amplitude', 'FontSize', Ftsz)
title('Extracted usable signal', 'FontSize', Ftsz)
if flag_plot==1
print(hfig_3,'-dpng', '-r300', [path4figure 'simu_usable_sig.png']);
print(hfig_3,'-dpdf', '-r300', [path4figure 'simu_usable_sig.pdf'],'-opengl');
end

%% Range profile
sig_len = length(sig_Rx_trc);
Num_fft_RP = 2^(nextpow2(sig_len)+1);

d = (0:Num_fft_RP-1)/Num_fft_RP * f_s /sweep_slope * c/2; 

RP_ref = ifft(sig_Rx_trc, Num_fft_RP);
RP_sig_full = ifft(sig_full_trc, Num_fft_RP);
RP_lowRaS = ifft(x_lowRaS, Num_fft_RP);



RP_max = max(abs([RP_ref, RP_sig_full, RP_lowRaS.']));

RP_ref_nor = db(abs(RP_ref)/RP_max);
RP_sig_full_nor = db(abs(RP_sig_full)/RP_max);
RP_lowRaS_nor = db(abs(RP_lowRaS)/RP_max);

Num_fft_disp = floor(Num_fft_RP/2);

d_1 = d(1:Num_fft_disp)/1e3;
hfig_5 = figure
plot(d_1, RP_ref_nor(1:Num_fft_disp), 'r-',...
     d_1, RP_sig_full_nor(1:Num_fft_disp), 'b--',...
     d_1, RP_lowRaS_nor(1:Num_fft_disp), 'g--')
grid on
axis([0 8 -60 2])
xlabel('Range [km]', 'Fontsize', Ftsz)
ylabel('Amplitude [dB]', 'Fontsize', Ftsz);
title('Range profiles', 'Fontsize', Ftsz);
legend('ref', 'sig\_Interf', 'Proposed', 'Location','southeast')

axes('Position', [0.68 0.72 0.2 0.18])
Ind = d_1>3.494 & d_1<3.506;
plot(d_1(Ind), RP_ref_nor(Ind), 'r-',...
     d_1(Ind), RP_sig_full_nor(Ind), 'b--',...
     d_1(Ind), RP_lowRaS_nor(Ind), 'g--')
grid on


axes('Position', [0.38 0.72 0.2 0.18])
Ind = d_1>1.994 & d_1<2.006;
plot(d_1(Ind), RP_ref_nor(Ind), 'r-',...
     d_1(Ind), RP_sig_full_nor(Ind), 'b--',...
     d_1(Ind), RP_lowRaS_nor(Ind), 'g--')
grid on

set(gcf, 'PaperUnits','inches', 'PaperSize', [4 3], 'PaperPosition', [0 0 4 3])
if flag_plot==1
    print(hfig_5, '-dpng', '-r300', [path4figure 'simu_PtTar_RP_SNR_' num2str(SNR) '.png']);
    print(hfig_5, '-dpdf', '-r300', [path4figure 'simu_PtTar_RP_SNR_' num2str(SNR) '.pdf']);
    saveas(hfig_5, [path4figure 'simu_PtTar_RP_SNR_' num2str(SNR) '.fig'])
end

%%
