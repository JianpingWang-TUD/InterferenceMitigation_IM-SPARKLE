flag_plot=0;
path = './figs/fig_20211108/';
%% FMCW
c = 3e8;
f_c = 3e9;
T_sw = 400e-6;
BW = 40e6;
sweep_slope = BW/T_sw;
Range_max = 8e3;          % max detection range
tau_max = 2*Range_max/c;   % max time delay
fb_max = sweep_slope * tau_max; % max beat frequency

BW_I = BW;
T_sw_I = T_sw;
sweep_slope_I = BW_I/T_sw_I;

f_s = 12e6;                 % sampling frequency
Npoint = floor(f_s*T_sw);

P_tx = 1;
SNR = 15; 
% SNR = inf;
%% Transmitter
t = (0 : 1 : Npoint-1)/f_s;
N_t = length(t);
amp     = sqrt(P_tx);
% sig_Tx  =  FMCW_sweep(t, T_sw, 0, f_c, -sweep_slope, amp);
%% Target
d_tar = [2e3,  3.5e3,  5e3 ];
N_tar = length(d_tar);
scat_coeff_tar = [0.5, 0.05, 0.4] .* exp( 1i*2*pi*rand(1, N_tar) );

sig_Rx = beatSig_FMCW(scat_coeff_tar(1:3), d_tar(1:3), t, f_c, T_sw, sweep_slope, c);
sigN_rx = awgn(sig_Rx,SNR,'measured');
%% FMCW Interference
amp_intf = [5, 3.6, 4.5];  %5, 3.6, 4.5
fr_intf = [3*sweep_slope_I, -2*sweep_slope_I, -1.5*sweep_slope_I ]; %3 -2 -1.5
fc_intf = [f_c, f_c, 1*f_c];
T_sw_intf = [T_sw_I, 1*T_sw_I, 1.2*T_sw_I];
t_d_intf = [10e-6, 150e-6,  -160e-6];

sig_I = beatInterfer_FMCW(amp_intf, fc_intf, fr_intf, T_sw_intf, t_d_intf,...
    t, f_c, sweep_slope, T_sw, fb_max);

%% Full signal (usable signal + interference + noise)
% sig_full = sigN_rx + sig_I;

%% Time truncation 
ind_t_eff = rectpuls(t-tau_max-T_sw/2, T_sw)>0.5;
t = t(ind_t_eff);
t_us = t*1e6;

sig_Rx_trc = sig_Rx(ind_t_eff);
sig_I_trc = sig_I(ind_t_eff);
% sig_full_trc = sig_full(ind_t_eff);

%====================================    
sigN_rx_trc = awgn(sig_Rx_trc,SNR,'measured');
% sigN_rx_trc = sig_Rx_trc;
sig_full_trc = sigN_rx_trc + sig_I_trc;

%=======================================================
len_sig = length(sig_full_trc);

%===================================
% figure
% plot(t_us,abs(sig_I_trc))
% figure
% plot(t_us,real(sig_full_trc))
len_interf = sum(abs(sig_I_trc)>0);
dur_pct_interf = len_interf/len_sig
%===================================

save('../data/Data4Demo.mat', 'sig_full_trc', 'sig_Rx_trc', 't', ...
    'f_c', 'f_s', 'sweep_slope', 'c')