load("C:\Users\Jesse Li\TU Delft Year 1\Q2\msd\dat\part2_B7_C_tf_fit.mat");

%%
C_sys = frd(C_tf_fit, omega_fit);
% C_sys.TimeUnit = 'seconds'; % default
% C_sys.FrequencyUnit = 'rad/TimeUnit's; % default

%%
n_poles = 4;
n_zeros = 3;
opt = tfestOptions('EnforceStability', true, 'InitializeMethod', 'all');
C_est = tfest(C_sys, n_poles, n_zeros, opt);
bode(C_sys, C_est);
% Compare the estimated transfer function with the original
legend('Original', 'AutoEstimated');
grid on;

% good enough lol
%%  Save results
num = C_est.Numerator;
den = C_est.Denominator;
save('dat/part2_B7_C_tf_fit.done.mat');

% %%
% % Manually guess
% Kp = 1500;
% wd = 750 * 2 * pi;
% zeta_d = 0.12;
% wl = 35000;
% zeta_l = 0.707;
% phi_lag = 20 * pi / 180;
% omega_lag = 700;
% s = tf('s');
% 
% alpha_lag = (1 - sin(phi_lag)) / (1 + sin(phi_lag));
% T_lag = 1 / (omega_lag * sqrt(alpha_lag));
% 
% integrator = tf([Kp], [1, 0]);
% double_zero = tf([1, 2 * zeta_d * wd, wd^2], wd^2);
% lpf = tf(wl^2, [1, 2 * zeta_l * wl, wl^2]);
% lag = tf([alpha_lag * T_lag, 1], [T_lag, 1]);
% 
% manual_guess = integrator * double_zero * lpf * lag;
% 
% % % Manually guess, assuming notch filter
% % Kp = 2000;
% % omega_n = 750 * 2 * pi; % notch filter frequency
% % q2_n = 10;
% % gain_n = 2;
% % wl = 40000;
% % zeta_l = 0.707;
% % phi_lag = 30 * pi / 180;
% % omega_lag = 1000;
% % s = tf('s');
% % 
% % alpha_lag = (1 - sin(phi_lag)) / (1 + sin(phi_lag));
% % T_lag = 1 / (omega_lag * sqrt(alpha_lag));
% % q1_n = gain_n * q2_n;
% % 
% % integrator = tf([Kp], [1, 0]);
% % % double_zero = tf([1, 2 * zeta_d * wd, wd^2], wd^2);
% % lpf = tf(wl^2, [1, 2 * zeta_l * wl, wl^2]);
% % notch = tf([1, omega_n / q1_n, omega_n^2], [1, omega_n / q2_n, omega_n^2]);
% % lag = tf([alpha_lag * T_lag, 1], [T_lag, 1]);
% % 
% % manual_guess = integrator * notch;
% 
% 
% 
% bode(C_sys, C_est, manual_guess);
% % Compare the estimated transfer function with the original
% legend('Original', 'AutoEstimated', 'Manual');
% grid on;
% 
% %%
% 
% % Manually initialize
% opt = tfestOptions('EnforceStability', true);
% C_est2 = tfest(C_sys, manual_guess);
% 
% bode(C_sys, C_est, C_est2, manual_guess);
% % Compare the estimated transfer function with the original
% legend('Original', 'AutoEstimated', 'ManInit', 'Manual');
% grid on;
