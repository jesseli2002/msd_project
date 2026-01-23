%% Define Variables
load('MSD2025_P3_Plant.mat')
s = tf('s');

load('MSD2025_P3_Plant.mat')
load('MSD2025_P3_Signals.mat')

ts = 30e-6;
fs = 1.0/ts;
L = length(d);
% t = ts * (num_samples-1);
t = (0 : L-1)*ts;
f = fs*(0:(L/2))/L;

C1 = 2750 / s; % integrator
C2 = (1 / make_double_pole(4610, 0.01)) * make_double_pole(8000, 1);
C3 = make_notch(6300, 5);
C4 = make_lpf_butter(28000, 2);
C5 = make_ll(300*2*pi, 20*pi/180);

C = ss(C1) * ss(C2) * ss(C3) * ss(C4) * ss(C5);
C = balreal(C);
C_tf = tf(C);

%% A.3 - Add an anti-notch

C6 = make_anti_notch(w_a, 10, 5);
C = C * ss(C6);

%% A.4 - Add an integrator
C7 = make_PI(5,2*pi*200);
C_adj = C * ss(C7);
Loop_adj = C_adj * G;
S_adj = feedback(1, Loop_adj);
T_adj = feedback(Loop_adj, Loop_adj);

%% Define loop
Loop = C * G;
S = feedback(1, Loop);   % 1 / (1 + Loop)
T = feedback(Loop, Loop);   % Loop / (1 + Loop)

% H_d = -T;
H_d = ss(G)*S;   % y / d
H_n = S;       % y / n

%% A.1 - Compute error
w_a = 2*pi*10;

S_mag = abs(freqresp(S, w_a));
T_mag = abs(freqresp(T, w_a));
PS_mag = abs(freqresp(G*S, w_a));

r_a = 10;          % um
d_a = 20;          % V
n_a = 0.01;        % 10nm = 0.1um

% e_real = sqrt((S_mag * r)^2 + (PS_mag * d)^2 + (T_mag * n)^2);
% fprintf('e_real = %.8f um \n', e_real)

e_r = S_mag * r_a;
e_d = PS_mag * d_a;
e_n = T_mag * n_a;
e_total = norm([e_r, e_d, e_n]);
fprintf('e_total = %.8f um \n', e_total)

%% B.1 Plot time domain signals
figure(1);
plot(t, d, 'b'); hold on;
plot(t, n, 'r');
xlabel('Time [s]');
ylabel('Amplitude');
legend('Disturbance d', 'Noise n');
grid on;

%% B.2 fft and such

D = fft(d);
N = fft(n);

f_b = f(2:end);
% f_b = f;
w_b = 2*pi*f_b;   % rad/s

D_S2 = (1/(fs*L)) * abs(D).^2;
D_S1 = D_S2(1:L/2+1);
D_S1(2:end-1) = 2*D_S1(2:end-1);

N_S2 = (1/(fs*L)) * abs(N).^2;
N_S1 = N_S2(1:L/2+1);
N_S1(2:end-1) = 2*N_S1(2:end-1);

% sanity check
d_var_time = mean(d.^2);
d_var_freq = sum(D_S1) * (fs/L);
n_var_time = mean(n.^2);
n_var_freq = sum(N_S1) * (fs/L);

D_S1 = D_S1(2:end);
N_S1 = N_S1(2:end);
Hd_mag = abs(squeeze(freqresp(H_d, w_b)));
Hn_mag = abs(squeeze(freqresp(H_n, w_b)));

% force everything into 1D column
D_S1   = D_S1(:);
N_S1   = N_S1(:);
Hd_mag = Hd_mag(:);
Hn_mag = Hn_mag(:);
f_b    = f_b(:);
w_b    = w_b(:);

PSD_d = (Hd_mag.^2) .* D_S1;   % (um^2/Hz)
PSD_n = (Hn_mag.^2) .* N_S1;   % (um^2/Hz)

PSD_y = PSD_d + PSD_n;

CPS_d = cumtrapz(f_b, PSD_d);
CPS_n = cumtrapz(f_b, PSD_n);
CPS_y = cumtrapz(f_b, PSD_y);


%% B.2 plot
figure(2);
loglog(f_b, D_S1); hold on;
loglog(f_b, N_S1);
xlabel('Frequency [Hz]');
ylabel('PSD [units^2/Hz]');
legend('d','n');
grid on;

%% B.3 plot (PSD y)
figure(3);
clf(3);
loglog(f_b, PSD_d, 'b', 'LineWidth', 1.5); hold on;
loglog(f_b, PSD_n, 'r', 'LineWidth', 1.5);
loglog(f_b, PSD_y, 'k--', 'LineWidth', 2);

xlabel('Frequency [Hz]');
ylabel('Output PSD S_y(f) [\mum^2/Hz]');
legend('From disturbance d', 'From noise n', 'Total');
grid on;
title('Output Displacement PSD');

%% B.4 plot (CPS y)
figure(4);
clf(4);
semilogx(f_b, CPS_d, 'b', 'LineWidth', 1.5); hold on;
semilogx(f_b, CPS_n, 'r', 'LineWidth', 1.5);
semilogx(f_b, CPS_y, 'k--', 'LineWidth', 2);

xlabel('Frequency [Hz]');
ylabel('Cumulative Power [\mum^2]');
legend('From disturbance d', 'From noise n', 'Total');
grid on;
title('Cumulative Power Spectrum (CPS)');


%% A.1 Bode plot (magnitude + phase) in one figure

% Frequency vector: 1 Hz to 20 kHz
f = logspace(0, log10(2e4), 500);  % Hz
w = 2*pi*f;                        % rad/s

% Ensure all systems are state-space (numerical stability)
G_ss = ss(G);      % plant
C_ss = C;          % controller
L_ss = C_ss * G_ss;

% Compute Bode responses
[magG, phaseG] = bode(G_ss, w);
[magC, phaseC] = bode(C_ss, w);
[magL, phaseL] = bode(L_ss, w);

% Squeeze to 1D arrays
magG   = squeeze(magG);   
phaseG = squeeze(phaseG);
magC   = squeeze(magC);   
phaseC = squeeze(phaseC);
magL   = squeeze(magL);   
phaseL = squeeze(phaseL);

% Plot magnitude and phase in subplots
figure;

% Magnitude subplot
subplot(2,1,1)
semilogx(f, 20*log10(magG), 'b', 'LineWidth', 1.5); hold on;
semilogx(f, 20*log10(magC), 'r', 'LineWidth', 1.5);
semilogx(f, 20*log10(magL), 'k', 'LineWidth', 2);
grid on;
ylabel('Magnitude [dB]');
title('Bode Plot: Plant G, Controller C, Open-loop L = C*G');
legend('G', 'C', 'L');

% Phase subplot
subplot(2,1,2)
semilogx(f, phaseG, 'b', 'LineWidth', 1.5); hold on;
semilogx(f, phaseC, 'r', 'LineWidth', 1.5);
semilogx(f, phaseL, 'k', 'LineWidth', 2);
grid on;
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
legend('G', 'C', 'L');


%% Functions

function C = make_double_pole(omega_n, zeta)
    C = tf(omega_n^2, [1, 2*zeta*omega_n, omega_n^2]);
end


function C = make_lpf_butter(omega_lpf, order)
    if nargin < 2
        order = 2;
    end
    [num, den] = butter(order, omega_lpf, 's');
    C = tf(num, den);
end

function C = make_notch(omega_n, gain, Q2)
%MAKE_NOTCH Create a second-order notch filter
%
%   omega_n : center frequency [rad/s]
%   gain    : depth/height of notch (default = 10)
%   Q2      : pole Q factor (controls width, default = 10)
%
%   Higher Q2  -> sharper notch
%   Higher gain -> deeper notch
    if nargin < 2
        gain = 10;
    end
    if nargin < 3
        Q2 = 10;
    end

    Q1 = gain * Q2;
    num = [1, omega_n / Q1, omega_n^2];
    den = [1, omega_n / Q2, omega_n^2];
    C = tf(num, den);
end

function C_anti = make_anti_notch(omega_n, gain, Q2)
%MAKE_ANTI_NOTCH Create an anti-notch filter to boost gain at omega_n
%
%   omega_n : center frequency [rad/s]
%   gain    : depth/height of anti-notch (default = 10)
%   Q2      : pole Q factor (controls width, default = 10)
%
%   Higher gain -> more boost at omega_n
%   Higher Q2   -> narrower frequency boost
    if nargin < 2
        gain = 10;
    end
    if nargin < 3
        Q2 = 10;
    end

    Q1 = gain * Q2;

    % Swap numerator & denominator of normal notch
    num = [1, omega_n / Q2, omega_n^2];
    den = [1, omega_n / Q1, omega_n^2];

    C_anti = tf(num, den);
end


function C = make_ll(omega_ll, phi_m)
    alpha = (1 - sin(phi_m)) / (1 + sin(phi_m));
    T = 1 / (omega_ll * sqrt(alpha));
    C = tf([T, 1], [alpha*T, 1]);
end

function Cpi = make_PI(Kp, wi)
%MAKE_PI Create a PI controller
%
%   Kp : proportional gain
%   wi : integrator corner frequency [rad/s]

    if nargin < 2
        error('Must provide both Kp and wi');
    end

    % PI transfer function
    Cpi = tf([Kp, Kp*wi], [1, 0]);   % Kp + Kp*wi/s
end
