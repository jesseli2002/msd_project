%% 
load("dat/MSD2025_P2_Plant.mat");

%% P
s = tf('s');

wb = 3000; % Cut-off frequency in rad/s
A = 1e-4;  % Desired disturbance attenuation at LF
M = 2;     % Max sensitivity peak
Wp = (s/M + wb) / (s + wb*A);

% Define weight for complementary sensitivity
M = 1.5;
omega_t = 6000;
omega_corner = omega_t / M; % corner frequency
WT_inv = M * omega_corner / (s + omega_corner);
WT = 1 / WT_inv;
Wu = tf(1/10);

G.u = {'u'};
G.y = {'y_plant'};

Wp.u = {'y'};
Wp.y = {'z_p'};

Sum = sumblk('y = y_plant + y_dist');

WT.u = {'y_plant'};
WT.y = {'z_T'};

Wu.u = {'u'};
Wu.y = {'z_u'};

blocks = {G, Wp, WT, Wu, Sum};
inputs = {'y_dist', 'u'};
outputs = {'z_p', 'z_T', 'z_u', 'y'};

P = connect(blocks{:}, inputs, outputs);
P = minreal(P);

%% hinfsyn 
[K, CL, gamma, info] = hinfsyn(P, 1, 1);
K = -K;

L = G * K;
S = 1 / (1 + L);
T = L / (1 + L);
%%

% Open loop behaviour
figure;
bode(K, G, L);
legend('K', 'G', 'L');
title("Open loop");

% Sensitivity
figure;
bode(S, 1/Wp);
legend('S', '1/Wp');
title('Sensitivity');

% Complementary sensitivity
figure;
bode(T, 1/WT);
legend('T', '1/WT');
title('Complementary sensitivity');

% Control Sensitivity
figure;
bode(K * S, 1/Wu);
legend('KS', '1/Wu');
title('Control sensitivity');

% Margins
figure;
margin(L);

gamma
