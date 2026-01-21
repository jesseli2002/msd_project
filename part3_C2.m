%%
load('MSD2025_P3_Plant.mat')

num = G.Numerator;
den = G.Denominator;
delay = G.InputDelay;
save('dat/MSD2025_P3_Plant_numdendelay.mat', 'num', 'delay', 'den');

%%
C1 = 1770 / tf('s');
C2 = make_double_pole(8000, 1) / make_double_pole(4610, 0.01);

% feedback controller
C_feedback = C1 * C2; 

% feedforward controller
plant_dc = dcgain(G);

% Plant transfer poles and zeros were decomposed in separate preprocessing step; these are just the numerical results

% G_den_tf0 = tf([1, 8.3578931e1, 2.15599929e+07],[21559992.9005093]);
% G_den_tf1 = tf([1, 1.39059457e+02, 3.99535798e+07], [39953579.83864359]);
% G_den_tf2 = tf([1, 8.65800866e+03], [8658.00865801]);

% add damping
G_den_tf0 = tf([1, 8.3578931e1, 2.15599929e+07],[21559992.9005093]);
G_den_tf1 = tf([1, 1.39059457e+02, 3.99535798e+07], [39953579.83864359]);
G_den_tf2 = tf([1, 8.65800866e+03], [8658.00865801]);
omega_zero = 6107.256118578557;

omega_lpf = 1000 * 2 * pi;
[lpf_b, lpf_a] = butter(4, omega_lpf, 's');
lpf = tf(lpf_b, lpf_a);

C_feedforward = 1 / plant_dc * G_den_tf0 * G_den_tf1 * G_den_tf2 * make_double_pole(omega_zero, 0.016) * lpf;

bode(C_feedforward)

%%

% prefilter is a delay of same length as plant, times the low pass filter used in feedforward
% Adding the low pass filter significantly reduces overshoot in a step response, and causes the r->y transfer function to simplfiy to just F(s) (where F is the low pass filter)
% prefilter = exp(-delay * tf('s')); 
prefilter = exp(-delay * tf('s')) * lpf; 

% Prepare for connect
C_feedforward.InputName = 'ra';
C_feedforward.OutputName = 'u_ff';

prefilter.InputName = 'ra';
prefilter.OutputName = 'r';

C_feedback.InputName = 'e';
C_feedback.OutputName = 'u_fb';

G.InputName = 'u';
G.OutputName = 'y';

sum_u = sumblk("u = u_fb + u_ff + d");
sum_e = sumblk("e = r - y");

inputs = {'ra', 'd'};
outputs = {'e', 'y'};

% [C]losed [L]oop, [f]eed[f]orward + [f]eed[b]ack
cl_fffb = connect(C_feedforward, prefilter, C_feedback, G, sum_u, sum_e, inputs, outputs);

% [C]losed [L]oop,  [f]eed[b]ack only
sum_u = sumblk('u = u_fb + d');
inputs = {'r', 'd'};
cl_fb = connect(C_feedback, G, sum_u, sum_e, inputs, outputs);

%%
step(cl_fb(2, 1), cl_fffb(2, 1))
    

%%
function [sys] = make_double_pole(omega_n, zeta)
    sys = tf(omega_n^2, [1, 2 * zeta * omega_n, omega_n^2]);
end
