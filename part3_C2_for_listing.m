%% Load plant
load('MSD2025_P3_Plant.mat')
delay = G.InputDelay;

%% Design feedback controller
C1 = 1770 / tf('s');
C2 = make_double_pole(8000, 1) / make_double_pole(4610, 0.01);

% feedback controller
C_feedback = C1 * C2; 

%% Design feedforward controller
plant_dc = dcgain(G);

% Plant transfer poles and zeros were decomposed in separate preprocessing step;
% these are just the numerical results

% Inverse of plant poles (denominator)
G_den_tf0 = tf([1, 8.3578931e1, 2.15599929e+07],[21559992.9005093]);
G_den_tf1 = tf([1, 1.39059457e+02, 3.99535798e+07], [39953579.83864359]);
G_den_tf2 = tf([1, 8.65800866e+03], [8658.00865801]);

% Invert plant zero (numerator), with extra damping:
omega_zero = 6107.256118578557;
G_num_maxsin = make_double_pole(omega_zero, 0.016); % Sin gain limited
G_num_worst = make_double_pole(omega_zero, 0.023); % Worst case limited

% Define low-pass filters
% Sinusoidal gain limited
[lpf_b, lpf_a] = butter(4, 1000 * 2 * pi, 's');
lpf_maxsin = tf(lpf_b, lpf_a);

% Worst case limited
[lpf_b, lpf_a] = butter(4, 700 * 2 * pi, 's');
lpf_worst = tf(lpf_b, lpf_a);

C_feedforward_maxsin = 1 / plant_dc * G_den_tf0 * G_den_tf1 ...
    * G_den_tf2 * G_num_maxsin * lpf_maxsin;

% Worst case limited
C_feedforward_worst = 1 / plant_dc * G_den_tf0 * G_den_tf1 ...
    * G_den_tf2 * G_num_worst * lpf_worst;

%% Connect all elements together

% prefilter is a delay of same length as plant, times the low pass 
%   filter used in feedforward
% Adding the low pass filter significantly reduces overshoot in a step 
% response, and causes the r->y transfer function to simplify to just F(s)
% (where F is the low pass filter)
prefilter_maxsin = exp(-delay * tf('s')) * lpf_maxsin; 
prefilter_worst = exp(-delay * tf('s')) * lpf_worst; 

% Define input/output names for `connect`
C_feedforward_maxsin.InputName = 'ra';
C_feedforward_maxsin.OutputName = 'u_ff';
C_feedforward_worst.InputName = 'ra';
C_feedforward_worst.OutputName = 'u_ff';

prefilter_maxsin.InputName = 'ra';
prefilter_maxsin.OutputName = 'r';
prefilter_worst.InputName = 'ra';
prefilter_worst.OutputName = 'r';

C_feedback.InputName = 'e';
C_feedback.OutputName = 'u_fb';

G.InputName = 'u';
G.OutputName = 'y';

sum_u = sumblk("u = u_fb + u_ff + d");
sum_e = sumblk("e = r - y");

inputs = {'ra', 'd'};
outputs = {'e', 'y'};

% [C]losed [L]oop, [f]eed[f]orward + [f]eed[b]ack
cl_fffb_maxsin = connect(C_feedforward_maxsin, prefilter_maxsin, ...
    C_feedback, G, sum_u, sum_e, inputs, outputs);
cl_fffb_worst = connect(C_feedforward_worst, prefilter_worst, ...
    C_feedback, G, sum_u, sum_e, inputs, outputs);

% [C]losed [L]oop,  [f]eed[b]ack only
sum_u = sumblk('u = u_fb + d');
inputs = {'r', 'd'};
cl_fb = connect(C_feedback, G, sum_u, sum_e, inputs, outputs);

%%
% Copy transfer function and delete input/output name to remove spurious
% subtitle in plot. Setting subtitle() doesn't work.

cl_fb_copy = cl_fb;
cl_fb_copy.InputName = ''; 
cl_fb_copy.OutputName = ''; 
cl_fffb_maxsin_copy = cl_fffb_maxsin;
cl_fffb_maxsin_copy.InputName = ''; 
cl_fffb_maxsin_copy.OutputName = ''; 
cl_fffb_worst_copy = cl_fffb_worst;
cl_fffb_worst_copy.InputName = ''; 
cl_fffb_worst_copy.OutputName = ''; 
step(cl_fb_copy(2, 1), cl_fffb_maxsin_copy(2, 1), cl_fffb_worst_copy(2, 1))
legend("FB only", "FB + FF sinusoidal gain limited", ...
    "FB + FF worst case limited", "Location", "southeast")
xlim([0, 0.01])
% subtitle('')
title('Step repsonse')

%% Get step info (printed to consol)
cl_fb_stepinfo = stepinfo(cl_fb(2, 1))
cl_fffb_maxsin_stepinfo = stepinfo(cl_fffb_maxsin(2, 1))
cl_fffb_worst_stepinfo = stepinfo(cl_fffb_worst(2, 1))

%% Get process sensitivity bode plot
% Process sensitivity is given by y/d
figure;

bp = bodeplot(cl_fb_copy(2, 2));
bp.FrequencyUnit = 'Hz';
bp.MagnitudeUnit = 'abs';
bp.MagnitudeScale = 'log';
bp.XLimits = {[10, 1e4]};
bp.YLimits(1) = {[1e-4, 1e2]}; % gain plot 
bp.YLimits(2) = {[-720, 180]}; % phase plot
bp.Title.String = "";
% bp.Subtitle.String = "";
grid on;
title('Process sensitivity')


%% Get process impulse response
figure;
impulse(pade(cl_fb_copy(2,2), 6))
title('Impulse response')

%% Functions
function [sys] = make_double_pole(omega_n, zeta)
    sys = tf(omega_n^2, [1, 2 * zeta * omega_n, omega_n^2]);
end
