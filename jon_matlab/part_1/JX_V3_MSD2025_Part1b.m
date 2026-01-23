%% MSD 2025 Assignment Part 1: Open-loop System Identification
% Script prepared by: Jonathan Xie

%% VARIABLES
fmin = 10;   
fmax = 1e4;
ts = 30e-6;              
fs = 1/ts;
%fs = 10*fmax; % sampling frequency (in Hz)
%ts = 1.0/fs;
ft = logspace(max(-1,log10(fmin)),max(3,log10(fmax)),1000);

%% Single Chirp ===============
totaltime = 10; % in seconds
t = 0:ts:totaltime; % Time vector

tmax = t(end); % Time of chirp signal
method = 'logarithmic'; % or 'logarithmic'
u = chirp(t, fmin, tmax, fmax, method); % Generate a chirp signal
y = Plant(u, t); % Process the chirp signal through the system

%% Multi Chirp
totaltime = 40;      % total time duration
t = 0:ts:totaltime;

n_chirp = 4;                                 % number of chirps
chirp_time = totaltime/n_chirp;              % each chirp duration
t_single = 0:ts:chirp_time;                  % time for one chirp
method = 'logarithmic';
u_single = chirp(t_single, fmin, t_single(end), fmax, method);

u = repmat(u_single, 1, n_chirp);            % repeat chirps
u = u(1:length(t));                          % ensure same length
y = Plant(u, t); % Process the chirp signal through the system

%% Amped Chirp
f_inst = fmin * (fmax/fmin).^(t/tmax);  % logarithmic chirp
% Smooth amplitude ramp
f_boost_start = 2e3;   % 1 kHz
f_boost_end   = 1e4;   % 10 kHz
amp_min = 1;           % low-frequency amplitude
amp_max = 1e2;           % high-frequency amplitude

% Logarithmic ramp
amp = amp_min + (amp_max-amp_min) * log10(f_inst/f_boost_start)/log10(f_boost_end/f_boost_start);
amp(f_inst<f_boost_start) = amp_min;
amp(f_inst>f_boost_end) = amp_max;

u = u .* amp;
y = Plant(u, t); % Process the chirp signal through the system

%% PLOTTING CHIRP
figure(1); clf(1);

subplot(2,1,1)
plot(t, u, 'LineWidth', 1.0);
grid on;
title('x(t), Input Chirp Signal');
ylabel('Amplitude');

subplot(2,1,2)
plot(t, y, 'LineWidth', 1.0);
grid on;
title('y(t), Output Signal from Plant.p');
ylabel('Amplitude');
xlabel('Time (s)')

% export
width = 900
height = 500
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB1_chirp.png', 'Resolution', 300);

%% Define the white noise (realisitic with normal distribution) ===============
totaltime = 30; % in seconds
t = 0:ts:totaltime; % Time vector
vr = 10; % Define variance
u = sqrt(vr)*randn(length(t),1)'; % Generate white noise with variance = vr
y = Plant(u, t); % Process the chirp signal through the system

%% PLOTTING WHITE NOISE
figure(2); clf(2);

subplot(2,1,1)
plot(t, u, 'LineWidth', 1.0);
grid on;
title(sprintf('x(t), Input White Noise Signal, variance = %.1f',vr));
ylabel('Amplitude');

subplot(2,1,2)
plot(t, y, 'LineWidth', 1.0);
grid on;
title('y(t), Output Signal from Plant.p');
ylabel('Amplitude');
xlabel('Time (s)')

% export
width = 900;
height = 500;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB2_noise.png', 'Resolution', 300);

%% NEW RECT
input = u';
output = y;
% tfestimate settings
n_win = [4, 1];
percent_overlap = [0.5, 0.5];
L = length(u);
% ---- Preallocate ----
win_length   = zeros(size(n_win));
n_overlap    = zeros(size(n_win));
Tresp        = cell(size(n_win));
Cresp        = cell(size(n_win));
% ---- Compute estimates ----
for i = 1:length(n_win)

    % compute window size and overlap
    win_length(i) = round(L / n_win(i));
    n_overlap(i)  = floor(win_length(i)*percent_overlap(i));

    % compute transfer function + coherence
    [Tresp{i}, f] = tfestimate(input,output,rectwin(win_length(i)),n_overlap(i),ft,fs);
    [Cresp{i},~]  = mscohere(input,output,rectwin(win_length(i)),n_overlap(i),ft,fs);
end


%% NEW TRANSPARENT PLOTTING (RECT) — scalable version

figure(3); clf(3);
legendStrings = compose('%d windows', n_win);
line_thickness = 1.2;
last_thickness = 1.2;
% ---------- COLORS ----------
% Create N distinguishable colors
baseColors = lines(length(n_win));

% ---------- MAGNITUDE ----------
subplot(3,1,1); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        % 50% transparent
        semilogx(f, mag2db(abs(Tresp{i})), 'LineWidth', line_thickness, ...
                 'Color', [baseColors(i,:) 0.5]);
    else
        % Last one fully opaque
        semilogx(f, mag2db(abs(Tresp{i})), 'LineWidth', last_thickness, ...
                 'Color', [baseColors(i,:) 1.0]);
    end
end
set(gca,'XScale','log');
title('Frequency Response with different Rectangular Windows');
ylabel('|G| (dB)');

% ---------- PHASE ----------
subplot(3,1,2); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        semilogx(f, rad2deg(angle(Tresp{i})), 'LineWidth', line_thickness, ...
                 'Color', [baseColors(i,:) 0.5]);
    else
        semilogx(f, rad2deg(angle(Tresp{i})), 'LineWidth', last_thickness, ...
                 'Color', [baseColors(i,:) 1.0]);
    end
end
set(gca,'XScale','log');
ylabel('\angleG (deg)');

% ---------- COHERENCE ----------
subplot(3,1,3); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        semilogx(f, Cresp{i}, 'LineWidth', line_thickness, ...
                 'Color', [baseColors(i,:) 0.5]);
    else
        semilogx(f, Cresp{i}, 'LineWidth', last_thickness, ...
                 'Color', [baseColors(i,:) 1.0]);
    end
end
set(gca,'XScale','log');
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend(legendStrings, 'Location', 'southwest');

% ---------- EXPORT ----------
width = 350;
height = 380;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB_rect.png', 'Resolution', 300);

%% STORING MULTICHIRP
T_multichirp = T_hann;
C_multichirp = C_hann;

%% LOADING MULTICHIRP
T_hann = T_multichirp;
C_hann = C_multichirp;

%% NEW HANN
input = u';
output = y;
% tfestimate settings
n_win = [4, 1];
percent_overlap = [0.5, 0.5];
L = length(u);
% ---- Preallocate ----
win_length   = zeros(size(n_win));
n_overlap    = zeros(size(n_win));
T_hann        = cell(size(n_win));
C_hann        = cell(size(n_win));
% ---- Compute estimates ----
for i = 1:length(n_win)

    % compute window size and overlap
    win_length(i) = round(L / n_win(i));
    n_overlap(i)  = floor(win_length(i)*percent_overlap(i));

    % compute transfer function + coherence
    [T_hann{i}, f] = tfestimate(input,output,hann(win_length(i)),n_overlap(i),ft,fs);
    [C_hann{i},~]  = mscohere(input,output,hann(win_length(i)),n_overlap(i),ft,fs);
end

%% NEW TRANSPARENT PLOTTING (HANN) — scalable version

figure(4); clf(4);
legendStrings = compose('%d windows', n_win);
line_thickness = 1.2;
last_thickness = 1.2;
% ---------- COLORS ----------
% Create N distinguishable colors
baseColors = lines(length(n_win));

% ---------- MAGNITUDE ----------
subplot(3,1,1); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        % 50% transparent
        semilogx(f, mag2db(abs(T_hann{i})), 'LineWidth', line_thickness, ...
                 'Color', [baseColors(i,:) 0.5]);
    else
        % Last one fully opaque
        semilogx(f, mag2db(abs(T_hann{i})), 'LineWidth', last_thickness, ...
                 'Color', [baseColors(i,:) 0.5]);
    end
end
set(gca,'XScale','log');
title('Frequency Response with different Hann Windows');
ylabel('|G| (dB)');

% ---------- PHASE ----------
subplot(3,1,2); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        semilogx(f, rad2deg(angle(T_hann{i})), 'LineWidth', line_thickness, ...
                 'Color', [baseColors(i,:) 0.5]);
    else
        semilogx(f, rad2deg(angle(T_hann{i})), 'LineWidth', last_thickness, ...
                 'Color', [baseColors(i,:) 0.5]);
    end
end
set(gca,'XScale','log');
ylabel('\angleG (deg)');

% ---------- COHERENCE ----------
subplot(3,1,3); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        semilogx(f, C_hann{i}, 'LineWidth', line_thickness, ...
                 'Color', [baseColors(i,:) 0.5]);
    else
        semilogx(f, C_hann{i}, 'LineWidth', last_thickness, ...
                 'Color', [baseColors(i,:) 1.0]);
    end
end
set(gca,'XScale','log');
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend(legendStrings, 'Location', 'southwest');

% ---------- EXPORT ----------
width = 700;
height = 400;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB_hann.png', 'Resolution', 300);


%% Obtain Continuous Time Transfer Function
T = T_hann{length(n_win)};
data = frd(T,2*pi*f,ts); %Make FRD Data

% Use tfest() function to obtain transfer function from data
% General Configuration: sys = tfest(data,np,nz,iodelay);
% Important to model delay for controller design

np = 5; % Tune number of poles
nz = 2; % Tune number of zeros
iodelay = 2e-4; % Tune delay 
sys = tfest(data,np,nz,iodelay);
Pnump = sys.Numerator;
Pdenp = sys.Denominator;
Ptf = tf(Pnump,Pdenp);

Ptf_corr = Ptf;       
Ptf_corr.InputDelay = iodelay;

%% Plot Continuous Time Transfer Function

figure(5); clf;
h = bodeplot(data,Ptf_corr,Ptf);   % overlays original TF, delayed TF and measured FRD
grid on;
legend('FRD data input','tfest with delay','tfest without delay');
%legend('tfest model','tfest + delay','FRD');
title(sprintf('Bode Plot with Time Delay = %.1e s',iodelay));

% ---------- EXPORT ----------
width = 900;
height = 500;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB_hann_transparent.png', 'Resolution', 300);

%% COMPUTE HIGHER QUALITY TRANSFER FUNCTION

fs = 10*fmax; % sampling frequency (in Hz)
ts = 1.0/fs;
totaltime = 10; % in seconds
t = 0:ts:totaltime; % Time vector

tmax = t(end); % Time of chirp signal
method = 'logarithmic'; % or 'logarithmic'
a = chirp(t, fmin, tmax, fmax, method); % Generate a chirp signal
b = Plant(a, t); % Process the chirp signal through the system

ft = logspace(max(-1,log10(fmin)),max(3,log10(fmax)),1000);

n_win = 1;
percent_overlap = 0;

win_length = round(length(a)/n_win);
window = rectwin(length(a));             % rect window
n_overlap = floor(win_length*percent_overlap);        % 50% overlap

[T,f] = tfestimate(input,output,window,n_overlap,ft,fs);
[C,f] = mscohere(input,output,window,n_overlap,ft,fs);

data = frd(T,2*pi*f,ts); %Make FRD Data

np = 5; % Tune number of poles
nz = 2; % Tune number of zeros
iodelay = 2e-4; % Tune delay 

sys = tfest(data,np,nz,iodelay);
Pnump = sys.Numerator;
Pdenp = sys.Denominator;
Ptf = tf(Pnump,Pdenp);

Ptf_corr = Ptf;       
Ptf_corr.InputDelay = iodelay;

