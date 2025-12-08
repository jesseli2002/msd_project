%% MSD 2025 Assignment Part 1: Open-loop System Identification
% Script prepared by: Jonathan Xie

%% VARIABLES
fmin = 10;   
fmax = 1e4;
ts = 30e-6;              
fs = 1/ts;
%fs = 10*fmax; % sampling frequency (in Hz)
%ts = 1.0/fs;

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

subplot(2,1,2)
plot(t, u, 'LineWidth', 1.0);
grid on;
title('x(t), Input White Noise Signal');
ylabel('Amplitude');

subplot(2,1,2)
plot(t, y, 'LineWidth', 1.0);
grid on;
title('y(t), Output Signal from Plant.p');
ylabel('Amplitude');

% export
width = 900
height = 500
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB2_noise.png', 'Resolution', 300);

%% NEW RECT
input = u';
output = y;
% tfestimate settings
n_win = [16, 8, 4, 1];
percent_overlap = [0.5, 0.5, 0.5, 0.5];
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
line_thickness = 1.2
last_thickness = 1.2
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
legend(legendStrings, 'Location', 'best');

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

% ---------- EXPORT ----------
width = 900;
height = 700;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB_rect_transparent.png', 'Resolution', 300);


%% NEW NEW PLOTTING (RECT) — scalable version

figure(3); clf(3);
legendStrings = compose('%d windows', n_win);

subplot(3,1,1); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        semilogx(f, mag2db(abs(Tresp{i})), 'LineWidth', 1.2, 'LineStyle', '--');
    else
        semilogx(f, mag2db(abs(Tresp{i})), 'LineWidth', 1.8);  % last one solid
    end
end
set(gca,'XScale','log');
title('Frequency Response with different number of Rectangular Windows');
ylabel('|G| (dB)');
legend(legendStrings, 'Location', 'best');


subplot(3,1,2); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        semilogx(f, rad2deg(angle(Tresp{i})), 'LineWidth', 1.2, 'LineStyle', '--');
    else
        semilogx(f, rad2deg(angle(Tresp{i})), 'LineWidth', 1.8);
    end
end
set(gca,'XScale','log');
ylabel('\angleG (deg)');


subplot(3,1,3); hold on; grid on;
for i = 1:length(n_win)
    if i < length(n_win)
        semilogx(f, Cresp{i}, 'LineWidth', 1.2, 'LineStyle', '--');
    else
        semilogx(f, Cresp{i}, 'LineWidth', 1.8);
    end
end
set(gca,'XScale','log');
xlabel('Frequency (Hz)');
ylabel('Coherence');

% --- EXPORT ---
width = 900;
height = 700;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB_rect.png', 'Resolution', 300);

%% NEW PLOTTING (RECT) — scalable version

figure(3); clf(3);
legendStrings = compose('%d windows', n_win);

subplot(3,1,1); hold on; grid on;
for i = 1:length(n_win)
    semilogx(f, mag2db(abs(Tresp{i})), 'LineWidth', 1.2);
end
set(gca,'XScale','log');     % <<<<<<<<<<
title('Frequency Response with different number of Rectangular Windows');
ylabel('|G| (dB)');
legend(legendStrings, 'Location', 'best');

subplot(3,1,2); hold on; grid on;
for i = 1:length(n_win)
    semilogx(f, rad2deg(angle(Tresp{i})), 'LineWidth', 1.2);
end
set(gca,'XScale','log');     % <<<<<<<<<<
ylabel('\angleG (deg)');

subplot(3,1,3); hold on; grid on;
for i = 1:length(n_win)
    semilogx(f, Cresp{i}, 'LineWidth', 1.2);
end
set(gca,'XScale','log');     % <<<<<<<<<<
xlabel('Frequency (Hz)');
ylabel('Coherence');

% --- EXPORT ---
width = 900;
height = 700;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB_rect.png', 'Resolution', 300);

%% NEW HANNING
input = u';
output = y;
% tfestimate settings
n_win = [16, 8, 4, 1];
percent_overlap = [0.5, 0.5, 0.5, 0.5];
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

%% NEW PLOTTING (HANNING) — scalable version

figure(4); clf(4);
legendStrings = compose('%d windows', n_win);

subplot(3,1,1); hold on; grid on;
for i = 1:length(n_win)
    semilogx(f, mag2db(abs(T_hann{i})), 'LineWidth', 1.2);
end
set(gca,'XScale','log');     % <<<<<<<<<<
title('Frequency Response with different number of Rectangular Windows');
ylabel('|G| (dB)');
legend(legendStrings, 'Location', 'best');

subplot(3,1,2); hold on; grid on;
for i = 1:length(n_win)
    semilogx(f, rad2deg(angle(T_hann{i})), 'LineWidth', 1.2);
end
set(gca,'XScale','log');     % <<<<<<<<<<
ylabel('\angleG (deg)');

subplot(3,1,3); hold on; grid on;
for i = 1:length(n_win)
    semilogx(f, C_hann{i}, 'LineWidth', 1.2);
end
set(gca,'XScale','log');     % <<<<<<<<<<
xlabel('Frequency (Hz)');
ylabel('Coherence');

% --- EXPORT ---
width = 900;
height = 700;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB_hann.png', 'Resolution', 300);

%% RECT tfestimate ============================
input = u';
output = y;

%tfestimate
n_win = [8,4,1];            
percent_overlap = [0.5, 0.5, 0.5];
L = length(u);

for i = 1:length(n_win)
    win_length(i) = round( L / n_win(i) );
    n_overlap(i) = floor( win_length(i) * percent_overlap(i) );
end

[T1r,f] = tfestimate(input,output,rectwin(win_length(1)),n_overlap(1),ft,fs);
[C1r,f] = mscohere(input,output,rectwin(win_length(1)),n_overlap(1),ft,fs);

[T2r,f] = tfestimate(input,output,rectwin(win_length(2)),n_overlap(2),ft,fs);
[C2r,f] = mscohere(input,output,rectwin(win_length(2)),n_overlap(2),ft,fs);

[T3r,f] = tfestimate(input,output,rectwin(win_length(3)),n_overlap(3),ft,fs);
[C3r,f] = mscohere(input,output,rectwin(win_length(3)),n_overlap(3),ft,fs);

%% PLOTTING RECT
figure(3);clf(3);
subplot(3,1,1)
semilogx(f, mag2db(abs(T1r)), 'LineWidth', 1.2); hold on;
semilogx(f, mag2db(abs(T2r)), 'LineWidth', 1.2);
semilogx(f, mag2db(abs(T3r)), 'LineWidth', 1.2);

grid on;
title('Frequency Response with different number of Rectangular Windows');
ylabel('|G| (dB)');
legend(sprintf('%d windows',n_win(1)),sprintf('%d windows',n_win(2)),sprintf('%d windows',n_win(3)),'Location','best');

subplot(3,1,2)
semilogx(f, rad2deg(angle(T1r)), 'LineWidth', 1.2); hold on;
semilogx(f, rad2deg(angle(T2r)), 'LineWidth', 1.2);
semilogx(f, rad2deg(angle(T3r)), 'LineWidth', 1.2);
grid on;
ylabel('\angleG (deg)');

subplot(3,1,3)
semilogx(f, C1r, 'LineWidth', 1.2); hold on;
semilogx(f, C2r, 'LineWidth', 1.2);
semilogx(f, C3r, 'LineWidth', 1.2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Coherence');

% export
width = 900
height = 700
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB1_rect_unscaled.png', 'Resolution', 300);

%% HANNING tfestimate
input = u';
output = y;

%tfestimate
n_win = [16,4,8];            
percent_overlap = [0.5, 0.5, 0.5];
L = length(u);

for i = 1:length(n_win)
    win_length(i) = round( L / n_win(i) );
    n_overlap(i) = floor( win_length(i) * percent_overlap(i) );
end

[T1h,f] = tfestimate(input,output,hann(win_length(1)),n_overlap(1),ft,fs);
[C1h,f] = mscohere(input,output,hann(win_length(1)),n_overlap(1),ft,fs);

[T2h,f] = tfestimate(input,output,hann(win_length(2)),n_overlap(2),ft,fs);
[C2h,f] = mscohere(input,output,hann(win_length(2)),n_overlap(2),ft,fs);

[T3h,f] = tfestimate(input,output,hann(win_length(3)),n_overlap(3),ft,fs);
[C3h,f] = mscohere(input,output,hann(win_length(3)),n_overlap(3),ft,fs);

%% PLOTTING HANNING
figure(4);clf(4);
subplot(3,1,1)
semilogx(f, mag2db(abs(T1h)), 'LineWidth', 1.2); hold on;
semilogx(f, mag2db(abs(T2h)), 'LineWidth', 1.2);
semilogx(f, mag2db(abs(T3h)), 'LineWidth', 1.2);
grid on;
title('Frequency Response with different number of Hanning Windows');
ylabel('|G| (dB)');
legend(sprintf('%d windows',n_win(1)),sprintf('%d windows',n_win(2)),sprintf('%d windows',n_win(3)),'Location','best');

subplot(3,1,2)
semilogx(f, rad2deg(angle(T1h)), 'LineWidth', 1.2); hold on;
semilogx(f, rad2deg(angle(T2h)), 'LineWidth', 1.2);
semilogx(f, rad2deg(angle(T3h)), 'LineWidth', 1.2);
grid on;
ylabel('\angleG (deg)');

subplot(3,1,3)
semilogx(f, C1h, 'LineWidth', 1.2); hold on;
semilogx(f, C2h, 'LineWidth', 1.2);
semilogx(f, C3h, 'LineWidth', 1.2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Coherence');

% export
width = 900
height = 700
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB1_hann_unscaled.png', 'Resolution', 300);

%% Obtain Continuous Time Transfer Function

% Generate data file to be used in tfest() function
%data = iddata(output,input,ts); %Make Input-Output Data
%OR
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

figure(5); clf;
h = bodeplot(data,Ptf_corr,Ptf);   % overlays original TF, delayed TF and measured FRD
grid on;
legend('FRD data input','tfest with delay','tfest without delay');
%legend('tfest model','tfest + delay','FRD');
title(sprintf('Bode Plot with Time Delay = %.1e s',iodelay));