%% Define parameters for chirp signal

fmin = 10; % Start Frequency (in Hz)
fmax = 1e4; % Max Frequency in Chirp (in Hz)

fs_alt = 10*fmax; % sampling frequency (in Hz)
ts_alt = 1.0/fs_alt;

%ts = 30e-6; % Sampling Time (in seconds)
%fs = 1.0/ts;

totaltime = 10; % in seconds
t = 0:ts_alt:totaltime; % Time vector

% Define the chirp signal
tmax = t(end); % Time of chirp signal
% If tmax < t(end), append the signals to generate multiple chirp signals
% Total time for appended signal should be the specified totaltime (variable)

% Define chirp method
method = 'logarithmic'; % or 'logarithmic'
a = chirp(t, fmin, tmax, fmax, method); % Generate a chirp signal\

%% Call the obfuscated Plant.p function
b = Plant(a, t); % Process the chirp signal through the system

% Plot the input and output signals
figure(11);clf(11);
subplot(2, 1, 1);
plot(t, a);
% title(sprintf('Input White Noise Signal, fmin = %d , fmax = %.e, vr = %d', fmin, fmax,vr));
title(sprintf('Input Chirp Signal, fmin = %d , fmax = %.e', fmin, fmax));
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t, b);
title('Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');

input = a';
output = b;

%% Hanning Windowed System Identification
ft = logspace(max(-1,log10(fmin)),max(3,log10(fmax)),1000); % Define frequency vector data set
% Estimate Frequency Response using tfestimate()

n_win = 1;            
percent_overlap = 0.5;

win_length = round(length(a)/n_win);
window = hann(win_length);             % Hann window
n_overlap = floor(win_length*percent_overlap);        % 50% overlap

[T,f] = tfestimate(input,output,window,n_overlap,ft,fs_alt);
[C,f] = mscohere(input,output,window,n_overlap,ft,fs_alt);
  
figure(12);clf(12);
subplot(3,1,1);semilogx(f,mag2db(abs(T)));
title(sprintf('tfestimate from chirp, hanning window, num windows = %d, overlap = %0.2f', n_win, percent_overlap));
xlabel('Freq [Hz]'); ylabel('|G| [dB]'); grid on; hold on;
subplot(3,1,2);semilogx(f,rad2deg(angle(T)))
xlabel('Freq [Hz]'); ylabel('\angleG [dB]'); grid on; hold on;
subplot(3,1,3);semilogx(f,C)
xlabel('Freq [Hz]'); ylabel('Coherence'); grid on; hold on;

%% Obtain Continuous Time Transfer Function

data = frd(T,2*pi*f,ts_alt); %Make FRD Data

np = 5; % Tune number of poles
nz = 2; % Tune number of zeros
iodelay = 2e-4; % Tune delay 

sys = tfest(data,np,nz,iodelay)

%% New Post-process

% Normalize
den_const = Pdenp(end);
Pnump_norm = Pnump / den_const;
Pdenp_norm = Pdenp / den_const;
Ptf = tf(Pnump_norm, Pdenp_norm);

Ptf_corr = Ptf
Ptf_corr.InputDelay = iodelay;

%% Plot Continuous Time Transfer Function

figure(13); clf;
h = bodeplot(data,Ptf_corr,Ptf);   % overlays original TF, delayed TF and measured FRD
grid on;
legend('FRD data input','tfest with delay','tfest without delay','Location', 'best');
%legend('tfest model','tfest + delay','FRD');
title(sprintf('Bode Plot of tfest Transfer Function, Time Delay = %.1e s',iodelay));
ylim([-360*2.5,90])

% ---------- EXPORT ----------
width = 1000/1.5;
height = 600/1.5;
set(gcf, 'Position', [20, 20, width, height]);
exportgraphics(gcf, 'PartB4.png', 'Resolution', 300);

%% Post-process

Pnump = sys.Numerator;
Pdenp = sys.Denominator;
Ptf = tf(Pnump,Pdenp);

Ptf_corr = Ptf;       
Ptf_corr.InputDelay = iodelay;
