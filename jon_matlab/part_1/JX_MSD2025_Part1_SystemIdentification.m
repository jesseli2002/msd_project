%% MSD 2025 Assignment Part 1: Open-loop System Identification
% Script prepared by: Aditya Natu, PhD Candidate, Mechatronics System Design

%% Define parameters for chirp signal

fmin = 10; % Start Frequency (in Hz)
fmax = 1e4; % Max Frequency in Chirp (in Hz)

fs = 10*fmax; % sampling frequency (in Hz)
ts = 1.0/fs;

%ts = 30e-6; % Sampling Time (in seconds)
%fs = 1.0/ts;

totaltime = 10; % in seconds
t = 0:ts:totaltime; % Time vector

% Define the chirp signal
tmax = t(end); % Time of chirp signal
% If tmax < t(end), append the signals to generate multiple chirp signals
% Total time for appended signal should be the specified totaltime (variable)

% Define chirp method
method = 'logarithmic'; % or 'logarithmic'
u = chirp(t, fmin, tmax, fmax, method); % Generate a chirp signal

%% Multi Chirp
fmin = 10;   
fmax = 1e4;
ts = 30e-6;              
fs = 1/ts;

totaltime = 10;      % total time duration
t = 0:ts:totaltime;

n_chirp = 1;                                 % number of chirps
chirp_time = totaltime/n_chirp;              % each chirp duration
t_single = 0:ts:chirp_time;                  % time for one chirp

method = 'logarithmic';
u_single = chirp(t_single, fmin, t_single(end), fmax, method);

u = repmat(u_single, 1, n_chirp);            % repeat chirps
u = u(1:length(t));                          % ensure same length

%% Define the white noise (realisitic with normal distribution)
vr = 10; % Define variance
u = sqrt(vr)*randn(length(t),1)'; % Generate white noise with variance = vr
%% HIGH FREQ
% instantaneous frequency of chirp
f_inst = fmin * (fmax/fmin).^(t/tmax);  % logarithmic chirp
% Smooth amplitude ramp
f_boost_start = 2e3;   % 1 kHz
f_boost_end   = 1e4;   % 10 kHz
amp_min = 1;           % low-frequency amplitude
amp_max = 1e2;           % high-frequency amplitude
% hard ramp
%amp = ones(size(t));          % default 1
%amp(f_inst > f_boost_start) = amp_max;             % boost 10x for f > 10 kHz
% Linear ramp
%amp = amp_min + (amp_max-amp_min) .* ...
%      min(max((f_inst-f_boost_start)/(f_boost_end-f_boost_start),0),1);
% Apply envelope

% Logarithmic ramp
amp = amp_min + (amp_max-amp_min) * log10(f_inst/f_boost_start)/log10(f_boost_end/f_boost_start);
amp(f_inst<f_boost_start) = amp_min;
amp(f_inst>f_boost_end) = amp_max;

u = u .* amp;

%% Call the obfuscated Plant.p function
y = Plant(u, t); % Process the chirp signal through the system

% Plot the input and output signals
figure(1);clf(1);
subplot(2, 1, 1);
plot(t, u);
% title(sprintf('Input White Noise Signal, fmin = %d , fmax = %.e, vr = %d', fmin, fmax,vr));
title(sprintf('Input Chirp Signal, fmin = %d , fmax = %.e', fmin, fmax));
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t, y);
title('Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');

input = u';
output = y;

%% Hanning Windowed System Identification
ft = logspace(max(-1,log10(fmin)),max(3,log10(fmax)),1000); % Define frequency vector data set
% Estimate Frequency Response using tfestimate()

n_win = 1;            
percent_overlap = 0.5;

win_length = round(length(u)/n_win);
window = hann(win_length);             % Hann window
n_overlap = floor(win_length*percent_overlap);        % 50% overlap

[T,f] = tfestimate(input,output,window,n_overlap,ft,fs);
[C,f] = mscohere(input,output,window,n_overlap,ft,fs);
  
figure(2);clf(2);
subplot(3,1,1);semilogx(f,mag2db(abs(T)));
title(sprintf('tfestimate from chirp, hanning window, num windows = %d, overlap = %0.2f', n_win, percent_overlap));
xlabel('Freq [Hz]'); ylabel('|G| [dB]'); grid on; hold on;
subplot(3,1,2);semilogx(f,rad2deg(angle(T)))
xlabel('Freq [Hz]'); ylabel('\angleG [dB]'); grid on; hold on;
subplot(3,1,3);semilogx(f,C)
xlabel('Freq [Hz]'); ylabel('Coherence'); grid on; hold on;
%ylim([0.95 1.05]);

%% Regular System Identification
ft = logspace(max(-1,log10(fmin)),max(3,log10(fmax)),1000); % Define frequency vector data set
% Estimate Frequency Response using tfestimate()
L = length(u);
wind = [];% hann(), or rectwin() etc.
[T,f] = tfestimate(input,output,wind,[],ft,fs); 
[C,f] = mscohere(input,output,wind,[],ft,fs);
  
figure(2);clf(2);
subplot(3,1,1);semilogx(f,mag2db(abs(T)));
title('tfestimate');
xlabel('Freq [Hz]'); ylabel('|G| [dB]'); grid on; hold on;
subplot(3,1,2);semilogx(f,rad2deg(angle(T)))
xlabel('Freq [Hz]'); ylabel('\angleG [dB]'); grid on; hold on;
subplot(3,1,3);semilogx(f,C)
xlabel('Freq [Hz]'); ylabel('Coherence'); grid on; hold on;

%% Rectangular Windowed System Identification
ft = logspace(max(-1,log10(fmin)),max(3,log10(fmax)),1000); % Define frequency vector data set
% Estimate Frequency Response using tfestimate()

n_win = 32;
percent_overlap = 0.5;

win_length = round(length(u)/n_win);
window = hann(win_length);             % Hann window
n_overlap = floor(win_length*percent_overlap);        % 50% overlap

[T,f] = tfestimate(input,output,window,n_overlap,ft,fs);
[C,f] = mscohere(input,output,window,n_overlap,ft,fs);
  
figure(2);clf(2);
subplot(3,1,1);semilogx(f,mag2db(abs(T)));
title(sprintf('tfestimate, rectangular window, num windows = %d, overlap = %0.2f', n_win, percent_overlap));
xlabel('Freq [Hz]'); ylabel('|G| [dB]'); grid on; hold on;
subplot(3,1,2);semilogx(f,rad2deg(angle(T)))
xlabel('Freq [Hz]'); ylabel('\angleG [dB]'); grid on; hold on;
subplot(3,1,3);semilogx(f,C)
xlabel('Freq [Hz]'); ylabel('Coherence'); grid on; hold on;


%% OLD Windowed System Identification
ft = logspace(max(-1,log10(fmin)),max(3,log10(fmax)),1000); % Define frequency vector data set
% Estimate Frequency Response using tfestimate()
L = length(u);

n_win = 2048*8;                     % Window length (samples)
window = hann(n_win);             % Hann window
n_overlap = floor(n_win/2);        % 50% overlap
Nfft = n_win;                     % FFT size

[T,f] = tfestimate(input,output,window,n_overlap,Nfft,fs);
[C,f] = mscohere(input,output,window,n_overlap,Nfft,fs);
  
figure(2);clf(2);
subplot(3,1,1);semilogx(f,mag2db(abs(T)));
title(sprintf('tfestimate, hanning window, N win = %d, N overlap = %d', n_win, n_overlap));
xlabel('Freq [Hz]'); ylabel('|G| [dB]'); grid on; hold on;
subplot(3,1,2);semilogx(f,rad2deg(angle(T)))
xlabel('Freq [Hz]'); ylabel('\angleG [dB]'); grid on; hold on;
subplot(3,1,3);semilogx(f,C)
xlabel('Freq [Hz]'); ylabel('Coherence'); grid on; hold on;

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

Ptf_del = Ptf;       
Ptf_del.InputDelay = iodelay;

figure(3); clf;
h = bodeplot(data,Ptf_del,Ptf);   % overlays original TF, delayed TF and measured FRD
grid on;
legend('FRD data input','tfest with delay','tfest without delay');
%legend('tfest model','tfest + delay','FRD');
title(sprintf('Bode Plot with Time Delay = %.1e s',iodelay));

%% OLD Continuous Time Transfer Function

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
figure(3);clf(3);
bode(Ptf,data) %Plant Transfer Function from Identification
title('Bode Plot using tfest, from Chirp')
legend('tfest Model','FRD from Chirp');

%% thanks GPT
% 1. Create iddata from time-domain signals
dataID = iddata(output, input, ts);

% 2. Estimate transfer function
np = 5;
nz = 2;
iodelay = 1;
sys = tfest(dataID, np, nz, iodelay);

% 3. Compute frequency-domain estimate (optional)
[T,f] = tfestimate(input, output, [], [], ft, fs);
dataFRD = frd(T, 2*pi*f, ts);

% 4. Plot Bode
figure;
bode(sys, dataFRD);  % overlay fitted model and measured data
grid on;
