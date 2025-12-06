%% MSD 2025 Assignment Part 1: Open-loop System Identification
% Script prepared by: Aditya Natu, PhD Candidate, Mechatronics System Design

%% Define parameters for chirp signal
ts = 30e-6; % Sampling Time (in seconds)
fmin = 1; % Start Frequency (in Hz)
fmax = 10; % Max Frequency in Chirp (in Hz)

totaltime = 1; % in seconds
t = 0:ts:totaltime; % Time vector

% Define the chirp signal
tmax = t(end); % Time of chirp signal
% If tmax < t(end), append the signals to generate multiple chirp signals
% Total time for appended signal should be the specified totaltime (variable)

% Define chirp method
method = 'logarithmic'; % or 'logarithmic'
u = chirp(t, fmin, tmax, fmax, method); % Generate a chirp signal

%% Define the white noise (realisitic with normal distribution)
vr = 1; % Define variance
u = sqrt(vr)*randn(length(t),1)'; % Generate white noise with variance = vr

%% Call the obfuscated Plant.p function
y = Plant(u, t); % Process the chirp signal through the system

% Plot the input and output signals
figure(1);clf(1);
subplot(2, 1, 1);
plot(t, u);
title('Input Chirp Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t, y);
title('Output Signal');
xlabel('Time (s)');
ylabel('Amplitude');


%% System Identification
% This section to be run after generating the input-output signals

ts = 30e-6; %Sampling Time (seconds)
fs = 1/ts; %Sampling Frequency (Hz)

input = u';
output = y;

ft = logspace(-1,4,1000); % Define frequency vector data set
% Estimate Frequency Response using tfestimate()
L = length(u);
wind = [];% hann(), or rectwin() etc.
[T,f] = tfestimate(input,output,wind,[],ft,fs); 
[C,f] = mscohere(input,output,wind,[],ft,fs);

figure(2);clf(2);
subplot(3,1,1);semilogx(f,mag2db(abs(T))); grid on; hold on;
subplot(3,1,2);semilogx(f,rad2deg(angle(T))); grid on; hold on;
subplot(3,1,3);semilogx(f,C); grid on; hold on;

%% Obtain Continuous Time Transfer Function

% Generate data file to be used in tfest() function
data = iddata(output,input,ts); %Make Input-Output Data
%OR
data = frd(T,2*pi*f,ts); %Make FRD Data

% Use tfest() function to obtain transfer function from data
% General Configuration: sys = tfest(data,np,nz,iodelay);
% Important to model delay for controller design

np = 1; % Tune number of poles
nz = 1; % Tune number of zeros
iodelay = 1; % Tune delay 
sys = tfest(data,np,nz,iodelay);
Pnump = sys.Numerator;
Pdenp = sys.Denominator;
Ptf = tf(Pnump,Pdenp);
figure(3);clf(3);
bode(Ptf,data) %Plant Transfer Function from Identification

