
%% load data
load('Collocated_spyder.mat');

%% scaling
% ----- Parameters -----
amp_min = 1;      % amplitude at start of ramp
amp_max = 1e3;    % amplitude at end of ramp
ramp_start_frac = 0.3;  % fraction of vector where ramp starts (30%)
% ----------------------

N = length(in);

% Define ramp start and end indices
start_idx = floor(ramp_start_frac * N);   % start of ramp
end_idx   = N;                             % end of ramp

% Initialize amplitude vector
amp = ones(N,1) * amp_min;  % default amplitude = amp_min

% Logarithmic ramp
ramp_len = end_idx - start_idx + 1;
ramp_idx = 0:(ramp_len-1);  % indices for ramp
amp_ramp = 10.^( log10(amp_min) + (log10(amp_max)-log10(amp_min)) * (ramp_idx/(ramp_len-1)) );

% Apply ramp
amp(start_idx:end_idx) = amp_ramp;

% Scale input
in = in .* amp;
%out = out .* amp;
figure(1)
plot(in)

%% no window
Ts = 1e-4;
Fs = 1/Ts;
f=logspace(0,3,1001);
[H,f] = tfestimate(in,out,[],[],f,Fs);
[C,f] = mscohere(in,out,[],[],f,Fs);

figure(2);
subplot(3,1,1);semilogx(f,20*log10(abs(H)));grid on;
ylabel('G(dB)');
subplot(3,1,2);semilogx(f,rad2deg(angle(H)));grid on;
ylabel('\angleG(deg)');
subplot(3,1,3);semilogx(f,C);grid on;
xlabel('Frequency(Hz)');
ylabel('Coherence');

%% hanning
load('Collocated_spyder.mat');
Ts = 1e-4;
Fs = 1/Ts;
f=logspace(0,3,1001);
LW=round(length(in)/3);
Window=rectwin(LW);
% Window=hann(LW); % uncomment this line if you want to check hanning window
n=100; % choose any number you want to check the effect of overlap, 
% you should not see much difference in this case since the input is chirp, this will be relevant if the multisine input is used
[H,f] = tfestimate(in,out,Window,n,f,Fs);
[C,f] = mscohere(in,out, Window,n,f,Fs);
figure(1);
subplot(3,1,1);semilogx(f,20*log10(abs(H)));grid on;
ylabel('G(dB)');
subplot(3,1,2);semilogx(f,rad2deg(angle(H)));grid on;
ylabel('\angleG(deg)');
subplot(3,1,3);semilogx(f,C);grid on;
xlabel('Frequency(Hz)');
ylabel('Coherence')