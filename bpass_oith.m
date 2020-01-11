clear all 
close all

% bandpass filter the data at a cutoff frequency and look at the
% xcorrelations and data plotted against one another

%% 
%%%%% OPEN DATA %%%%%

pt_oith = '/Users/Orenstein/Desktop/corrected_counts_091817.txt';

[raw, tt] = load_oith(pt_oith, '03-11-2015 12:00:00 PM', '08-01-2015 12:00:00 AM');

%% 
%%%%% BANDPASS FILTER %%%%%

fs = 1/24; % sampling freqency (24 cycles/day for hourly data)
f0 = 1/(24*7); % center of passband
fn = fs/2; % nyquist frequency assuming sampling of 1 hour
f0n = f0/fn; % normalize the frequency for the notch
n_width = 0.001; % passband width
con = 1.25; %stopband multiplier

bdesign = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    f0n-(con*n_width), f0n-n_width, f0n+n_width, f0n+(con*n_width), 60, 1, 60,fs);

bpass = design(bdesign);

freqz(bpass)

%test = filter(bpass, (raw.oith-mean(raw.oith))./std(raw.oith));
test = abs(hilbert(filter(bpass, (raw.oith-mean(raw.oith)))));
figure; plot(test)
hold on
plot(abs(hilbert((raw.oith-mean(raw.oith))))) %./std(raw.oith))

%% 
%%%%% BUTTERWORTH IIR %%%%%
Fs = 1/24; % Sampling frequency
Fn = Fs/2; % Nyquist
Wp = [0.005 0.007]/Fn; % passband
Ws = [0.0045 0.0075]/Fn; % stopband
Rp = 1; % Ppassband ripple
Rs = 25; % stopband ripple
[n, Wn] = buttord(Wp, Ws, Rp, Rs);
[b, a] = butter(n, Wn);
[sos, g] = tf2sos(b, a);

test = filtfilt(sos, g, (raw.oith - mean(raw.oith)));

freqz(sos)

figure; plot(test)
hold on 
plot(raw.oith - mean(raw.oith))
hold off

[pxx, freq] = pwelch(test, 1024, [], [], Fs);
figure; semilogx(freq, 10*log10(pxx))

[pxx, freq] = pwelch(raw.oith - mean(raw.oith), 1024, [], [], Fs);
figure; semilogx(freq, 10*log10(pxx))

%%

fs = 1/24; % sampling freqency (24 cycles/day for hourly data)
f0 = (1/(24*7)*(1/24)); % center of passband
fn = fs/2; % nyquist frequency assuming sampling of 1 hour
f0n = f0/fn; % normalize the frequency for the notch
n_width = 0.001; % passband width
con = 1.25; %stopband multiplier

bdesign = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    f0n-(con*n_width), f0n-n_width, f0n+n_width, f0n+(con*n_width), 60, 1, 60, fs);

bpass = design(bdesign, 'equiripple');

freqz(bpass, 512, fs)

test = abs(hilbert(filter(bpass, (raw.oith - mean(raw.oith)))));
figure; plot(test)
hold on
plot(abs(hilbert((raw.oith - mean(raw.oith)))))


[pxx, freq] = pwelch(test, 1024, [], [], fs);
figure; semilogx(freq, 10*log10(pxx))

[pxx, freq] = pwelch(raw.oith - mean(raw.oith), 1024, [], [], fs);
figure; semilogx(freq, 10*log10(pxx))

%%

filt.oith = abs(hilbert(filter(bpass, (raw.oith - mean(raw.oith)))));
filt.para = abs(hilbert(filter(bpass, (raw.para - mean(raw.para)))));
filt.egg = abs(hilbert(filter(bpass, (raw.egg - mean(raw.egg)))));

[C, lags] = xcorr(filt.oith, filt.para, 'coeff');
figure; plot(lags, C)