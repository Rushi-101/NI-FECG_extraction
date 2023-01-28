%% Variables
[fecg_signal, fs, tm] = rdsamp('sub01_snr12dB_l1_fecg1.dat');
[mecg_signal, fs, tm] = rdsamp('sub01_snr12dB_l1_mecg.dat');
channels = [6,10,20,30,34];

for i=1:length(channels)   
    signal = fecg_signal(:,1:channels(i)) + mecg_signal(:,1:channels(i));
    %new_signal = awgn(signal,10,'measured');
    A = signal;
end
%datar01 = readtable('data_r01.csv');
%A = [datar01.Var2,datar01.Var3,datar01.Var4,datar01.Var5];
Fs = 250;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 3e5;             % Length of signal
t = (0:L-1)*T;        % Time vector

%% Time series plot and notch filtering
tsin = timeseries(A);
interval = [0.0498,0.0502];
tsoutnotch = idealfilter(tsin,interval,'notch');
tsoutnotchmean = tsoutnotch + 500;
figure;
plot(tsin,'-.');
hold on;
plot(tsoutnotchmean);
xlim([0 5000])
xlabel('Time'), ylabel('ECG signal');
title('ECG signal v/s Time (Datar01.Var1)');
legend('Original Data','Filtered Data',...
       'Location','NorthWest');
hold off;

%% FFT
f = Fs*(0:(L/2))/L;
Y = fft(A);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("Frequency (Hz)")
ylabel("Magnitude")

X = fft(tsoutnotch.Data);
P4 = abs(X/L);
P3 = P4(1:L/2+1);
P3(2:end-1) = 2*P3(2:end-1);

figure;
plot(f,P3)
title("Single-Sided Amplitude Spectrum of S(t) after notch filtering at 50Hz")
xlabel("Frequency (Hz)")
ylabel("Magnitude")

%% LPF
blo = fir1(1000,0.3);
lps = filter(blo,1,A);
lps_0phase = filtfilt(blo,1,A);

figure;
freqz(blo,1);

figure;
subplot(2,1,1);
plot(lps);
xlim([0 L/60])
title("LPF signal filtered at 150Hz")
xlabel("Time")
ylabel("ECG signal")

subplot(2,1,2);
plot(lps_0phase);
xlim([0 L/60])
title("LPF signal filtered at 150Hz (Zero phase filtering)")
xlabel("Time")
ylabel("ECG signal")

X = fft(lps);
P4 = abs(X/L);
P3 = P4(1:L/2+1);
P3(2:end-1) = 2*P3(2:end-1);

figure;
plot(f,P3) 
title("Single-Sided Amplitude Spectrum of S(t) after LPF at 150Hz")
xlabel("Frequency (Hz)")
ylabel("Magnitude")

X = fft(lps_0phase);
P4 = abs(X/L);
P3 = P4(1:L/2+1);
P3(2:end-1) = 2*P3(2:end-1);

figure;
plot(f,P3) 
title("Single-Sided Amplitude Spectrum of S(t) after LPF at 150Hz (Zero phase filtering)")
xlabel("Frequency (Hz)")
ylabel("Magnitude")

figure;
pspectrum([lps,A],Fs)

% HPF

% hpfreq = [0.1 0.2 0.3 0.5 0.75 1 2 3 4 5];
% for i=1:length(hpfreq)
%     [hps,d2] = highpass(lps,hpfreq(i),Fs,ImpulseResponse="fir",Steepness=0.7);
%     figure;
%     plot(hps);
%     xlim([0 L/60])
%     title("HPF signal filtered at " + hpfreq(i) + "Hz")
%     xlabel("Time")
%     ylabel("ECG signal")
% end
bhi = fir1(1000,0.0015,'high');
hps = filter(bhi,1,lps_0phase);
hps_0phase = filtfilt(bhi,1,lps_0phase);

figure;
freqz(bhi,1);

figure;
subplot(2,1,1);
plot(hps);
xlim([0 L/60])
title("HPF signal filtered at 0.75Hz")
xlabel("Time")
ylabel("ECG signal")

subplot(2,1,2);
plot(hps_0phase);
xlim([0 L/60])
title("HPF signal filtered at 0.75Hz (Zero phase filtering)")
xlabel("Time")
ylabel("ECG signal")

Principal Component Analysis
[out_comps,qrsmethod,W,norm_data] = FECGSYN_bss_extraction(A,'PCA',Fs,300000,1);

figure;
for i=1:4
    subplot(4,1,i)
    plot(A(1:end,i),'-');
    xlim([0 Fs*5])
    xlabel('Time'), ylabel('ECG signal');
    title("ECG signal v/s Time (Var"+i+")");
end
out_comps = out_comps';

figure;
for i=1:4
    subplot(4,1,i)
    plot(out_comps(1:end,i),'-');
    xlim([0 Fs*5])
    xlabel('Time'), ylabel('ECG signal');
    title("Prinicipal component"+i);
end

figure;
for i=1:4
    subplot(4,1,i)
    plot(out_comps(1:end,i+4),'-');
    xlim([0 Fs*5])
    xlabel('Time'), ylabel('ECG signal');
    title("Prinicipal component"+i);
end


figure;
for i=1:4
    fecg = norm_data(1:end,i) - out_comps(1:end,2);
    subplot(4,1,i)
    plot(fecg,'-');
    xlim([0 Fs*5])
    xlabel('Time'), ylabel('ECG signal');
    title("fECG"+i);
end
figure;
for i=1:4
    subplot(4,1,i)
    plot(fecg_signal(1:end,i),'-');
    xlim([0 Fs*5])
    xlabel('Time'), ylabel('FECG signal');
    title("fECG"+i);
end
figure;
for i=1:4
    subplot(4,1,i)
    plot(mecg_signal(1:end,i),'-');
    xlim([0 Fs*5])
    xlabel('Time'), ylabel('MECG signal');
    title("mECG"+i);
end

[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(out_comps(1:end,4),fs,1);
[qrs_amp_raw1,qrs_i_raw1,delay1] = pan_tompkin(fecg_signal(1:end,1),fs,1);

qrs_i_raw = qrs_i_raw';
qrs_i_raw1 = qrs_i_raw1';
x = 0;
for i=1:762
    for j=1:695
        if (qrs_i_raw(i)>=qrs_i_raw1(j)-10)&&(qrs_i_raw(i)<=qrs_i_raw1(j)+10)
            x = x+1;
        end
    end
end
