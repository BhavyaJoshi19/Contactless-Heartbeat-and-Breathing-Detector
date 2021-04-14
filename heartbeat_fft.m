close all
X = val';
Fs = 3600/10; %  [samples/s] sampling frequency 
T  = 1/Fs;      % [s] sampling period       
N  = 3600;       % [samples] Length of signal
t  = (0:N-1)*T; % [s] Time vector
deltaF = Fs/N; % [1/s]) frequency intervalue of discrete signal
%figure('color','white','position',[70 100 600 900]); 
%subplot(3,1,1);
figure
plot(1e3*t,X)
title({'Heartbeat data'})
xlabel('t (milliseconds)')
ylabel('X(t)')
% compute the fast fourier transform
Y1 = fft(X);
% manually shifting the FFT
Y2 = abs(Y1/N);
Y2(1)=1;
Amp = [Y2(ceil(end/2)+1:end)' Y2(end) Y2(2:ceil(end/2))']';
if (mod(N,2) == 0)
sampleIndex = -N/2:1:N/2-1; %raw index for FFT plot
else
sampleIndex = -(N-1)/2:1:(N-1)/2; %raw index for FFT plot
end
%subplot(3,1,2);
figure
plot(deltaF*sampleIndex, Amp);
hold on;
idx = find(Amp > 3.5);
idx(sampleIndex(idx) < 0) = [];
plot(deltaF*sampleIndex(idx), Amp(idx), '+');
for k = 1:length(idx)
    if (idx(k) > (N-1)/2 && Amp(idx(k))>3.5)
        h = text(deltaF*sampleIndex(idx(k)), Amp(idx(k))+0.15,...
            ['f=' num2str(deltaF*sampleIndex(idx(k))) ' Hz']);
        set(h,'rotation',90)
    end
end
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title(['Heartbeat = ' num2str(deltaF*sampleIndex(idx(1))) ' Hz = ' ...
    num2str(60.0*(deltaF*sampleIndex(idx(1)))) ' BPM']);
xlim([0 max(deltaF*sampleIndex)/4]);
%subplot(3,1,3);
figure
half_f = deltaF*(0:(N/2));
plot(fftshift([half_f -fliplr(half_f(2:end+mod(N,2)-1))]), ...
    abs(fftshift(Y2)/N));
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Using fftshift');