% =========================================================
% SISO OFDM SYSTEM MODEL using BPSK modulation
%
% 'N' is the number of subcarriers
% 'symbols' is a (-1,1) sequence of desired length
% 'num_symbols' is the number of symbols to be transmitted
% 'num_OFDM' is the number of OFDM symbols to produce
% 'num_cyclic' is the length of the cyclic prefix
% 'h' is the channel impulse response
% 'H' is the channel frequency response
% 'OFDM_symbol' is the OFDM symbol received
%==========================================================

SAMPLING_RATE = 20000000; % Symbol Sampling Rate = 20 MHz
sampling_period = 1/SAMPLING_RATE; % Sampling Period = 50ns
size_tap = 3;
N = 64;
num_cyclic = 5; % Length of the prefix same as the length of channel response
symbols = sign(randn(640,1));
num_symbols = length(symbols);
num_OFDM = num_symbols / N;
freq_offset = SAMPLING_RATE/N * 0.10; %Frequency Offset 5% of carrier spacing
alpha = 0.04;
phase_noise = zeros(1,N);
amp = 10;

x = [];
serial_output = [];
IFFTblock_input = [];
transmitted_signal = [];
count = 0;
BER = zeros(amp*2 - 1);
Error_symbols = zeros(amp*2 - 1);

s = N;
l = size_tap;
est = 100;
eps = freq_offset/SAMPLING_RATE; % normalized freq. offset
CFO_est = zeros(est,s);
L_eps = zeros(est,1);
phi = zeros(s,s);


for amplitude = 1:0.5:amp
    count = count + 1;
    for a=1:num_OFDM
       y = count;
       x = num_OFDM;
       
       serial_output = [];
       n = (a-1)*N+1:(a-1)*N+N;
       IFFTblock_input = symbols(n); % IFFT block input
       x = ifft(IFFTblock_input); % IFFT block output
       for m = 1:N
           serial_output = [serial_output,x(m)]; % P/S output
       end
 % Adding the cyclic prefix
       serial_output = [serial_output(length(serial_output)-num_cyclic + 1:length(serial_output)),serial_output];
       transmitted_signal = serial_output * amplitude; 
       h = ((randn(size_tap,1) + randn(size_tap,1)*j)/sqrt(2));       
       
     %************** MODELING THE WEINER PHASE NOISE **************************
     AWGN_noise = randn(1,N)*alpha;
     phase_noise(1) = AWGN_noise(1);
     for n = 2:N
         phase_noise(1,n) = phase_noise(1,n-1) + AWGN_noise(n);
     end
     %*************************************************************************
       
       channel_signal = conv(transmitted_signal,h); % Linear convolution with channel
       noise = (randn(length(channel_signal),1) + randn(length(channel_signal),1)*j)/sqrt(2);
       received_signal = channel_signal + noise.';
 % Get current OFDM symbol. Remove the cyclic prefix
       FFTblock_input = received_signal(num_cyclic + 1:length(received_signal)- size_tap + 1);
       FFTblock_input_freqoffset = FFTblock_input.*exp(j*2*pi*freq_offset*(1:N)*sampling_period);
       FFTblock_input_phaseoffset = FFTblock_input_freqoffset.*exp(j*phase_noise(1,1:N));
       FFTblock_output = fft(FFTblock_input_phaseoffset);
       H = fft(h,length(FFTblock_input)); % Channel Frequency response
       OFDM_symbol = FFTblock_output.*(H)';
       for k = 1:N
           received_symbol = sign(real(OFDM_symbol(k)));
           if(received_symbol~=symbols((a-1)*N+k))
               Error_symbols(count) = Error_symbols(count) + 1;
           end
       end
    end
   BER(count) = Error_symbols(count)/num_symbols;
end

amplitude = 1:0.5:amp;
SNR = 20*log10(amplitude);
semilogy(SNR, BER);
title('BER vs SNR for a SISO OFDM System');
xlabel('SNR (db)');
ylabel('BER (log scale)');

save ('test2.mat', 'BER', 'SNR', '-ascii', '-tabs')