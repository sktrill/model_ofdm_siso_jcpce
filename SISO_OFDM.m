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

size_tap = 3;
N = 64;
num_cyclic = 5; % Length of the prefix same as the length of channel response
symbols = sign(randn(1280,1));
num_symbols = length(symbols);
num_OFDM = num_symbols / N;

x = [];
serial_output = [];
IFFTblock_input = [];
transmitted_signal = [];
count = 0;
amp = 40;
BER = zeros(amp*2 - 1);
Error_symbols = zeros(amp*2 - 1);


for amplitude = 1:0.5:amp
    count = count + 1;
    for a=1:num_OFDM
       serial_output = [];
       n = (a-1)*N+1:(a-1)*N+N;
       IFFTblock_input = symbols(n); % IFFT block input
       x = ifft(IFFTblock_input); % IFFT block output
       for l = 1:N
           serial_output = [serial_output,x(l)]; % P/S output
       end
 % Adding the cyclic prefix
       serial_output = [serial_output(length(serial_output)-num_cyclic + 1:length(serial_output)),serial_output];
       transmitted_signal = serial_output * amplitude; 
       h = ((randn(size_tap,1) + randn(size_tap,1)*j)/sqrt(2));       
       channel_signal = conv(transmitted_signal,h); % Linear convolution with channel
       noise = (randn(length(channel_signal),1) + randn(length(channel_signal),1)*j)/sqrt(2);
       received_signal = channel_signal + noise.';
 % Get current OFDM symbol. Remove the cyclic prefix
       FFTblock_input = received_signal(num_cyclic + 1:length(received_signal)- size_tap + 1);
       FFTblock_output = fft(FFTblock_input);
       plot(FFTblock_output)
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

save ('test.mat', 'BER', 'SNR', '-ascii', '-tabs')