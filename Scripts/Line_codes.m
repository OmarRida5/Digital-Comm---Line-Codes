A = 4;                          % Line code amplitude
samples = 7;                    % Number of samples
realization_num = 500;          % Number of realizations for each line code
ts = 0.01;                      % Sampling time of the DAC = 10 ms
tb = 0.07;                      % Bit width
tot_t = 7;                      % Total bit time
check_zero = 0;                 % A counter for the elements that will be set to zero
check_amplitude = 0;            % A counter for the elements that will not change
t = 0.01 : ts : tot_t;

%% Create a random stream of bits for 500 realizations
data = randi([0 1], 500, 101);
%% Generate a random delay for each realization
delay = randi([0 6], 500, 1); 

%% Map the binary values to amplitudes
UNIPOLAR_NRZ = data * A;
POLAR_RZ = ((2*data)-1)*A;
POLAR_NRZ = ((2*data)-1)*A;

%% Generate the samples
uni_NRZ_sampled = repelem(UNIPOLAR_NRZ, 1, samples);    
polar_NRZ_sampled = repelem(POLAR_NRZ, 1, samples);    
polar_RZ_sampled = repelem(POLAR_RZ, 1, samples);    

%% Polar return to zero procedure

Tx_out = polar_RZ_sampled;
for counter_realization = 1 : realization_num
    for counter_samples = 5 : 707
        if check_zero ~= 3 %%checking to find the elements that will be set to zero
            Tx_out(counter_realization,counter_samples) = 0;
            check_zero = check_zero + 1;
        else
            check_amplitude = check_amplitude + 1;
        end
        if check_amplitude == 4 % indicating that the next number will be set to zero
            check_amplitude = 0;
            check_zero = 0;
        end
    end
    check_amplitude = 0;
    check_zero = 0;

end

% Tx_delay=polar_RZ_sampled;
% for counter_realization=1:realization_num
%     for counter_samples=1:delay(counter_realization)
%       polar_RZ_sampled(counter_realization,counter_samples)=0;
%     end
% end
%% Delay generation
uni_NRZ_delayed = zeros(500, 700);  % Initialize matrix to store delayed uni_NRZ_sampled
polar_NRZ_delayed = zeros(500, 700);   % Initialize matrix to store delayed uni_RZ_sampled
polar_RZ_delayed = zeros(500, 700);
for i = 1 : realization_num
    % Assign delayed samples to unipolar NRZ
    uni_NRZ_delayed(i, :) = uni_NRZ_sampled(i, delay(i) + 1 : delay(i) + 700);
    
    % Assign delayed samples to polar NRZ
    polar_NRZ_delayed(i, :) = polar_NRZ_sampled(i, delay(i) + 1 : delay(i) + 700);
   
    % Assign delayed samples to polar RZ
    polar_RZ_delayed(i, :) = Tx_out(i, delay(i) + 1 : delay(i) + 700);
end
% Tx_delay=Tx_out;
% for counter_realization=1:realization
%     for counter_samples=1:delay(counter_realization)
%       Tx_delay(counter_realization,counter_samples)=0;
%     end
% end

%% Plot uni_NRZ_delayed for the first 4 realizations
figure;
for i = 1:4
    subplot(4,1,i);
    plot(t, uni_NRZ_delayed(i,:));
    title(['Delayed unipolar NRZ (Realization ', num2str(i), ')']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

%% Plot polar_NRZ_delayed for each realization
figure;
for i = 1:4
    subplot(4,1,i);
    plot(t, polar_NRZ_delayed(i,:));
    title(['Delayed polar NRZ (Realization ', num2str(i), ')']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

%% Plot delayed Polar RZ for four realizations
figure;
for i = 1:4
    subplot(4,1,i);
    plot(t, polar_RZ_delayed(i,:));
    title(['Delayed polar RZ (Realization ', num2str(i), ')']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end
%% Time auto correlation for a single waveform
%Rx(n) = sum on time{waveform_i[n1]*waveform_i[n1+m]}/700
% First we preallocate empty arrays 
uniNRZac = zeros (1, 300); 
polNRZac = zeros (1, 300); 
polRZac = zeros (1, 300); 
taw = -150: 149;
%outer loop, time shift as the variable
for m = taw
M = m + 151; %Positive array index, starts from 1
%inner loop sample position as the variable starting from 1 to 300
for sample = 151: 550 %Time samples vector, 400 time samples
uniNRZac(M) = uniNRZac(M) + uni_NRZ_delayed(1,sample)* uni_NRZ_delayed(1,sample+m);
polNRZac(M) = polNRZac(M) + polar_NRZ_delayed(1,sample)* polar_NRZ_delayed(1,sample+m);
polRZac(M) = polRZac(M) + polar_RZ_delayed(1,sample)* polar_RZ_delayed(1,sample+m);
end
uniNRZac(M) = uniNRZac(M)/400;
polNRZac(M) = polNRZac(M)/400;
polRZac(M) = polRZac(M)/400;
end

figure;
freq = -150:149;
plot(freq,polNRZac);
 title('Polar NRZ time autocorrelation');
 xlabel('f');
 ylabel('Amplitude');
 
figure;
plot(freq,uniNRZac);
 title('unipolar NRZ time autocorrelation');
 xlabel('f');
 ylabel('Amplitude');

figure;
plot(freq,polRZac);
 title('Polar RZ time autocorrelation');
 xlabel('f');
 ylabel('Amplitude');






%% Calculating Ensemble Autocorrelation
% Rx(m) = sum{waveform_i[n1]*waveform_i[n1+m]}/500

% Define the starting index for autocorrelation calculation
n1 = 351;

% Define the lag range
taw = -350:349;

% Preallocate arrays for autocorrelation results
UNI_NRZ = zeros(1, 700); % Unipolar NRZ autocorrelation preallocation
POLAR_NRZ = zeros(1, 700); % Polar NRZ autocorrelation preallocation
POLAR_RZ = zeros(1, 700); % Polar RZ autocorrelation preallocation

% Loop through each lag value
for m = taw
    % Shift the index to ensure positive array indices, starting with M=1
    M = m + n1;
    
    % Loop through each waveform
    for i = 1:500
        % Accumulate the product of delayed waveforms for each type
        % of line code (Uni_NRZ, Polar_NRZ, Polar_RZ)
        UNI_NRZ(M) = UNI_NRZ(M) + uni_NRZ_delayed(i, n1) * uni_NRZ_delayed(i, n1 + m);
        POLAR_NRZ(M) = POLAR_NRZ(M) + polar_NRZ_delayed(i, n1) * polar_NRZ_delayed(i, n1 + m);
        POLAR_RZ(M) = POLAR_RZ(M) + polar_RZ_delayed(i, n1) * polar_RZ_delayed(i, n1 + m);
    end
    
    % Average the accumulated values over all waveforms
    UNI_NRZ(M) = UNI_NRZ(M) / 500;
    POLAR_NRZ(M) = POLAR_NRZ(M) / 500;
    POLAR_RZ(M) = POLAR_RZ(M) / 500;
end

 %Plotting unipolar NRZ ensemble autocorrelation
figure;

freq = -350:349; % Define the frequency range
plot(freq, UNI_NRZ); 
title('Unipolar NRZ Ensemble Autocorrelation'); 
xlabel('Frequency (f)'); 
ylabel('Amplitude'); % 

% Plotting polar NRZ ensemble autocorrelation
figure;
plot(freq, POLAR_NRZ); % Plot autocorrelation against frequency
title('Polar NRZ Ensemble Autocorrelation'); 
xlabel('Frequency (f)'); 
ylabel('Amplitude'); 

% Plotting polar RZ ensemble autocorrelation
figure;
plot(freq, POLAR_RZ); 
title('Polar RZ Ensemble Autocorrelation'); 
xlabel('Frequency (f)'); 
ylabel('Amplitude'); 


% Calculate the statistical mean for unipolar NRZ
statistical_mean_unipolar_NRZ = zeros(1, size(uni_NRZ_delayed, 2));
for i = 1:size(uni_NRZ_delayed, 2)
    statistical_mean_unipolar_NRZ(i) = sum(uni_NRZ_delayed(:, i)) / size(uni_NRZ_delayed, 1);
end

% Calculate the statistical mean for polar NRZ
statistical_mean_polar_NRZ = zeros(1, size(polar_NRZ_delayed, 2));
for i = 1:size(polar_NRZ_delayed, 2)
    statistical_mean_polar_NRZ(i) = sum(polar_NRZ_delayed(:, i)) / size(polar_NRZ_delayed, 1);
end

% Calculate the statistical mean for polar RZ
statistical_mean_polar_RZ = zeros(1, size(polar_RZ_delayed, 2));
for i = 1:size(polar_RZ_delayed, 2)
    statistical_mean_polar_RZ(i) = sum(polar_RZ_delayed(:, i)) / size(polar_RZ_delayed, 1);
end

% Plotting the statistical mean for each signaling scheme
figure;
plot(t, statistical_mean_unipolar_NRZ, 'r', 'LineWidth', 2);
title('Statistical Mean for Unipolar NRZ');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
plot(t, statistical_mean_polar_NRZ, 'g', 'LineWidth', 2);
title('Statistical Mean for Polar NRZ');
xlabel('Time (s)');
ylabel('Amplitude');

figure;
plot(t, statistical_mean_polar_RZ, 'b', 'LineWidth', 2);
title('Statistical Mean for Polar RZ');
xlabel('Time (s)');
ylabel('Amplitude');

% Select one waveform index (for example, the first waveform)
waveform_index = 1;

% Calculate the time mean for the selected waveform for each signaling scheme
time_mean_unipolar_NRZ = calculate_time_mean(uni_NRZ_delayed(waveform_index, :));
time_mean_polar_NRZ = calculate_time_mean(polar_NRZ_delayed(waveform_index, :));
time_mean_polar_RZ = calculate_time_mean(polar_RZ_delayed(waveform_index, :));

% Plotting the time mean across realizations for each signaling scheme
figure;
plot(1:realization_num, repmat(time_mean_unipolar_NRZ, realization_num, 1), 'r', 'LineWidth', 2);
title('Time Mean Across Realizations for Unipolar NRZ');
xlabel('Realization');
ylabel('Amplitude');

figure;
plot(1:realization_num, repmat(time_mean_polar_NRZ, realization_num, 1), 'g', 'LineWidth', 2);
title('Time Mean Across Realizations for Polar NRZ');
xlabel('Realization');
ylabel('Amplitude');

figure;
plot(1:realization_num, repmat(time_mean_polar_RZ, realization_num, 1), 'b', 'LineWidth', 2);
title('Time Mean Across Realizations for Polar RZ');
xlabel('Realization');
ylabel('Amplitude');


% Define the frequency axis
Fs = 1 / ts; % Sampling frequency
N = length(UNI_NRZ); % Length of the autocorrelation function
f = linspace(-Fs/2, Fs/2, N); % Frequency axis 

% Calculate the power spectral density using the autocorrelation of the ensemble
psd_unipolar_NRZ =  abs( fftshift( fft( UNI_NRZ ) /(Fs) ) ) ;
psd_polar_NRZ = abs( fftshift(fft(POLAR_NRZ) /(Fs) ) );
psd_polar_RZ = abs( fftshift(fft(POLAR_RZ) /(Fs) ) );

figure;

% Subplot for Unipolar NRZ
subplot(3,1,1);
plot(f, psd_unipolar_NRZ, 'r', 'LineWidth', 2);
title('Power Spectral Density for Unipolar NRZ');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency');

% Subplot for Polar NRZ
subplot(3,1,2);
plot(f, psd_polar_NRZ, 'g', 'LineWidth', 2);
title('Power Spectral Density for Polar NRZ');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency');

% Subplot for Polar RZ
subplot(3,1,3);
plot(f, psd_polar_RZ, 'b', 'LineWidth', 2);
title('Power Spectral Density for Polar RZ');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency');


function mean_value = calculate_time_mean(signal)
    mean_value = sum(signal) / length(signal);
end

