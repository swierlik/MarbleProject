% wavelet_denoise_optimization.m

% Load the original deterministic and noisy data
load('Data/lorenzData.mat', 'sol', 't', 'dt');
load('Data/systemData.mat', 'V_x', 'A_x', 'B_x', 'xReg', 'r', 'tspan', 'dt', 'S_x');

load('Data/lorenzDataStochastic.mat', 'sol2', 't', 'dt');
load('Data/systemDataStochastic.mat', 'V_x2', 'A_x2', 'B_x2', 'xReg2', 'r2', 'tspan2', 'dt2', 'S_x2');

% Extract x from original and noisy data
x_original = sol(:, 1);
x_noisy = sol2(:, 1);

% Define L (exclude initial and final transients)
L = 300:length(xReg)-300;

% Set reproducibility seed
rng(117);

% List of wavelets to test
wavelet_families = {'haar', 'db2', 'db4', 'sym2', 'sym4', 'coif2', 'bior1.3', 'bior3.5', 'rbio1.3', 'rbio3.5', 'dmey'};
levels = 1:9; % Test wavelet denoising at multiple levels

% Initialize best parameters
best_wavelet = '';
best_level = 1;
min_mse = Inf;
mse_results = zeros(length(wavelet_families), length(levels));

% Loop over all wavelet families and levels
for w = 1:length(wavelet_families)
    for l = 1:length(levels)
        wavelet = wavelet_families{w};
        level = levels(l);
        
        % Apply wavelet denoising
        x_denoised = wdenoise(x_noisy, level, 'Wavelet', wavelet);
        
        % Compute MSE between original and denoised x
        mse = mean((x_original(L) - x_denoised(L)).^2);
        mse_results(w, l) = mse;
        
        % Track best parameters
        if mse < min_mse
            min_mse = mse;
            best_wavelet = wavelet;
            best_level = level;
        end
    end
end

% Display best wavelet and level
disp(['Best Wavelet: ', best_wavelet]);
disp(['Best Level: ', num2str(best_level)]);
disp(['Minimum MSE: ', num2str(min_mse)]);

% Plot comparison of original, noisy, and best denoised signal
x_best_denoised = wdenoise(x_noisy, best_level, 'Wavelet', best_wavelet);

figure;
plot(tspan(L), x_original(L), 'b', 'LineWidth', 1.5);
hold on;
plot(tspan(L), x_noisy(L), 'r', 'LineWidth', 1.2);
plot(tspan(L), x_best_denoised(L), 'g', 'LineWidth', 1.2);
title('Comparison of Original, Noisy, and Optimized Wavelet Denoised Signals');
xlabel('Time');
ylabel('x');
legend('Original', 'Noisy', 'Best Wavelet Denoised');
grid on;

% Plot MSE heatmap
figure;
imagesc(levels, 1:length(wavelet_families), mse_results);
colorbar;
yticks(1:length(wavelet_families));
yticklabels(wavelet_families);
xlabel('Wavelet Level');
ylabel('Wavelet Family');
title('MSE Comparison of Different Wavelet Denoising Techniques');

grid on;
