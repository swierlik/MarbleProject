% noise_optimization.m (Updated for SVWE)

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

% Compute singular value weighted error for noisy system
sys_x2 = ss(A_x2, B_x2, eye(r-1), zeros(r-1, 1));
[y_sim_x_noisy, ~] = lsim(sys_x2, xReg2(L, r2), dt2 * (L - 1), xReg2(1, 1:r2-1));
error_noisy = sqrt(sum(S_x2(1:r-1) .* (V_x2(L,1:r-1) - y_sim_x_noisy(:,1:r-1)).^2, 2));
mean_error_noisy = mean(error_noisy);

% ------------------ Optimization for Smoothing (Using SVWE) ------------------

% Moving Average Filter Optimization
best_ma_window_size = 5;
min_error_ma = Inf;
for window_size = 3:2:21
    x_smooth_ma = movmean(x_noisy, window_size);
    [V_ma, A_x_ma, B_x_ma, xReg_ma, S_ma] = getSystem(x_smooth_ma, 100, r, dt, tspan);
    sys_x_ma = ss(A_x_ma, B_x_ma, eye(r-1), zeros(r-1, 1));
    [y_sim_x_ma, ~] = lsim(sys_x_ma, xReg_ma(L, r), dt * (L - 1), xReg_ma(1, 1:r-1));
    error_ma = sqrt(sum(S_ma(1:r-1) .* (V_ma(L,1:r-1) - y_sim_x_ma(:,1:r-1)).^2, 2));
    mean_error_ma = mean(error_ma);
    if mean_error_ma < min_error_ma
        min_error_ma = mean_error_ma;
        best_ma_window_size = window_size;
    end
end

% Savitzky-Golay Filter Optimization
best_sg_order = 3;
best_sg_framelen = 7;
min_error_sg = Inf;
for order = 2:4
    for framelen = max(5, order+2):2:50
        if mod(framelen, 2) == 1
            x_smooth_sg = sgolayfilt(x_noisy, order, framelen);
            [V_sg, A_x_sg, B_x_sg, xReg_sg, S_sg] = getSystem(x_smooth_sg, 100, r, dt, tspan);
            sys_x_sg = ss(A_x_sg, B_x_sg, eye(r-1), zeros(r-1, 1));
            [y_sim_x_sg, ~] = lsim(sys_x_sg, xReg_sg(L, r), dt * (L - 1), xReg_sg(1, 1:r-1));
            error_sg = sqrt(sum(S_sg(1:r-1) .* (V_sg(L,1:r-1) - y_sim_x_sg(:,1:r-1)).^2, 2));
            mean_error_sg = mean(error_sg);
            if mean_error_sg < min_error_sg
                min_error_sg = mean_error_sg;
                best_sg_order = order;
                best_sg_framelen = framelen;
            end
        end
    end
end

% Wavelet Denoising Optimization
best_wd_level = 1;
min_error_wd = Inf;
for level = 1:9
    x_smooth_wd = wdenoise(x_noisy, level);
    [V_wd, A_x_wd, B_x_wd, xReg_wd, S_wd] = getSystem(x_smooth_wd, 100, r, dt, tspan);
    sys_x_wd = ss(A_x_wd, B_x_wd, eye(r-1), zeros(r-1, 1));
    [y_sim_x_wd, ~] = lsim(sys_x_wd, xReg_wd(L, r), dt * (L - 1), xReg_wd(1, 1:r-1));
    error_wd = sqrt(sum(S_wd(1:r-1) .* (V_wd(L,1:r-1) - y_sim_x_wd(:,1:r-1)).^2, 2));
    mean_error_wd = mean(error_wd);
    if mean_error_wd < min_error_wd
        min_error_wd = mean_error_wd;
        best_wd_level = level;
    end
end

% ------------------ Final Error Comparison ------------------

% Compute best smoothing results
x_best_ma = movmean(x_noisy, best_ma_window_size);
x_best_sg = sgolayfilt(x_noisy, best_sg_order, best_sg_framelen);
x_best_wd = wdenoise(x_noisy, best_wd_level);

[V_best_ma, ~, ~, ~, S_best_ma] = getSystem(x_best_ma, 100, r, dt, tspan);
[V_best_sg, ~, ~, ~, S_best_sg] = getSystem(x_best_sg, 100, r, dt, tspan);
[V_best_wd, ~, ~, ~, S_best_wd] = getSystem(x_best_wd, 100, r, dt, tspan);

error_ma_best = sqrt(sum(S_best_ma(1:r-1) .* (V_x(L,1:r-1) - V_best_ma(L,1:r-1)).^2, 2));
error_sg_best = sqrt(sum(S_best_sg(1:r-1) .* (V_x(L,1:r-1) - V_best_sg(L,1:r-1)).^2, 2));
error_wd_best = sqrt(sum(S_best_wd(1:r-1) .* (V_x(L,1:r-1) - V_best_wd(L,1:r-1)).^2, 2));

% ------------------ Plot Results ------------------

figure;
bar([mean_error_noisy, min_error_ma, min_error_sg, min_error_wd]);
set(gca, 'XTickLabel', {'Noisy', 'MA Optimized', 'SG Optimized', 'WD Optimized'});
title('Comparison of Minimum Errors After Optimization');
xlabel('Method');
ylabel('Mean SVWE Error');
grid on;

figure;
plot(tspan(L), error_noisy, 'r', 'LineWidth', 1);
hold on;
plot(tspan(L), error_ma_best, 'g', 'LineWidth', 1);
plot(tspan(L), error_sg_best, 'm', 'LineWidth', 1);
plot(tspan(L), error_wd_best, 'k', 'LineWidth', 1);
legend('Noisy', 'MA Optimized', 'SG Optimized', 'WD Optimized');
title('Error Evolution Over Time (SVWE)');
xlabel('Time');
ylabel('Error');
grid on;

disp(['Best Moving Average Window: ', num2str(best_ma_window_size)]);
disp(['Best Savitzky-Golay Order: ', num2str(best_sg_order), ', Frame Length: ', num2str(best_sg_framelen)]);
disp(['Best Wavelet Denoising Level: ', num2str(best_wd_level)]);

% ------------------ End of Script ------------------
