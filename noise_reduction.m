% noise_reduction.m

% Load the data
load('Data/lorenzData.mat') % Contains 'sol', 't', 'dt'
load('Data/systemData.mat') % Contains 'V_x', 'V_y', 'A_x', 'B_x', 'A_y', 'B_y', 'xReg', 'yReg', 'r', 'tspan', 'dt', 'S_x'

load('Data/lorenzDataStochastic.mat') % Contains 'sol2', 't2', 'dt2'
load('Data/systemDataStochastic.mat') % Contains 'V_x2', 'V_y2', 'A_x2', 'B_x2', 'A_y2', 'B_y2', 'xReg2', 'yReg2', 'r2', 'tspan2', 'dt2', 'S_x2'

% Extract x, y, z from sol
x_original = sol(:,1);
y_original = sol(:,2);
z_original = sol(:,3);

% Extract x, y, z from sol2
x_noisy = sol2(:,1);
y_noisy = sol2(:,2);
z_noisy = sol2(:,3);

% Define L (exclude initial and final transients)
L = 300:length(xReg)-300;

% ------------------ Noise Reduction on Stochastic Data ------------------

%Universal parameters (from optimization results)
%for 0.0001 noise
window_size = 3; % Moving Average
order = 2; % Savitzky-Golay Order
framelen = 45; % Savitzky-Golay Frame Length
wd_level = 1; % Wavelet Denoising Level
wd_wavelet ='bior3.5'; %Wavelet Denoising Wavelet
% %for 0.02 noise
% window_size = 7; % Moving Average
% order = 2; % Savitzky-Golay Order
% framelen = 49; % Savitzky-Golay Frame Length
% wd_level = 1; % Wavelet Denoising Level
% wd_wavelet ='dmey'; %Wavelet Denoising Wavelet
% %for 0.05 noise
% window_size = 13; % Moving Average
% order = 2; % Savitzky-Golay Order
% framelen = 23; % Savitzky-Golay Frame Length
% wd_level = 9; % Wavelet Denoising Level
% wd_wavelet ='sym4'; %Wavelet Denoising Wavelet
%for 0.3 noise
% window_size = 13; % Moving Average
% order = 3; % Savitzky-Golay Order
% framelen = 39; % Savitzky-Golay Frame Length
% wd_level = 6; % Wavelet Denoising Level


% Apply noise reduction to x's
% Apply noise reduction to x's
x_movmean = movmean(x_noisy, window_size); % Moving Average
%x_movmean(1:window_size) = x_original(1:window_size); % Set first values to original

x_sg = sgolayfilt(x_noisy, order, framelen); % Savitzky-Golay
%x_sg(1:framelen) = x_original(1:framelen); % Set first values to original

x_wd = wdenoise(x_noisy, wd_level, 'Wavelet', wd_wavelet); % Wavelet Denoising

% ------------------ Generate HAVOK System Data ------------------
[V_movmean, A_x_movmean, B_x_movmean, xReg_movmean, S_movmean] = getSystem(x_movmean, 100, r, dt, tspan);
[V_sg, A_x_sg, B_x_sg, xReg_sg, S_sg] = getSystem(x_sg, 100, r, dt, tspan);
[V_wd, A_x_wd, B_x_wd, xReg_wd, S_wd] = getSystem(x_wd, 100, r, dt, tspan);

% ------------------ Reconstruct and simulate all systems ------------------
L = 1:min(length(tspan), size(xReg, 1));

% Original
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);  % System matrices for x
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% Noisy
sys_x2 = ss(A_x2, B_x2, eye(r-1), 0*B_x2);  % System matrices for x
[y_sim_x2, t_sim_x2] = lsim(sys_x2, xReg2(L, r2), dt2*(L-1), xReg2(1, 1:r2-1));

% Moving Average
sys_x_movmean = ss(A_x_movmean, B_x_movmean, eye(r-1), 0*B_x_movmean);  % System matrices for x
[y_sim_x_movmean, t_sim_x_movmean] = lsim(sys_x_movmean, xReg_movmean(L, r), dt*(L-1), xReg_movmean(1, 1:r-1));

% Savitzky-Golay
sys_x_sg = ss(A_x_sg, B_x_sg, eye(r-1), 0*B_x_sg);  % System matrices for x
[y_sim_x_sg, t_sim_x_sg] = lsim(sys_x_sg, xReg_sg(L, r), dt*(L-1), xReg_sg(1, 1:r-1));

% Wavelet Denoising
sys_x_wd = ss(A_x_wd, B_x_wd, eye(r-1), 0*B_x_wd);  % System matrices for x
[y_sim_x_wd, t_sim_x_wd] = lsim(sys_x_wd, xReg_wd(L, r), dt*(L-1), xReg_wd(1, 1:r-1));

% % ------------------ Compute Errors ------------------

% ------------------ Compute Singular Value Weighted Error (SVWE) ------------------
S_norm = S_x;
S_noisy_norm = S_x2;
S_movmean_norm = S_movmean;
S_sg_norm = S_sg;
S_wd_norm = S_wd;

% Weight errors using singular values
error_x = sqrt(sum(S_norm(1:r-1) .* (V_x(L,1:r-1) - y_sim_x(L,1:r-1)).^2, 2));
error_x_noisy = sqrt(sum(S_noisy_norm(1:r-1) .* (V_x2(L,1:r-1) - y_sim_x2(L,1:r-1)).^2, 2));
error_x_movmean = sqrt(sum(S_movmean_norm(1:r-1) .* (V_movmean(L,1:r-1) - y_sim_x_movmean(L,1:r-1)).^2, 2));
error_x_sg = sqrt(sum(S_sg_norm(1:r-1) .* (V_sg(L,1:r-1) - y_sim_x_sg(L,1:r-1)).^2, 2));
error_x_wd = sqrt(sum(S_wd_norm(1:r-1) .* (V_wd(L,1:r-1) - y_sim_x_wd(L,1:r-1)).^2, 2));


% error_x = sqrt(sum((V_x(L,1:r-1) - y_sim_x(L,1:r-1)).^2, 2));
% error_x_noisy = sqrt(sum((V_x(L,1:r-1) - y_sim_x2(L,1:r-1)).^2, 2));
% error_x_movmean = sqrt(sum((V_x(L,1:r-1) - y_sim_x_movmean(L,1:r-1)).^2, 2));
% error_x_sg = sqrt(sum((V_x(L,1:r-1) - y_sim_x_sg(L,1:r-1)).^2, 2));
% error_x_wd = sqrt(sum((V_x(L,1:r-1) - y_sim_x_wd(L,1:r-1)).^2, 2));

%Compute errors weighted based off te weights in the system

% %RMSE for each x compare to original
% rmse_noisy = sqrt(mean((x_original(L) - x_noisy(L)).^2));
% rmse_movmean = sqrt(mean((x_original(L) - x_movmean(L)).^2));
% rmse_sg = sqrt(mean((x_original(L) - x_sg(L)).^2));
% rmse_wd = sqrt(mean((x_original(L) - x_wd(L)).^2));

% disp(['RMSE for Noisy: ', num2str(rmse_noisy)]);
% disp(['RMSE for Moving Average: ', num2str(rmse_movmean)]);
% disp(['RMSE for Savitzky-Golay: ', num2str(rmse_sg)]);
% disp(['RMSE for Wavelet Denoising: ', num2str(rmse_wd)]);

% %errors based off x at spots 1-3
% error_x = sqrt(sum((V_x(L,1:3) - y_sim_x(L,1:3)).^2, 2));
% error_x_noisy = sqrt(sum((V_x(L,1:3) - y_sim_x2(L,1:3)).^2, 2));
% error_x_movmean = sqrt(sum((V_x(L,1:3) - y_sim_x_movmean(L,1:3)).^2, 2));
% error_x_sg = sqrt(sum((V_x(L,1:3) - y_sim_x_sg(L,1:3)).^2, 2));
% error_x_wd = sqrt(sum((V_x(L,1:3) - y_sim_x_wd(L,1:3)).^2, 2));

% ------------------ Plot Results ------------------

% %PLot the RMSEs in a bar chart
% figure;
% bar([rmse_noisy, rmse_movmean, rmse_sg, rmse_wd]);
% title('RMSE for each denoising method');
% xlabel('Denoising Method');
% ylabel('RMSE');
% set(gca, 'xticklabel', {'Noisy', 'Moving Average', 'Savitzky-Golay', 'Wavelet'});
% grid on;

% % Figure 1: Original and noisy forcing vectors
% figure;
% subplot(2, 2, 1);
% plot3(x_original, y_original, z_original, 'b', 'LineWidth', 1);
% hold on;
% plot3(x_noisy, y_noisy, z_noisy, 'r', 'LineWidth', 1);
% title('Original vs Noisy System');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% legend('Original', 'Noisy');
% grid on;

% subplot(2, 2, 2);
% plot(x_noisy, 'r', 'LineWidth', 1);
% title('Noisy Data');
% xlabel('Time');
% ylabel('X');
% grid on;

% subplot(2, 2, 3);
% plot(x_movmean, 'g', 'LineWidth', 1);
% title('Denoised (Moving Average)');
% xlabel('Time');
% ylabel('X');
% grid on;

% subplot(2, 2, 4);
% plot(x_sg, 'm', 'LineWidth', 1);
% hold on;
% plot(x_wd, 'k', 'LineWidth', 1);
% title('Denoised (SG & Wavelet)');
% xlabel('Time');
% ylabel('X');
% legend('Savitzky-Golay', 'Wavelet');
% grid on;

% % Figure: Original vs Noisy vs each of the denoised x's on one plot
% figure;
% plot(tspan, x_original, 'b', 'LineWidth', 1);
% hold on;
% plot(tspan, x_noisy, 'r', 'LineWidth', 1);
% plot(tspan, x_movmean, 'g', 'LineWidth', 1);
% plot(tspan, x_sg, 'm', 'LineWidth', 1);
% plot(tspan, x_wd, 'k', 'LineWidth', 1);
% title('Original vs Noisy vs Denoised x');
% xlabel('Time');
% ylabel('X');
% legend('Original', 'Noisy', 'Moving Average', 'Savitzky-Golay', 'Wavelet');
% grid on;

L = 1000:min(length(tspan), size(xReg, 1));

% Figure 2: Delay Embedded Attractors
figure;
subplot(2, 2, 1);
%plot3(V_x(L,1), V_x(L,2), V_x(L,3), 'b', 'LineWidth', 1);
% hold on;
plot3(y_sim_x2(L,1), y_sim_x2(L,2), y_sim_x2(L,3), 'r', 'LineWidth', 1);
title('Noisy Delay Embedded Attractor');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
% legend('Original', 'Noisy');
view(-15,65)
grid on;

subplot(2, 2, 2);
plot3(y_sim_x_movmean(L,1), y_sim_x_movmean(L,2), y_sim_x_movmean(L,3), 'g', 'LineWidth', 1);
title('Reconstructed (Moving Average)');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
view(-15,65)
grid on;

subplot(2, 2, 3);
plot3(y_sim_x_sg(:,1), y_sim_x_sg(:,2), y_sim_x_sg(:,3), 'm', 'LineWidth', 1);
title('Reconstructed (Savitzky-Golay)');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
view(-15,65)
grid on;

subplot(2, 2, 4);
plot3(y_sim_x_wd(:,1), y_sim_x_wd(:,2), y_sim_x_wd(:,3), 'k', 'LineWidth', 1);
title('Reconstructed (Wavelet)');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
view(-15,65)
grid on;

% Figure 3: Error Comparison
figure;
plot(tspan(L), error_x(L), 'b', 'LineWidth', 1);
hold on;
plot(tspan(L), error_x_noisy(L), 'r', 'LineWidth', 1);
plot(tspan(L), error_x_movmean(L), 'g', 'LineWidth', 1);
plot(tspan(L), error_x_sg(L), 'm', 'LineWidth', 1);
plot(tspan(L), error_x_wd(L), 'k', 'LineWidth', 1);
title('Error Comparison');
xlabel('Time');
ylabel('Error');
legend('Original','Noisy' , 'Moving Average', 'Savitzky-Golay', 'Wavelet');
grid on;

%Figure 4: Delay embedded attractors
figure;

subplot(2, 2, 1);
plot3(V_x(L,1), V_x(L,2), V_x(L,3), 'b', 'LineWidth', 1);
title('Original Delay Embedded Attractor');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
%legend('Original');
view(-15,65)
grid on;

subplot(2, 2, 2);
%plot3(y_sim_x2(L,1), y_sim_x2(L,2), y_sim_x2(L,3), 'r', 'LineWidth', 1);
plot3(V_x2(L,1), V_x2(L,2), V_x2(L,3), 'r', 'LineWidth', 1);
title('Noisy Delay Embedded Attractor');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
%legend('Noisy');
view(-15,65)
grid on;

subplot(2, 2, 3);
plot3(y_sim_x(L,1), y_sim_x(L,2), y_sim_x(L,3), 'b', 'LineWidth', 1);
title('Reconstructed Original');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
%legend('Original');
view(-15,65)
grid on;

subplot(2, 2, 4);
plot3(y_sim_x2(L,1), y_sim_x2(L,2), y_sim_x2(L,3), 'r', 'LineWidth', 1);
%plot3(V_x2(L,1), V_x2(L,2), V_x2(L,3), 'r', 'LineWidth', 1);
title('Reconstructed Noisy');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
%legend('Noisy');
view(-15,65)
grid on;



%--------------- Plot all of the different reconstructed x's on one plot ---------------
figure;
plot(tspan(L), V_x(L,1), 'b', 'LineWidth', 1);
hold on;
plot(tspan(L), y_sim_x2(L,1), 'r', 'LineWidth', 1);
plot(tspan(L), y_sim_x_movmean(L,1), 'g', 'LineWidth', 1);
plot(tspan(L), y_sim_x_sg(L,1), 'm', 'LineWidth', 1);
plot(tspan(L), y_sim_x_wd(L,1), 'k', 'LineWidth', 1);
title('Reconstructed x');
xlabel('Time');
ylabel('x');
legend('Original', 'Noisy', 'Moving Average', 'Savitzky-Golay', 'Wavelet');



%Plot separate figures for each noise reduction method

% Moving Average

Rs = [1,5,10];
figure;
for j = 1:3
    i = Rs(j);
    subplot(3, 1, j);
    plot(tspan(L), V_x(L,i), 'b', 'LineWidth', 1);
    hold on;
    plot(tspan(L), y_sim_x2(L,i), 'r', 'LineWidth', 1);
    plot(tspan(L), y_sim_x_movmean(L,i), 'g', 'LineWidth', 1);
    title(['Moving Average: Reconstructed v', num2str(i)]);
    xlabel('Time');
    ylabel(['x', num2str(i)]);
    legend ('Original', 'Noisy', 'Moving Average');
end

% Savitzky-Golay
figure;
for j = 1:3
    i = Rs(j);
    subplot(3, 1, j);
    plot(tspan(L), V_x(L,i), 'b', 'LineWidth', 1);
    hold on;
    plot(tspan(L), y_sim_x2(L,i), 'r', 'LineWidth', 1);
    plot(tspan(L), y_sim_x_sg(L,i), 'm', 'LineWidth', 1);
    title(['Savitzky-Golay: Reconstructed v', num2str(i)]);
    xlabel('Time');
    ylabel(['x', num2str(i)]);
    legend ('Original', 'Noisy', 'Savitzky-Golay');
end

% Wavelet Denoising
figure;
for j = 1:3
    i = Rs(j);
    subplot(3, 1, j);
    plot(tspan(L), V_x(L,i), 'b', 'LineWidth', 1);
    hold on;
    plot(tspan(L), y_sim_x2(L,i), 'r', 'LineWidth', 1);
    plot(tspan(L), y_sim_x_wd(L,i), 'k', 'LineWidth', 1);
    title(['Wavelet: Reconstructed v', num2str(i)]);
    xlabel('Time');
    ylabel(['x', num2str(i)]);
    legend ('Original', 'Noisy', 'Wavelet');
end

% %Plot all sgolayfilt x's multiple subplots
% figure;

% j=1;
% for i = [1,5,10]
%     subplot(3, 1, j);
%     plot(tspan(L), y_sim_x2(L,i), 'r', 'LineWidth', 1);
%     hold on;
%     plot(tspan(L), y_sim_x_movmean(L,i), 'g', 'LineWidth', 1);
%     plot(tspan(L), y_sim_x_sg(L,i), 'm', 'LineWidth', 1);
%     plot(tspan(L), y_sim_x_wd(L,i), 'k', 'LineWidth', 1);
%     plot(tspan(L), V_x(L,i), 'b', 'LineWidth', 1);
%     title(['Reconstructed v', num2str(i)]);
%     xlabel('Time');
%     ylabel(['x', num2str(i)]);
%     legend ('Noisy', 'Moving Average', 'Savitzky-Golay', 'Wavelet','Original');
%     j=j+1;
% end


% ------------------ End of Script ------------------
