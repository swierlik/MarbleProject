% noise_reduction.m

% Load the data
load('Data/lorenzData.mat') % Contains 'sol', 't', 'dt'

% Extract x from sol
x_original = sol(:,1);

% %Add a weird point
% x_original(10000) = 150;

%--------Add Gaussian Noise to x_original--------
variance = 0.01; % Noise level
x_noisy = x_original + variance*randn(size(x_original)); % Add Gaussian noise to x0

% ------------------ Noise Reduction on Stochastic Data ------------------

%Universal parameters (from optimization results)

%r=5;-------------

% %for 0.01 noise
% window_size = 11; % Moving Average
% order = 3; % Savitzky-Golay Order
% framelen = 95; % Savitzky-Golay Frame Length
% wd_level = 10; % Wavelet Denoising Level
% wd_wavelet ='dmey'; %Wavelet Denoising Wavelet

% %for 0.1 noise
% window_size = 15; % Moving Average
% order = 2; % Savitzky-Golay Order
% framelen = 87; % Savitzky-Golay Frame Length
% wd_level = 8; % Wavelet Denoising Level
% wd_wavelet ='dmey'; %Wavelet Denoising Wavelet

% %for 1 noise
% window_size = 15; % Moving Average
% order = 2; % Savitzky-Golay Order
% framelen = 83; % Savitzky-Golay Frame Length
% wd_level = 5; % Wavelet Denoising Level
% wd_wavelet ='coif3'; %Wavelet Denoising Wavelet

%r=10;-----------------

% %for 0.01 noise
% window_size = 15; % Moving Average
% order = 4; % Savitzky-Golay Order
% framelen = 79; % Savitzky-Golay Frame Length
% wd_level = 4; % Wavelet Denoising Level
% wd_wavelet ='coif3'; %Wavelet Denoising Wavelet


%for 0.15 noise
window_size = 15; % Moving Average
order = 4; % Savitzky-Golay Order
framelen = 79; % Savitzky-Golay Frame Length
wd_level = 4; % Wavelet Denoising Level
wd_wavelet ='coif3'; %Wavelet Denoising Wavelet

% ----------------Apply noise reduction to x's-----------------
x_movmean = movmean(x_noisy, window_size); % Moving Average
x_sg = sgolayfilt(x_noisy, order, framelen); % Savitzky-Golay
x_wd = wdenoise(x_noisy, wd_level, 'Wavelet', wd_wavelet); % Wavelet Denoising

% ------------------ Generate HAVOK System Data ------------------
r=7; % Koopman rank
tspan = dt:dt:50;
[V_x, A_x, B_x, xReg, U_x, E_x] = getSystem(x_original, 100, r, dt, tspan);
[V_x2, A_x2, B_x2, xReg2,U_x2, E_x2] = getSystem(x_noisy, 100, r, dt, tspan);
[V_movmean, A_x_movmean, B_x_movmean, xReg_movmean, U_movmean, E_movmean] = getSystem(x_movmean, 100, r, dt, tspan);
[V_sg, A_x_sg, B_x_sg, xReg_sg, U_sg, E_sg] = getSystem(x_sg, 100, r, dt, tspan);
[V_wd, A_x_wd, B_x_wd, xReg_wd, U_wd, E_wd] = getSystem(x_wd, 100, r, dt, tspan);

% ------------------ Reconstruct and simulate all systems ------------------
L = 1:min(length(tspan), size(xReg, 1));

% Original
sys_x = ss(A_x, B_x, eye(r-1), 0*B_x);  % System matrices for x
[y_sim_x, t_sim_x] = lsim(sys_x, xReg(L, r), dt*(L-1), xReg(1, 1:r-1));

% Noisy
sys_x2 = ss(A_x2, B_x2, eye(r-1), 0*B_x2);  % System matrices for x
[y_sim_x2, t_sim_x2] = lsim(sys_x2, xReg2(L, r), dt*(L-1), xReg2(1, 1:r-1));

% Moving Average
sys_x_movmean = ss(A_x_movmean, B_x_movmean, eye(r-1), 0*B_x_movmean);  % System matrices for x
[y_sim_x_movmean, t_sim_x_movmean] = lsim(sys_x_movmean, xReg_movmean(L, r), dt*(L-1), xReg_movmean(1, 1:r-1));

% Savitzky-Golay
sys_x_sg = ss(A_x_sg, B_x_sg, eye(r-1), 0*B_x_sg);  % System matrices for x
[y_sim_x_sg, t_sim_x_sg] = lsim(sys_x_sg, xReg_sg(L, r), dt*(L-1), xReg_sg(1, 1:r-1));

% Wavelet Denoising
sys_x_wd = ss(A_x_wd, B_x_wd, eye(r-1), 0*B_x_wd);  % System matrices for x
[y_sim_x_wd, t_sim_x_wd] = lsim(sys_x_wd, xReg_wd(L, r), dt*(L-1), xReg_wd(1, 1:r-1));

% -------------- Reconstruct x(t) for Each Method --------------
% Ensure U matches Koopman mode rank
U_x = U_x(:, 1:r-1);  
U_x2 = U_x2(:, 1:r-1);  
U_movmean = U_movmean(:, 1:r-1);  
U_sg = U_sg(:, 1:r-1);  
U_wd = U_wd(:, 1:r-1);  

% Ensure E matches Koopman mode rank
E_x = E_x(1:r-1, 1:r-1);
E_x2 = E_x2(1:r-1, 1:r-1);
E_movmean = E_movmean(1:r-1, 1:r-1);
E_sg = E_sg(1:r-1, 1:r-1);
E_wd = E_wd(1:r-1, 1:r-1);

% Project simulated results back to original space using SVD modes U
x_reconstructed = U_x * E_x * y_sim_x.';
x_reconstructed_noisy = U_x2 * E_x2 * y_sim_x2.';
x_reconstructed_movmean = U_movmean * E_movmean * y_sim_x_movmean.';
x_reconstructed_sg = U_sg * E_sg * y_sim_x_sg.';
x_reconstructed_wd = U_wd * E_wd * y_sim_x_wd.';

x_reconstructed = dehankelize(x_reconstructed);
x_reconstructed_noisy = dehankelize(x_reconstructed_noisy);
x_reconstructed_movmean = dehankelize(x_reconstructed_movmean);
x_reconstructed_sg = dehankelize(x_reconstructed_sg);
x_reconstructed_wd = dehankelize(x_reconstructed_wd);

% -------------- Error Calculation for Reconstructed x(t) --------------
error_x = (x_original(L) - x_reconstructed(L)).^2;
error_x_noisy = (x_original(L) - x_reconstructed_noisy(L)).^2;
error_x_movmean = (x_original(L) - x_reconstructed_movmean(L)).^2;
error_x_sg = (x_original(L) - x_reconstructed_sg(L)).^2;
error_x_wd = (x_original(L) - x_reconstructed_wd(L)).^2;

% -------------- Plot Results for Comparison --------------

figure;
hold on;
%plot(t, x_original, 'k', 'LineWidth', 1.5); % Original
plot(t_sim_x, x_reconstructed(L), 'b', 'LineWidth', 1.2); % HAVOK Original
plot(t_sim_x2, x_reconstructed_noisy(L), 'r', 'LineWidth', 1.2); % Noisy HAVOK
plot(t_sim_x_movmean, x_reconstructed_movmean(L), 'g', 'LineWidth', 1.2); % Moving Avg
plot(t_sim_x_sg, x_reconstructed_sg(L), 'm', 'LineWidth', 1.2); % Savitzky-Golay
plot(t_sim_x_wd, x_reconstructed_wd(L), 'c', 'LineWidth', 1.2); % Wavelet
%legend('Original', 'HAVOK Original', 'Noisy HAVOK', 'Moving Avg', 'Savitzky-Golay', 'Wavelet');
legend('HAVOK Original', 'Noisy HAVOK', 'Moving Avg', 'Savitzky-Golay', 'Wavelet');

xlabel('Time');
ylabel('x(t)');
title('Reconstruction of x(t) using HAVOK for Different Methods');
hold off;


% ------------------ Plot Results ------------------




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

%Noisy attractor plot
figure;
plot3(V_x2(L,1), V_x2(L,2), V_x2(L,3), 'r', 'LineWidth', 1);
title('Noisy Delay Embedded Attractor');
xlabel('v_1'); ylabel('v_2'); zlabel('v_3');
grid on;


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



% %--------------- Plot all of the different reconstructed x's on one plot ---------------
% figure;
% plot(tspan(L), V_x(L,1), 'b', 'LineWidth', 1);
% hold on;
% plot(tspan(L), y_sim_x2(L,1), 'r', 'LineWidth', 1);
% plot(tspan(L), y_sim_x_movmean(L,1), 'g', 'LineWidth', 1);
% plot(tspan(L), y_sim_x_sg(L,1), 'm', 'LineWidth', 1);
% plot(tspan(L), y_sim_x_wd(L,1), 'k', 'LineWidth', 1);
% title('Reconstructed x');
% xlabel('Time');
% ylabel('x');
% legend('Original', 'Noisy', 'Moving Average', 'Savitzky-Golay', 'Wavelet');



%Plot separate figures for each noise reduction method

% Moving Average

Rs = [1, int8(r/2), r-1];
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

% -------------- Bar Chart: Sum of Errors --------------
% Compute the sum of squared errors for each method over the index range L
errorMeans = [mean(error_x(L)), mean(error_x_noisy(L)), mean(error_x_movmean(L)), ...
             mean(error_x_sg(L)), mean(error_x_wd(L))];
methods = {'HAVOK Original', 'Noisy HAVOK', 'Moving Average', 'Savitzky-Golay', 'Wavelet'};

figure;
bar(errorMeans);
set(gca, 'XTickLabel', methods, 'XTick',1:numel(methods));
ylabel('Mean of Squared Errors');
title('Mean of Errors for Each Method');

% -------------- Bar Chart: Percentage of Time in the Same Lobe --------------
% For each method, calculate the percentage of time that the sign of the 
% reconstructed x matches the sign of the original x.
sameLobe_orig = sum( ((x_original(L) > 0 & x_reconstructed(L) > 0) | (x_original(L) < 0 & x_reconstructed(L) < 0)) );
sameLobe_noisy = sum( ((x_original(L) > 0 & x_reconstructed_noisy(L) > 0) | (x_original(L) < 0 & x_reconstructed_noisy(L) < 0)) );
sameLobe_movmean = sum( ((x_original(L) > 0 & x_reconstructed_movmean(L) > 0) | (x_original(L) < 0 & x_reconstructed_movmean(L) < 0)) );
sameLobe_sg = sum( ((x_original(L) > 0 & x_reconstructed_sg(L) > 0) | (x_original(L) < 0 & x_reconstructed_sg(L) < 0)) );
sameLobe_wd = sum( ((x_original(L) > 0 & x_reconstructed_wd(L) > 0) | (x_original(L) < 0 & x_reconstructed_wd(L) < 0)) );

percentageSameLobe = 100 * [sameLobe_orig, sameLobe_noisy, sameLobe_movmean, sameLobe_sg, sameLobe_wd] / length(L);

figure;
bar(percentageSameLobe);
set(gca, 'XTickLabel', methods, 'XTick',1:numel(methods));
ylabel('Percentage (%)');
title('Percentage of Time Vectors Are in the Same Lobe');



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
