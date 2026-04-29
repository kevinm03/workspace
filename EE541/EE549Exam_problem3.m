N = 10;
n = (-(N-1)/2 : (N-1)/2)';
d_lam = 0.5;
u = linspace(-1, 1, 3000);

%interference @ broadside u = 0;
ui = 0;

%spatial correlation coefficient between 0-1
% ps1 = vs^H * v1 / N
%steering vector function
sv = @(uv) exp(1j * 2*pi * n * uv);  

v1 = sv(ui);
vs = sv(u);

ps1 = vs' * v1 / N;
ps1_sqr = abs(ps1) .^2;

INR_dB = [-20 -10 0 10 20 30];
INR = 10.^(INR_dB/10);

A = N * (1 + INR.^2);
A = A .* (1 + N*INR.^2 .* (1 - ps1_sqr)) ./ (1 + N*INR.^2);
A_dB = 10 * log10(A);

% ---- Plotting ----
for k = 1:length(INR_dB)
    plot(1-ps1_sqr, A_dB(:,k), 'LineWidth', 2);
    hold on;
    grid on;
end
title('A_0 for varying spatial correlation coefficient');
legend('INR = -20dB', 'INR = -10 dB','INR = 0dB', 'INR = 10 dB', 'INR = 20 dB', 'INR = 30 dB');
xlabel('1-|ρ_s_1|^2');
ylabel('A_o dB');