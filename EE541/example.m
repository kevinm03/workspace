clear; clc; close all;

% Array parameters 
N      = 10;   % number of elements
d      = 0.5;   % half-wavelength
n      = (0:N-1)';   % element index

% Signal & Interference parameters
u_s = 0.0;   % desired signal 
u_I = [0.18, 0.3];   % interferer

SNR_dB = [0,10];
INR_dB = [0,10];
c=0;
for i=1:2
    for a=1:2
        for b=1:2
            c = c+1;
            sigma_s2(i) = 10^(SNR_dB(i)/10);   % signal power  (sigma_w^2 normalised to 1)
            sigma_I2(a) = 10^(INR_dB(a)/10);   % interferer power
            sigma_w2 = 1.0;               % white noise variance
            
            %  Steering vectors 
            v_s = exp(1j*2*pi*d*u_s*n);   % signal manifold
            v_I = exp(1j*2*pi*d*u_I(b)*n);   % interferer manifold
            
            % Covariance matrix and cross-correlation 
            S_x(:,:,c) = sigma_s2(i)*(v_s*v_s') + sigma_I2(a)*(v_I(b)*v_I(b)') + sigma_w2*eye(N);
            p   = v_s * sigma_s2(i);          % p = S_{xd*} = v_s * sigma_s^2
            
            % Optimal (Wiener-Hopf) solution 
            w_o(:,c) = S_x(:,:,c) \ p;  % w_o = S_x^{-1} p
            
            xi_o(c) = sigma_s2(i) - w_o(:,c)'*p;     % MMSE floor  (Eq. 7.324, real part)
            xi_o(c) = real(xi_o(c));
            
            % Step size 
            lambda(:,c) = eig(S_x(:,:,c));
            lambda_max(c) = max(real(lambda(:,c)));
            lambda_min(c) = min(real(lambda(:,c)));
            
            alpha(c) = 0.05 / lambda_max(c);    % conservative choice (~5% of stability limit)
            fprintf('lambda_max = %.4f,  lambda_min = %.4f\n', lambda_max(c), lambda_min(c));
            fprintf('alpha      = %.6f  (stability limit = %.6f)\n', alpha(c), 2/lambda_max(c));
            
            % Initial weight vector
            w = v_s / N;
            
            % Steepest descent iterations 
            K_max = 100;
            xi_sd = zeros(1, K_max+1);          % transient MSE
            
            % MSE at K=0  (Eq. 7.366)
            xi_sd(1) = real(xi_o(c) + (w - w_o(:,c))' * S_x(:,:,c) * (w - w_o(:,c)));
            
            for K = 1:K_max
                gradient = -p + S_x(:,:,c)*w;          % gradient of xi w.r.t. w^H (Eq. 7.322)
                w = w - alpha * gradient;        % steepest descent step (Eq. 7.325)
                xi_sd(K+1) = real(xi_o(c) + (w - w_o(:,c))' * S_x(:,:,c) * (w - w_o(:,c)));
            end
            
            w_final = w;   % converged weight vector
            
            % Time constants
            tau(:,c) = -1 ./ log(abs(1 - alpha*real(lambda(:,c))));   % exact
            tau_approx(:,c) = 1 ./ (alpha * real(lambda(:,c)));        % small-alpha approx
            fprintf('tau_min ≈ %.1f,  tau_max ≈ %.1f  (iterations)\n', ...
                    min(tau_approx(:,c)), max(tau_approx(:,c)));
            
            % Beam pattern 
            u_grid = linspace(-1, 1, 1001);
            B_o    = zeros(size(u_grid));
            B_K    = zeros(size(u_grid));
            
            for k = 1:numel(u_grid)
                va       = exp(1j*2*pi*d*u_grid(k)*n);
                B_o(k)   = w_o(:,c)'    * va;
                B_K(k)   = w_final' * va;
            end
            
            %% ── Figure 1: MSE convergence ────────────────────────────────────────────
            figure('Name','MSE Convergence','NumberTitle','off','Position',[100 500 700 320]);
            semilogy(0:K_max, xi_sd, 'b-', 'LineWidth', 1.5); hold on;
            semilogy([0 K_max], [xi_o(c) xi_o(c)], 'r--', 'LineWidth', 1.2);
            xlabel('Iteration K');
            ylabel('\xi(K)  (MSE)');
            title('Steepest Descent Convergence  (u_s=0, u_I=0.3, SNR=INR=10 dB)');
            legend('\xi_{sd}(K)', '\xi_o(c)  (MMSE floor)', 'Location','northeast');
            grid on;
            xlim([0 K_max]);
            
            %% ── Figure 2: Beam patterns ──────────────────────────────────────────────
            figure('Name','Beam Patterns','NumberTitle','off','Position',[100 150 700 380]);
            plot(u_grid, 20*log10(abs(B_o)+1e-10), 'r--', 'LineWidth', 1.5); hold on;
            plot(u_grid, 20*log10(abs(B_K)+1e-10), 'b-',  'LineWidth', 1.5);
            xline(u_s, 'k:', 'u_s', 'LabelVerticalAlignment','bottom');
            xline(u_I, 'm:', 'u_I', 'LabelVerticalAlignment','bottom');
            xlabel('u = sin(\theta)');
            ylabel('Beam pattern (dB)');
            title('Beam Pattern: Optimal vs. Steepest Descent (converged)');
            legend('w_o(:,c)  (Wiener-Hopf)', 'w_{SD}  (converged)', 'Location','southwest');
            ylim([-60 15]);
            grid on;
            
            %% ── Figure 3: Eigenvalue spectrum ────────────────────────────────────────
            figure('Name','Eigenvalues','NumberTitle','off','Position',[820 500 420 280]);
            stem(sort(real(lambda),'descend'), 'b', 'filled', 'MarkerSize', 5);
            xlabel('Eigenvalue index');  ylabel('\lambda_n');
            title('Eigenvalues of S_x(:,:,c)');
            grid on;
            
            %% ── Figure 4: Weight comparison (magnitude) ──────────────────────────────
            figure('Name','Weights','NumberTitle','off','Position',[820 150 420 280]);
            plot(1:N, abs(w_o(:,c)),    'rs-', 'LineWidth', 1.4, 'MarkerSize', 6); hold on;
            plot(1:N, abs(w_final),'bo-', 'LineWidth', 1.4, 'MarkerSize', 6);
            xlabel('Element n');  ylabel('|w_n|');
            title('Weight magnitudes');
            legend('w_o(:,c)','w_{SD}','Location','best');
            grid on;
            
            %% ── Summary ──────────────────────────────────────────────────────────────
            fprintf('\n=== Summary ===\n');
            fprintf('  MMSE floor  xi_o(c)       = %.6f\n', xi_o(c));
            fprintf('  Final MSE   xi_sd(end) = %.6f\n', xi_sd(end));
            fprintf('  ||w_SD - w_o(:,c)||         = %.2e\n', norm(w_final - w_o(:,c)));
            fprintf('  Eigenvalue spread      = %.2f  (lambda_max/lambda_min)\n', lambda_max/lambda_min);
        end
    end
end