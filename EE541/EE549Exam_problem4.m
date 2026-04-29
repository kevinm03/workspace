clc; clear; close all;

N   = 10;
d   = 0.5;
n   = (-(N-1)/2 : (N-1)/2)';
u   = linspace(-1, 1, 1000);

vs = @(u0) exp(1j*2*pi*d*n*u0);
BP = @(w0, u_arr) abs(exp(-1j*2*pi*d * u_arr(:) * n') * w0);

% Signal correlation matrix
S_F = [1.0 0.5 0.2;
       0.5 1.0 0.5;
       0.2 0.5 1.0];

% Scenarios: each row = [u_f1, u_f2, u_f3]
% f3 is always broadside (u=0), f1/f2 are interferers
ufi = [-0.3,  0.3,  0.0;
       -0.18, 0.18, 0.0;
       -0.02, 0.02, 0.0];

% Noise levels
sigma_w2 = [0, 10^(-20/10)]; 

% Build manifold matrices V(:,:,i) for each scenario
V = zeros(N, 3, 3);
for i = 1:3
    for k = 1:3
        V(:,k,i) = vs(ufi(i,k));
    end
end

colors    = {'b','r','k'};
noise_lbl = {'sigma^2=0', 'sigma^2=-20dB'};

for k = 1:2   % noise case
    S_n = sigma_w2(k) * eye(N);
    figure(); hold on; grid on;
    title(sprintf('H_3 Beam Pattern | %s', noise_lbl{k}));
    xlabel('u-space'); ylabel('|Beam-Pattern| dB');

    for i = 1:3   % scenario u_fi = .3 or .18 or .02
        R_x = V(:,:,i) * S_F * V(:,:,i)' + S_n;   % [N x N]
        H0  = S_F * V(:,:,i)' / R_x;            % [3 x N]
        H3  = H0(3,:).';                      % [N x 1] weights for f3

        bp = BP(H3, u);
        bp_n = bp / max(bp);
        bp_dB = 20*log10(bp_n);

        plot(u, bp_dB, colors{i}, 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Beam-Pattern for u_{f1,2}=\\pm%.2f', abs(ufi(i,1))));

        % Mark interferer and desired locations
        xline(ufi(i,1), '--', 'Color', colors{i}, 'Alpha', 0.4, ...
            'DisplayName', sprintf('u_{f1}=%.2f', ufi(i,1)));
        xline(ufi(i,2), '--', 'Color', colors{i}, 'Alpha', 0.4, ...
            'DisplayName', sprintf('u_{f2}=%.2f', ufi(i,2)));
    end
    xline(0, 'k:', 'LineWidth', 1.5, 'DisplayName', 'u_0=0 (desired)');
    ylim([-60 0]); legend('Location','southeast');
end