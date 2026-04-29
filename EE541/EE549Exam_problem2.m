clear; clc; close all;
N     = 10;
d_lam = 0.5;
u0    = 0.0; % broadside u = 0, theta = 90 degrees
n = (-(N-1)/2 : (N-1)/2)';
u = linspace(-1, 1, 6001);
theta = linspace(0, pi, 6001);

% Steering vector function
sv = @(uv) exp(1j * 2*pi * d_lam * n * uv);   
% Beam Pattern 
BP = @(w0, u_arr) abs(exp(-1j*2*pi*d_lam * u_arr(:) * n') * w0);
% Theta space Beam Pattern function 
BP_theta = @(w0, theta_arr) abs(exp(-1j*2*pi*d_lam * cos(theta_arr(:)) * n') * w0);

SLL_dB = -25;
SLL = 10^(-25/20);
p = (-(N-1)/2 : (N-1)/2)'; 

vT = sv(u0);     

% [A]_mn = sinc(2π/λ |p_m - p_n|) = sinc(2 |p_m - p_n|) 
% Using matlab sinc def: sinc(x) = sin(πx)/(πx)
A = zeros(N, N);
for m = 1:N
    for nn = 1:N
        A(m,nn) = sinc(2 *d_lam* abs(p(m) - p(nn)));   % eq. (3.302)
    end
end

% Split into r small sectors for tight control
r       = 40;                          % number of sectors per side
%u_sl +- 0.26 for main beam. 
u_sl    = [linspace(-1, -0.24, r+1)', linspace(0.24, 1, r+1)'];
sectors = [];
for s = 1:2
    for k = 1:r
        ui   = (u_sl(k,s) + u_sl(k+1,s)) / 2;   % sector center
        Di   = (u_sl(k+1,s) - u_sl(k,s)) / 2;   % sector delta
        sectors = [sectors; ui, Di];
    end
end

num_sectors = size(sectors, 1);

% [Q_i]_mn 
Q = zeros(N, N, num_sectors);
% Li initialization
L = zeros(80,1);
for i = 1:num_sectors
    ui = sectors(i,1);
    Di = sectors(i,2);

    %Side lobe level per sector
    L(i) = SLL^2 * 2 * Di;
    for m = 1:N
        for nn = 1:N
            dp = d_lam *(p(m) - p(nn));
            Q(m,nn,i) = exp(1j*2*pi*dp*ui) * 2*Di * sinc(2*Di*dp); % Q formed
        end
    end
end

% eq. (3.315)
lambda     = ones([80,1]);           % initial factor lambda = 1
delta_lambda = zeros([r*2,1]);

% Build A_Q eq. (3.311)
AQ = A;
temp = zeros(80,1);

for iter = 1:20

    % Compute weights — eq. (3.315)
    w    = (AQ \ vT) / (vT' * (AQ \ vT));

    if iter==1 || iter == 5
        if iter == 1
            figure(1);
        elseif iter ==5
            figure(3);
        end

        bp = BP_theta(w, theta);
        bp_dB = 20*log10(bp/max(bp));
        plot(theta/pi*180,bp_dB); hold on;
        yline(SLL_dB,'r--','LineWidth',1,'Label','SLL=-25dB','LabelHorizontalAlignment','left');
        grid on;
        ylim([-50,0]);
        xlabel('theta (degrees)'); ylabel('|B|');

        if iter == 1
            title('Beam Pattern Algorithm  first iteration');
        elseif iter ==5
            title('Beam Pattern Algorithm  5th iteration');
        end
        legend('|B|','SLL=−26dB');

        figure(2);
    end

    B = zeros(num_sectors,1);

    for i = 1:num_sectors  % Beam-pattern magnitude in Sector i
        B(i) = w' * Q(:,:,i) * w;
        
        % delta update per sector based on Li
        if abs(B(i)) > L(i)           
            delta_lambda(i) = 0.3 * lambda(i);
        else             
            delta_lambda(i) = 0;
        end
    end

    AQ = AQ + sum(Q .* reshape(delta_lambda, 1, 1, []), 3); % AQ update.
    lambda = lambda + delta_lambda; % lambda update

    bp = BP_theta(w, theta);
    BP_dB = 20* log10(bp / max(bp));
    plot(theta/pi*180,BP_dB); hold on;

    if(sum(delta_lambda)==0) % Sidelobe levels in every sector (L_i)
                             % Dropped to the correct level
        fprintf('Ended after %ith iteration\n', iter);
        break;
    end
end
yline(SLL_dB,'r--','LineWidth',1,'Label','SLL=-25dB','LabelHorizontalAlignment','left');
grid on;
xlabel('theta (degrees)'); ylabel('|B|');
ylim([-50,0]);
title('Beam Pattern Algorithm all iterations');
legend('|B|','SLL=−26dB');


figure(4);
bp = BP_theta(w, theta);
bp_dB = 20*log10(bp/max(bp));
plot(theta/pi*180,bp_dB); hold on;
yline(SLL_dB,'r--','LineWidth',1,'Label','SLL=-25dB','LabelHorizontalAlignment','left');
ylim([-50,0]);
grid on;
xlabel('theta (degrees)'); ylabel('|B|');
title('Beam Pattern Algorithm  last iteration');
legend('|B|','SLL=−26dB');
%% c tabulate weight vector

fprintf('== Final Weight Vector W ==\n');
for i= 1:N
    fprintf('|      w%i = %.5f       |\n',i-1, w(i));
end
fprintf('===========================\n')

%% Part d

B_theta = exp(-1j*pi * cos(theta(:)) * n') * w;
B_theta_dB = 20*log10(B_theta/max(B_theta));

figure(5);

plot(theta/pi*180,bp_dB, 'g--'); hold on; grid on;
plot(theta/pi*180,B_theta_dB,LineStyle=':');
yline(SLL_dB,'r--','LineWidth',1,'Label','SLL=-25dB','LabelHorizontalAlignment','left');
ylim([-50,0]);
grid on;
xlabel('theta (degrees)'); ylabel('|B|');
title('Beam Pattern Algorithm vs Final Weight vector');
legend('|B_{algorithm}|','|B_{weight?}|','SLL=−26dB');