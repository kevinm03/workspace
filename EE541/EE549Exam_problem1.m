clc;clear;close all;

u = linspace(-1,1,1000)';
j = sqrt(-1);
lambda = 1;
d = lambda/2;
k = -2*pi/lambda * u;


location = [0 d 4*d 6*d];
centered_location = location - (3*d);

V =  exp(j*centered_location.*k)';

% Assuming Even weights
W = ones([4,1]) * 1/4;
AF = W' * V;

%normalize
AF = AF ./ max(abs(AF));

%dB conversion
AF_dB = 20 * log10(abs(AF));

plot(u,AF_dB);
title('Non-Uniform 4 Element Array Factor');
xlabel('u');ylabel('|AF_{dB}|'); grid on;

% ----------- part b ------------------
%Null-to-null beamwidth
idx = find(islocalmin(real(AF_dB)));
left_null = idx(length(idx)/2 );
right_null = idx(length(idx)/2+1);

BW_NN = u(right_null) - u(left_null);

%Half power beamwidth
left3dB = u(find(AF_dB > -3, 1, "first"));
right3dB = u(499 + find(AF_dB(500:1000) < -3, 1, "first"));

HPBW = right3dB - left3dB;


fprintf("BW_nn = %.3f in u-space\n", BW_NN)
fprintf("HPBW = %.3f in u-space\n", HPBW)

% ----------- part c -------------------
location2 = [0 d 2*d 3*d 4*d 5*d 6*d];
centered_location2 = location2 - 3.5*d;

V2 =  exp(j*centered_location2.*k)';

% Assuming Even weights
W2 = ones([7,1]) * 1/7;

AF2 = W2' * V2;
%normalize
AF2 = AF2 ./ max(abs(AF2));

AF2_dB = 20 * log10(abs(AF2));

figure();
hold on;
plot(u,AF2_dB);
plot(u,AF_dB);
title('Non-Uniform 4 Element vs Uniform 7 Element');
subtitle('Array Factor');
xlabel('u');ylabel('|AF_{dB}|'); grid on;
legend(["uniform 7-element", "non-uniform 4-element"]);