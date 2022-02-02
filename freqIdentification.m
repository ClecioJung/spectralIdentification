%--------------------------------------------------------------------------
% Frequency Response Identification of a linear circuit
% Apply a voltage signal consisting of N frequencies to a linear circuit.
% Provide the voltage and current signals measured to this function
% and it should estimate the coefficients of circuit transfer function.
%
% Clecio Jung
% 06/08/2020
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% INPUTS
% N       - Number of frequencies in the voltage and current signals
% m	      - FT's numerator order
% n	      - FT's denominator order
% Ta      - Sampling time (s)
% voltage - Measured voltage vector
% current - Measured current vector
%--------------------------------------------------------------------------
% OUTPUTS
% a       - FT's denominator coefficients
% b       - FT's numerator coefficients
% f       - Frequency for the estimated poles
% h	      - Amplitude vector
% v_est   - Votlage obtained by the estimated system
% i_est   - Current obtained by the estimated system
%--------------------------------------------------------------------------
function [a, b, f, h, v_est, i_est] = freqIdentification(N, m, n, Ta, voltage, current)
epslon = 0.01;

[zv, hv, ~, fv, ~, ~, v_est, ~, ~, ~, ~, ~] = spectralIdentification('matrix', N, Ta, voltage, []);
[~, hi, ~, fi, ~, ~, ~, ~, ~, ~, ~, ~] = spectralIdentification('matrix', N, Ta, current, []);

f = fv(fv>=0);
h = zeros(length(f),1);
hminus = zeros(length(f),1);

for i = 1:length(f)
    h(i) = hi(abs(fi-f(i)) < epslon)/hv(abs(fv-f(i)) < epslon);
    hminus(i) = hi(abs(fi+f(i)) < epslon)/hv(abs(fv+f(i)) < epslon);
end

[a,b] = levy(m, n, 2*pi*f, h);

hi_est = zeros(length(hv),1);
for i = 1:length(hv)
    w = 2*pi*fv(i);
    s = 1i*w;
    num = 0;
    den = 1;
    for j = 1:length(b)
        num = num + b(j)*s^(j-1);
    end
    for j = 1:length(a)
        den = den + a(j)*s^j;
    end
    Y = num/den;
    hi_est(i) = Y*hv(i);
end
clear w s num den Y;

V = zeros(length(v_est), N);
for i = 1:length(v_est)
    for j = 1:N
        V(i, j) = zv(j)^(i-1);
    end
end

i_est = real(V*hi_est);
end
%--------------------------------------------------------------------------
