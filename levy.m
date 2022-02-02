%--------------------------------------------------------------------------
% Levy's identification method
% Clecio Jung
% 06/08/2020
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% INPUTS
% m	- FT's numerator order
% n	- FT's denominator order
% w	- Frequencies vector (rad/s)
% h	- Amplitude vector
%--------------------------------------------------------------------------
% OUTPUTS
% a - FT's denominator coefficients
% b - FT's numerator coefficients
%--------------------------------------------------------------------------
function [a, b] = levy(m, n, w, h)

N = length(h);

f = zeros(2*N,1);
A1 = zeros(2*N,floor(n/2));
A2 = zeros(2*N,ceil(n/2));
B = zeros(2*N,m+1);

for i = 1:N
    f(2*i-1) = real(h(i));
    f(2*i) = imag(h(i));
    for j = 1:floor(n/2)
        A1(2*i-1,j) = (-1)^j*real(h(i))*w(i)^(2*j);
        A1(2*i,j) = (-1)^j*imag(h(i))*w(i)^(2*j);
    end
    for j = 1:ceil(n/2)
        A2(2*i-1,j) = (-1)^j*imag(h(i))*w(i)^(2*j-1);
        A2(2*i,j) = (-1)^(j+1)*real(h(i))*w(i)^(2*j-1);
    end
    for j = 1:ceil((m+1)/2)
        B(2*i-1,j) = (-1)^(j+1)*w(i)^(2*j-2);
    end
    for j = 1:floor((m+1)/2)
        B(2*i,j+ceil((m+1)/2)) = (-1)^(j+1)*w(i)^(2*j-1);
    end
end

K = [A1, A2, -B];

theta = -pinv(K)*f;
v = theta(1:floor(n/2));
w = theta(floor(n/2)+1:ceil(n/2)+floor(n/2));
y = theta(ceil(n/2)+floor(n/2)+1:end);

a = zeros(n,1);
b = zeros(m+1,1);
for i = 1:n
    if mod(i,2) == 0 % a2, a4, a6, ...  
        a(i) = v(floor((i+1)/2));
    else % a1, a3, a5, ...       
        a(i) = w(floor((i+1)/2));
    end
end
for i = 1:(m+1)
    if mod(i,2) == 0 % b1, b3, b5, ...  
        b(i) = y(floor(i/2)+ceil((m+1)/2));
    else % b0, b2, b4, ...        
        b(i) = y(floor((i+1)/2));
    end
end

end
%--------------------------------------------------------------------------
