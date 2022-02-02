%--------------------------------------------------------------------------
% PRONY METHOD
%
% Allows real and imaginary frequencies
% Does not work in the presence of noise
%--------------------------------------------------------------------------
% INPUTS
% P         - System order
% x         - Data obtained experimentally
%--------------------------------------------------------------------------
% OUTPUTS
% z         - Estimated complex poles for the system
%--------------------------------------------------------------------------
function z = prony_method(P, x)

%--------------------------------------------------------------------------
% INITIALIZATION

% Number of samples
N = length(x);

% Initialization of matrices and vectors
H = zeros(N - P, P);
xa = zeros(N - P, 1);

%--------------------------------------------------------------------------
% DETERMINATION OF THE COEFFICIENTS OF THE CHARACTERISTIC POLYNOMIUM

% Construction of matrices and vectors
for i = 1:(N-P)
    for j = 1:P
        H(i,j) = x(P + i - j);
    end
    
    xa(i) = x(P + i);
end

% Calculation of the coefficients of the characteristic polynomial
% a = -(inv(H'*H)*H')*xa;
a = -pinv(H)*xa;

%--------------------------------------------------------------------------
% DETERMINATION OF THE SYSTEM POLES

% Calculation of the poles (it is necessary to add the element a(0) = 1)
z = roots([1; a]);

end
%--------------------------------------------------------------------------
