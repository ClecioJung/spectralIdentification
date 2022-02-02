%--------------------------------------------------------------------------
% SPECTRAL IDENTIFICATION
% Author: Clecio Jung
% Date: 21/01/2018
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Inputs
% method    - Method chosen for determining the poles
% P         - System order
% Ts        - Sampling time
% x         - Data obtained experimentally
% z_k       - Previously known poles (if any)
%--------------------------------------------------------------------------
% Outputs
% z         - Estimated complex poles for the system
% h         - Estimated complex amplitude for the system
% tau       - Time constant for the estimated poles
% f         - Frequency for the estimated poles
% A         - Amplitude for the estimated poles
% theta     - Delay angle for the estimated pole
% x_est     - Data obtained by the estimated system
% y         - Data filtered by known poles
% err       - Error for the estimated parameters
% J         - Residual sum of squares (RSS)
% Jm        - Average RSS
%--------------------------------------------------------------------------

function [z, h, tau, f, A, theta, x_est, y, err, J, Jm] = spectralIdentification(method, P, Ts, x, z_k)

%--------------------------------------------------------------------------
% INITIALIZATION

% Tests if provided values are valid
if (P < 1)
    error('Error: P must be greater than or equal to one!');
end

% Number of samples and channels
[N, Q] = size(x);
if (N < 2*P)
    error('Error: N must be greater than or equal to 2P!');
end

% Number of previously known poles
z_k = z_k(:);
F = length(z_k);
if (F >= P)
    error('Error: All poles are already known!');
end

% Initialization of matrices and vectors
V = zeros(N, P);
alfa = zeros(P, 1);
tau = zeros(P, 1);
f = zeros(P, 1);
A = zeros(P, 1);
theta = zeros(P, 1);

%--------------------------------------------------------------------------
% STEP 1: PRE-FILTER

% Coefficients of the polynomial that defines the known poles
b = poly(z_k);

% Filters input data using polynomial coefficients
y = filter(b, 1, x);

%--------------------------------------------------------------------------
% STEP 2: DETERMINATION OF ''Z'' POLES BY THE SPECIFIED METHOD

if (strcmp(method, 'matrix'))
    % Operates with one or more channels
    z_est = matrix_method((P-F), y);
else
    if (Q ~= 1)
        error('Error: To use more than one channel, use the ''matrix'' method!');
    else
        if (strcmp(method, 'prony'))
            z_est = prony_method((P-F), y);
        elseif (strcmp(method, 'esprit'))
            z_est = esprit_method((P-F), y);
        else
            error('Error: Specified method not found!');
        end
    end
end

% Attach previously known poles
z = [z_k; z_est];

for i = 1:P
    % Determination of the time constant of the pole
    alfa(i) = log(abs(z(i)))/Ts;
    if (alfa(i) ~= 0)
        tau(i) = -1/alfa(i);
    end
    
    % Determination of the frequency of the pole
    f(i) = angle(z(i))/(2*pi()*Ts);
end

%--------------------------------------------------------------------------
% STEP 3: DETERMINATION OF RANGE AND ANGLE OF RESPONSE DELAY

% Construction of matrices and vectors
for i = 1:N
    for j = 1:P
        V(i, j) = z(j)^(i-1);
    end
end

% Calculation of complex amplitudes
% h = (inv(V'*V)*V')*x;
h = pinv(V)*x;

for q = 1:Q
    for i = 1:P
        % Determining the real amplitude
        A(i, q) = abs(h(i, q));
        
        % Determining the sinusoid delay angle
        theta(i, q) = angle(h(i, q));
    end
end

%--------------------------------------------------------------------------
% STEP 4: DETERMINATION OF ESTIMATED RESPONSE AND PERFORMANCE INDICATORS

if (isreal(x))
    x_est = real(V*h);
else
    x_est = V*h;
end

err = (x - x_est);
J = err'*err;
Jm = J/N;

end
%--------------------------------------------------------------------------