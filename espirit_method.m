%--------------------------------------------------------------------------
% ESPRIT METHOD
%
% Allows only imaginary frequencies
% Works in the presence of noise
%--------------------------------------------------------------------------
% INPUTS
% P         - System order
% x         - Data obtained experimentally
%--------------------------------------------------------------------------
% OUTPUTS
% z         - Estimated complex poles for the system
%--------------------------------------------------------------------------
function z = esprit_method(P, x)

%--------------------------------------------------------------------------
% INITIALIZATION

% Number of samples
N = length(x);

% Dimension of the covariance matrix
M = N;

% Initialization of matrices and vectors
R = zeros(M, M);

%--------------------------------------------------------------------------
% DETERMINATION OF THE COVARIANCE MATRIX

% Construction of matrices and vectors
for k = M:N
    R = R + (1/N)*x(k:-1:(k-M+1))*x(k:-1:(k-M+1))';
end

%--------------------------------------------------------------------------
% DETERMINATION OF THE DECOMPOSITION OF R

% Determine the eigenvalues of R in decreasing form
[U, ~, ~] = svd(R);     % R = U*D*G'

% Select the first N eigenvectors of U
S = U(:, 1:P);

% Determination of matrices S1 and S2
S1 = S(1:M-1, :);
S2 = S(2:M, :);

%--------------------------------------------------------------------------
% DETERMINATION OF THE SYSTEM POLES

% phi determination
% phi = (inv(S1'*S1)*S1')*S2;
phi = pinv(S1)*S2;

% Calculation of the poles
z = eig(phi);

end
%--------------------------------------------------------------------------