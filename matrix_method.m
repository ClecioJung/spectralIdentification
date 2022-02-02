%--------------------------------------------------------------------------
% MATRIX METHOD - SINGLE CHANNEL
%
% Allows real and imaginary frequencies
% Works in the presence of noise
% Operates with one or more channels
%--------------------------------------------------------------------------
% INPUTS
% P         - System order
% x         - Data obtained experimentally
%--------------------------------------------------------------------------
% OUTPUTS
% z         - Estimated complex poles for the system
%--------------------------------------------------------------------------
function z = matrix_method(P, x)

%--------------------------------------------------------------------------
% INITIALIZATION

% Number of samples and channels
[N, Q] = size(x);

% Determine the dimension of the Hankel matrix
L = ceil((N + 1)/2);
M = N - L + 1;

% Initialization of matrices and vectors
H = zeros(L, M*Q);

%--------------------------------------------------------------------------
% DETERMINATION OF THE HANKEL MATRIX

% Construction of the Hankel matrix
for q = 1:Q
    for i = 1:L
        for j = 1:M
            H(i, j + M*(q-1)) = x(i + j - 1, q);
        end
    end
end

%--------------------------------------------------------------------------
% DETERMINATION OF THE DECOMPOSITION OF H

% Determine the eigenvalues of H in descending form
[U, ~, ~] = svd(H);     % H = U*D*G'

% Select the first N eigenvectors of U
U = U(:, 1:P);

% Determination of matrices S1 and S2
U1 = U(1:L-1, :);
U2 = U(2:L, :);

%--------------------------------------------------------------------------
% DETERMINATION OF THE SYSTEM POLES

% phi determination
% phi = (inv(U1'*U1)*U1')*U2;
phi = pinv(U1)*U2;

% Calculation of the poles
z = eig(phi);

end
%--------------------------------------------------------------------------