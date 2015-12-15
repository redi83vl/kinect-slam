function line = MEst(X, Y, eps,varargin)
%MEST line interpolation by M-Estimation
%   line = MEst(X, Y, eps [,L0] )
%       X,Y - coordinates of points
%       eps - delta value for difference between lines
%       L0 - Not Required. Is first approximation of line
%            in form a*x + b*y + c = 0, a^2 + b^2 <> 0
%
%   RETURN: line - line in form: a*x + b*y + c = 0, a^2 + b^2 = 1


if size(X,1) ~= 1
    X = X';
end
if size(Y,1) ~= 1
    Y = Y';
end

global POINT;
POINT = [X; Y; ones([1 length(X)])];
global sigm;


% First approximation
if nargin == 4
    L(1,:) = varargin{1};
    if any(size(L(1,:)) ~= [1 3])
        L(1,:) = TLS(X,Y);
    else
        if sqrt(L(1,1) ^ 2 + L(1,2) ^ 2) ~= 0
            L(1,:) = L(1,:) /sqrt(L(1,1) ^ 2 + L(1,2) ^ 2);
        else
            L(1,:) = TLS(X,Y);
        end
    end
else
    L(1,:) = TLS(X,Y);
end
L(2,:) = L(1,:) + 1 + eps;

eps = abs(eps);

% Main Cycle
count = 1;
while( sum( (L(2,:) - L(1,:)) .^ 2 ) > eps && count < 100)
    
    L(2,:) = L(1,:);
    
    % Calculating SIGM for this iteration
    R = sort(L(1, :) * POINT);
    ind = round(length(R)/2);
    sigm = 1.4826 * R(ind)*2;

    % find new L(INE)
    L(1,:) = fminsearch(@Cond, L(1,:));
    L(1,:) = L(1,:) /sqrt(L(1,1) ^ 2 + L(1,2) ^ 2);
    
    count = count + 1;
end

fprintf('Iterations in MEst: %d\n', count);
line = L(1,:);
return;

function J = Cond(line)
global POINT;
global sigm;

% Normalization !!!
line = line / sqrt( line(1) ^ 2 + line(2) ^ 2 );
R2 = (line * POINT) .^ 2;
J = sum(R2 ./ (R2 + sigm^2));