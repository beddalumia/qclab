%%  EENTROPY Computes the von Neumann or Renyi entropy of a density matrix
%   This function has one required argument:
%     RHO: a density matrix
%
%   S = eentropy(RHO) is the (base 2) von Neumann entropy of RHO.
%
%   This function has two optional input arguments:
%     BASE (default 2)
%     ALPHA (default 1)
%
%   S = eentropy(RHO,BASE,ALPHA) is the entropy of RHO, computed with
%   logarithms in the base specified by BASE. If ALPHA = 1 then this is the
%   von Neumann entropy. If ALPHA <> 1 then this is the Renyi-ALPHA
%   entropy.
%
%   copyright: kindly provided by http://www.qetlab.com/Entropy under BSD-2-Clause

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: May 12, 2016

function ent = eentropy(rho,varargin)

% set optional argument defaults: base=2, alpha=1
% [base,alpha] = opt_args({ 2, 1 },varargin{:});
if nargin < 2
    base = 2;
else
    base = varargin{1};
end
if nargin < 3
    alpha = 1;
else
    alpha = varargin{2};
end 

% Input assertion
if not(qclab.math.is_rdm(rho))
    warning("The input variable is not a well defined density matrix")
end

lam = eig(full(rho));
lam = lam(lam>0); % handle zero entries better: we want 0*log(0) = 0, not NaN

% If alpha == 1, compute the von Neumann entropy
if(abs(alpha - 1) <= eps^(3/4))
    if(base == 2)
        ent = -sum(real(lam.*log2(lam)));
    else
        ent = -sum(real(lam.*log(lam)))/log(base);
    end
elseif(alpha >= 0)
    if(alpha < Inf) % Renyi-alpha entropy with ALPHA < Inf
        ent = log(sum(lam.^alpha))/(log(base)*(1-alpha));

        % Check whether or not we ran into numerical problems due to ALPHA
        % being large. If so, compute the infinity-entropy instead.
        if(ent == Inf)
            alpha = Inf;
            warning('Entropy:LargeAlpha','Numerical problems were encountered due to a large value of ALPHA. Computing the entropy with ALPHA = Inf instead.');
        end
    end
    
    % Do not merge the following if statement with the previous one: we
    % need them separate, since this one catches a warning from the
    % previous block.
    if(alpha == Inf) % Renyi-infinity entropy
        ent = -log(max(lam))/log(base);
    end
else
    error('Entropy:InvalidAlpha','ALPHA must be non-negative.');
end