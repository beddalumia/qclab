%%  RANDOM_RDM  Generates a random density matrix
%   This function has one required argument:
%     DIM: the number of rows (and columns) of the density matrix
%
%   RHO = random_rdm(DIM) generates a random DIM-by-DIM density
%   matrix, distributed according to the Hilbert-Schmidt measure.
%
%   This function has three optional arguments:
%     RE (default 0)
%     K (default DIM)
%     DIST (default 'haar')
%
%   RHO = random_rdm(DIM,RE,K,DIST) generates a random density
%   matrix of rank <= K, distributed according to the distribution DIST. If
%   RE = 1 then all of its entries will be real. DIST must be one of:
%     'haar' or 'hs' (default) - Generate a larger  pure state according to
%                      Haar measure and trace out the extra dimensions.
%                      Sometimes called the Hilbert-Schmidt measure when
%                      K = DIM.
%     'bures'        - the Bures measure
%
%   URL: http://www.qetlab.com/RandomDensityMatrix

%   requires: random_unitary.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 10, 2014

function rho = random_rdm(dim,varargin)

% set optional argument defaults: re=0, k=dim, dist='hs'
%[re,k,dist] = opt_args({ 0, dim, 'haar' },varargin{:});
if nargin<2
    re = false;
else
    re = varargin{1};
end
if nargin<3
    k = dim;
else
    k = varargin{2};
end
if nargin<4
    dist = 'hs';
else
    dist = varargin{3};
end

% Haar/Hilbert-Schmidt measure
gin = randn(dim,k);
if(~re)
    gin = gin + 1i*randn(dim,k);
end
if(strcmpi(dist,'bures')) % Bures measure
    rUn = random_unitary(dim);
    if re
        rUn = real(rUn);
    end
    gin = (random_unitary(dim) + eye(dim))*gin;
end

rho = gin*gin';
rho = rho/trace(rho);