function U = random_unitary(N)
	%% Draws a random unitary matrix of rank N
	%  according to the Haar measure, a matrix
	%  equivalent of uniform distributions.
	% 
	%  Copyright: Francesco Mezzadri, Gabriele Bellomia, 2025
	%
	%  See:
	%  https://mathoverflow.net/questions/76295/intuition-for-haar-measure-of-random-matrix
	%  https://arxiv.org/abs/math-ph/0609050 

	z = (randn(N) + 1i*randn(N))/sqrt(2); % Ginibre Ensemble
	[q,r] = qr(z);   % QR decomposition (algorithm-independent)
	d = diag(r);     % Extract the diagonal of R
	p = sign(d);     % Get the signs of D
	q = q*diag(p)*q; % Select the correct "gauge" for Q

	U = q; % Haar-measure (uniformly distributed U(N) matrix)

end