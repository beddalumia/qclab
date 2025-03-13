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
	[q,r] = qr(z); % QR decomposition (algorithm-independent)
	d = diag(r);   % Extract the diagonal of R
	p = sign(d);   % Get the signs of D
	p(p==0) = 1;   % Numerical guardrail

	% Haar-measure (uniformly distributed U(N) matrix)
	U = bsxfun(@times,q,p.'); 

	% Final assertion
	if any(any(abs(inv(U)-U')>1E-12))
		error("Critical failure: the generated matrix is not unitary!")
	end

end