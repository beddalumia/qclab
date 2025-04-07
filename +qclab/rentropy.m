function S = rentropy(A,B)
	%% Relative entropy of density matrix A with respect to density matrix B
	%  NB. It's a noncommutative function of A and B
	%
	%  S(A||B) = Tr[ A * ( log2(A) - log2(B) ) ]
	%
	%  Copyright: Gabriele Bellomia, 2025

	% Sanitize noisy zeros
	A(abs(A) < 1e-8) = 0; 
	B(abs(B) < 1e-8) = 0;

	if not(isequal(size(A),size(B)))
		error("The density matrices must have the same dimension!")
	end

	if not(qclab.math.is_rdm(A))
		error("The first argument is not a well defined density matrix")
	end

	if not(qclab.math.is_rdm(B))
		error("The second argument is not a well defined density matrix")
	end

	if and(isdiag(A),isdiag(B))
		% Fast exit for classical probabilities
		% > we just compute the Kullback-Leibler divergence
		a = diag(A); b = diag(B);
		S = KL_divergence(a,b);
		return
	end

	if isequal(A*B,B*A)
		% If A and B commute we can diagonalize simultaneously
		[U,a] = eig(A,'vector');
		b = diag(U*B*U');
		% and just compute the Kullback-Leibler divergence
		S = KL_divergence(a,b);
		return
	end

	% Otherwise we gotta diagonalize separately (slowest exit)
	%[U,a] = eig(A,'vector'); a = real(a);
	%[V,b] = eig(B,'vector'); b = real(b);
	%log_A = U'*diag(log2(a))*U;
	%log_B = V'*diag(log2(b))*V;
	%A_ = U'*diag(a)*U; % Eigenvector sorting can be a bitch!
	% And finally get our /_quantum_/ relative entropy
	%S = trace( A_ * (log_A - log_B) );
    S = trace(A*(logm(A)-logm(B)))/log(2);
	% Assertion
	if S<0
		warning("Something bad happened: relative entropy is negative!")
		S = NaN;
	end
	return

end

function D = KL_divergence(p,q)
	%% Kullback-Leibler divergence 
	%  of two classical probability
	%  distributions {p} and {q}
	p = p(p>0); % we need to assume
	q = q(p>0); % p*log(p)=0 if p=0
	D = sum(p.*(log2(p) - log2(q)));
end