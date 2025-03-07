function S = rentropy(A,B)
	%% Relative entropy of density matrix A with respect to density matrix B
	%  NB. It's a noncommutative function of A and B
	%
	%  S(A||B) = Tr[ A * ( log2(A) - log2(B) ) ]
	%
	%  Copyright: Gabriele Bellomia, 2025

	if not(isequal(size(A),size(B)))
		error("The density matrices must have the same dimension!")
	end

	if and(isdiag(A),isdiag(B))
		% Fast exit for classical probabilities
		% > we just compute the Kullback-Leibler divergence
		a = diag(A); b = diag(B);
		S = KL_divergence(a,b);
		%return
	end

	if isequal(A*B,B*A)
		% If A and B commute we can diagonalize simultaneously
		[U,a] = eig(A,'vector');
		b = diag(U*B*conj(U'));
		% and just compute the Kullback-Leibler divergence
		S = KL_divergence(a,b);
		%return
	end

	% Otherwise we gotta diagonalize separately (slowest exit)
	[U,a] = eig(A,'vector');
	[V,b] = eig(B,'vector');
	log_A = conj(U')*diag(log2(a))*U;
	log_B = conj(V')*diag(log2(b))*V;
	A = conj(U')*diag(a)*U; % Eigenvector sorting can be a bitch!
	% And finally get our /_quantum_/ relative entropy
	S = trace( A * (log_A - log_B) );
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