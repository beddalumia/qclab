function answer = is_rdm(rho)
	%% Checks if a given matrix is a proper (reduced) density matrix
	%  - it checks if the trace is unitary
	%  - it checks if the matrix has nonnegative eigenvalues
	%  - it checks if the matrix is hermitian
	%
	%  Copyright: Gabriele Bellomia, 2025

	NN = size(rho);
	if NN(1)~=NN(2)
		error("The matrix is not even square!")
	end

	if abs(trace(rho)-1)>1E-12
		error("The matrix is not normalized!")
	end

	if abs(rho-conj(rho'))>1E-12
		error("The matrix is not hermitian!")
	end

	p = eig(rho);
	if any(abs(p(p<0))>1E-12)
		error("The matrix has negative eigenvalues!")
	end

	answer = true;

end