function answer = is_rdm(rho)
	%% Checks if a given matrix is a proper (reduced) density matrix
	%  - it checks if the trace is unitary
	%  - it checks if the matrix has nonnegative eigenvalues
	%  - it checks if the matrix is hermitian
	%
	%  Copyright: Gabriele Bellomia, 2025

	answer = true;

	NN = size(rho);
	if NN(1)~=NN(2)
		warning("The matrix is not even square!")
		answer = false;
	end

	if abs(trace(rho)-1)>1E-12
		warning("The matrix is not normalized!")
		answer = false;
	end

	if abs(rho-conj(rho'))>1E-12
		warning("The matrix is not hermitian!")
		answer = false;
	end

	p = eig(rho);
	if any(abs(p(p<0))>1E-12)
		warning("The matrix has negative eigenvalues!")
		answer = false;
	end


end