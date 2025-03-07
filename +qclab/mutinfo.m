function I = mutinfo(rho,rho_A,rho_B)
	%% Mutual information of the state rho with respect to the A:B bipartion
	%
	%  I(A:B) = S(ρ_A) + S(ρ_B) - S(ρ) == S( ρ || ρ_A ⊗ ρ_B )
	%
	%  where ρ_A = trace_B ( ρ ), ρ_B = trace_A ( ρ )
	%  and S(ρ) denotes the von Neumann entropy of ρ
	%
	%  Copyright: Gabriele Bellomia, 2025
	%
	%  See also: qclab.eentropy
	%
	% NOTE: For now we are leaving out the issue of defining an ergonomic API
	%       to select the A and B bipartition, let's just assume that the user
	%       can do his traces and pass ρ, ρ_A and ρ_B. In the future we want 
	%       the function to perform the traces itself, of course.

	S = qclab.eentropy(rho);
	A = qclab.eentropy(rho_A);
	B = qclab.eentropy(rho_B);
	I = A + B - S;

end