function [Q,chi,C] = discord(rho,A,B)
	%% Compute Quantum Discord and Classical Correlations
	%  in the information geometry framework, assuming a
	%  metric induced by the Quantum Relative Entropy.
	%
	%  [Q,C,chi] = discord(rho)
	%
	%  with
	%
	%        rho: given density matrix (mixed or pure)
	%        A,B: bipartite reduced density matrices
	%        chi: closest pseudo-classical state
	%        Q:   quantum discord S(rho||chi)
	%        C:   classical correlations S(chi||A⊗B)
	%
	%  A and B are the reduced density matrices of rho, 
	%  associated to a given bipartition. NB: For now we
	%  assume a fixed bipartition, we'll see how to get
	%  an ergonomic generalized API in the future.
	%
	%  Copyright: Gabriele Bellomia, 2025
	%  
	%  Ref: Lexin Ding et al 2023 Quantum Sci. Technol. 8 015015

	N = size(rho,1); 
	if mod(N,4)~=0
		error("Is this a Fermionic density matrix?!")
	end

	try 
		qclab.math.is_rdm(rho);
	catch
		error("The given rho is not a valid density matrix!")
	end

	Nr = sqrt(N); % Reduced System Dimension

	% Initial guess for the local unitaries
	Ua = eye(Nr);  
	Ub = eye(Nr);

    % Initial guess for the chi
    chi = zeros(N,N);
    for i = 0:Nr-1
    	for j = 0:Nr-1
    		ixj = 1 + i + j*Nr;
    		chi(ixj,ixj) = rho(ixj,ixj);
    	end
    end
	try 
		qclab.math.is_rdm(chi);
	catch
		error("The initial guess for the closest pseudo-classical state is not a valid density matrix")
	end

    % Initial guess for the Discord
    Q = qclab.rentropy(rho,chi);

    M = 10; beta = 1;
    t = 0; converged = false;
    chi_ = zeros(N,N);

    converged = false;
    while not(converged)
    	% Sample local unitary matrices
    	Va = qclab.math.random_unitary(Nr);
    	Vb = qclab.math.random_unitary(Nr);
    	% Constrain them to be close to eye
    	dVa = mpower(Va,1/M); 
    	dVb = mpower(Vb,1/M); 
    	% Update local bases
    	Ua_ = dVa*Ua;
    	Ub_ = dVb*Ub;
    	% Map the density matrix to the new basis
    	Uaxb = kron(Ua_,Ub_);
    	rho_ = Uaxb*rho*Uaxb';
    	% Update guess for chi in the new basis
	    for i = 0:Nr-1
    		for j = 0:Nr-1
			    ixj = 1 + i + j*Nr;
    			chi_(ixj,ixj) = rho_(ixj,ixj);
    		end
    	end
		try 
			qclab.math.is_rdm(chi_);
		catch
			error("The updated candidate closest pseudo-classical state is not a valid density matrix")
		end
    	% Rotate back to the computational basis
    	% chi_ = Uaxb'*chi_*Uaxb; % DOING IT HERE CAUSES ALL SORTS OF NUMERICAL ISSUES
    	% Update guess for the Quantum Discord
    	Q_ = qclab.rentropy(rho_,chi_);
    	% Thermal shananingans
    	p = rand;
    	if exp(-beta*(Q_-Q)) > p
	    	% Convergence criteria
	    	if abs(Q_-Q)<1E-7
	    		converged = true;
	    	end
	    	if t > 10^4
	    		converged = true;
	    	end
    		Q = Q_;
    		chi = Uaxb'*chi_*Uaxb; % (rotate back to original basis)
    		Ua = Ua_;
    		Ub = Ub_;
    	end
    	t = t + 1;
    	beta = beta + 100; % Schedule :|

	end

	% Once we finish, we can compute the classical correlations
	% contained in the candidate closest pseudo-classical state
	C = qclab.mutinfo(chi,A,B); % More reliable for diagonal rho
	c = qclab.rentropy(chi,kron(A,B));
	if abs(C-c)>1E-12
		warning("The closest pseudo-classical state might be pathological.");
	end

	% Final sanity check and return
	assert(abs(imag(Q))<1E-12,"Imaginary discord!"); Q = real(Q);
	assert(abs(imag(C))<1E-12,"Imaginary classical correlations!"); C = real(C);

end