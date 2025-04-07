function [Q,chi,C] = discord(rho,A,B)
	%% Compute Quantum Discord and Classical Correlations
	%  in the information geometry framework, assuming a
	%  metric induced by the Quantum Relative Entropy.
	%
	%  [Q,chi,C] = discord(rho,A,B)
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
	%       (but we implement some additional tricks to minimize)

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
	Nc = 1000; % Number of candidates

	taboos = {};
	q_best = NaN;
	for instance = 1:Nc

        fprintf("Starting the walk #%d\n",instance)

		% Initializations
		chi = zeros(N,N);
		M = 1000+instance; beta = 1;

		% Initial guess for the local unitaries
        if instance==1

            % Always try the identity (other Haar starting points apparently
            %                          suck... to be investigated why)
            Ua = eye(Nr);  
            Ub = eye(Nr);
    
		    % Initial guess for the chi
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
        
        else

		    % In general try a random Haar-measure unitary
		    Ua = qclab.math.random_unitary(Nr);
		    Ub = qclab.math.random_unitary(Nr);
            %
		    % Map the density matrix to the new basis
		    Uaxb = kron(Ua,Ub);
		    rho_ = Uaxb*rho*Uaxb';
            %
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
            % 
		    % % Initial guess for the Discord
		    Q = qclab.rentropy(rho_,chi_);
            % 
		    % % Rotate back to the computational basis
		    chi = Uaxb'*chi_*Uaxb;

        end

		converged = false;
		t = 0; n = 0;
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
			% chi = Uaxb'*chi_*Uaxb; % DOING IT HERE CAUSES ALL SORTS OF NUMERICAL ISSUES
			% Update guess for the Quantum Discord
			Q_ = qclab.rentropy(rho_,chi_);
			% Rotate back to the computational basis
			chi = Uaxb'*chi_*Uaxb; 
			% Determine if the new chi is in the taboo list
			taboo = false;
			Dl = zeros(1,numel(taboos));
			Dr = zeros(1,numel(taboos));
			for i = 2:numel(taboos)
				Dl(i) = qclab.rentropy(chi,taboos{i}); % Rel. Entropy is 
				Dr(i) = qclab.rentropy(taboos{i},chi); % noncommutative!
			end
			taboo = any(Dl(i)<1e-8) || any(Dr(i)<1e-8);
			if taboo
				% We discard as if very high
                fprintf("TABOO: let's avoid entering a known false minimum\n")
				break
			else
				% Thermal shenanigans
				p = rand;
				if exp(-beta*(Q_-Q)) > p
					Q = Q_;
					Ua = Ua_;
					Ub = Ub_;
					n = n + 1;
					% Convergence criteria
					if abs(Q_-Q)<1E-7
						converged = true;
                        fprintf("Exiting walk #%d as we are converged!\n",instance)
					end
					if n > 10^4
						converged = true;
                        fprintf("Exiting walk #%d after too many steps!\n",instance)
					end
					if converged 
					% Taboo-enhanced criteria
						if q_best<Q
							% We have already a better candidate
							% so we discard Q and flag chi as a
							% taboo state
                            fprintf("TABOO: found a false minimum!\n")
                            fprintf("       > discarding the walk.\n")
							taboos = {taboos,chi};
                            break
							converged = false;
                            t = 0; n = 0;
						else
							% The new converged Q is a good candidate,
							% so we do not update the taboo list and
							% we record the new best candidate
							q_best = Q;
                            chi_best = chi;
							converged = true;
						end	
					end
				end
			end
			t = t + 1;
			beta = beta + 100; % Schedule :|
			if mod(t,100) == 0
				fprintf("Time step #%d\n---------------------\n",t)
				fprintf("Annealing temperature: T = %.2E\n",1/beta)
				fprintf("Total acceptance rate: %f%%\n\n", n/t*100)
			end
		end

	end

	% Once we finish, we can compute the classical correlations
	% contained in the candidate closest pseudo-classical state
	C = qclab.mutinfo(chi_best,A,B); % More reliable for diagonal rho
	c = qclab.rentropy(chi_best,kron(A,B));
	if abs(C-c)>1E-12
		warning("The closest pseudo-classical state might be pathological.");
		fprintf("χ* = \n\n")
		disp(chi_best)
		qclab.math.is_rdm(chi_best);
	end

	% Final sanity check and return
    Q = q_best;
	assert(abs(imag(Q))<1E-12,"Imaginary discord!"); Q = real(Q);
	assert(abs(imag(C))<1E-12,"Imaginary classical correlations!"); C = real(C);

end