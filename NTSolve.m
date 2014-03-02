function [ X, Z, mu, X_sens, X_sens_tau, Z_sens, Z_sens_tau, mu_sens, mu_sens_tau, iter ] = NTSolve(Q, C, lambda, A, a, tau, tol, X, Z , mu, Participation, display)

% Solve the NT system up to tolerance tol

    %Some dimensions
    Nvar   = size(X,1);
    Nconst = size(A,3);
    n = Nvar*(Nvar+1)/2;
    
    %Form NT system & compute residual
    [ rhs, NTMat, norm_residual ] = FormNTSystem(  X, Z, mu, lambda, Q, A, a, C, Participation, tau  );

    %Storage
    Record.res_norm = [norm_residual];
    Record.alpha    = [];
        
    iter = 0;
    while (norm_residual > tol) || (iter == 0)

        % Compute NT direction
        NTstep = NTMat\rhs;

        dmu = NTstep(1:Nconst);
        dX  = smat(NTstep(Nconst+1:Nconst+n));
        dZ  = smat(NTstep(Nconst+n+1:end));

        
        alphaX = alphaPos( X, dX );
        alphaZ = alphaPos( Z, dZ );
         

        

        X  = X  + alphaX*dX;
        Z  = Z  + alphaZ*dZ;
        mu = mu + alphaZ*dmu;
      
        %Update residual & NT system
        NTMatOLD = NTMat; %Keep previous NTMat for sensitivity computation (in case factorization is not updated)
        [ rhs, NTMat, norm_residual ] = FormNTSystem(  X, Z, mu, lambda, Q, A, a, C, Participation, tau  );

        %Storage
        Record.res_norm = [Record.res_norm;norm_residual];
        %Record.alpha    = [Record.alpha;alpha];
        
        iter = iter + 1;
        
        %Display
        if (display == 2)
            figure(1);clf
            subplot(1,2,1)
            semilogy(Record.res_norm,'linestyle','none','marker','.')
            line([1 length(Record.res_norm)],[tol tol],'color','r')
            title('Tol');grid on;axis tight
            subplot(1,2,2)
            plot(Record.alpha,'linestyle','none','marker','.')
            title('Step-size');grid on
        end
    

    end

    
    PosCheck( X, Z );
    
    %Compute Gradient of X w.r.t lambda
    for k = 1:length(Participation)
        rhs_sens(:,k) = [zeros(size(a));
                         svec(C(:,:,Participation(k)));
                         zeros(n,1)];
    end

    rhs_sens(:,k+1) = [zeros(size(a));
                       zeros(n,1);
                       svec(eye(Nvar))];
   
    Sol_sens = NTMatOLD\rhs_sens;

    %Extract Sensitivities
    for k = 1:length(Participation)
        mu_sens(:,Participation(k))  = Sol_sens(1:Nconst,k);
        X_sens(:,:,Participation(k)) = smat(Sol_sens(Nconst+1:Nconst+n,k));
        Z_sens(:,:,Participation(k)) = smat(Sol_sens(Nconst+n+1:end,k));
    end

    mu_sens_tau = Sol_sens(1:Nconst,k+1);
    X_sens_tau = smat(Sol_sens(Nconst+1:Nconst+n,k+1));
    Z_sens_tau = smat(Sol_sens(Nconst+n+1:end,k+1));
    
    
    
    if (display == 1)
        figure(1);clf
        subplot(1,2,1)
        semilogy(Record.res_norm,'linestyle','none','marker','.')
        line([1 length(Record.res_norm)],[tol tol],'color','r')
        title('Tol');grid on;axis tight
        subplot(1,2,2)
        plot(Record.alpha,'linestyle','none','marker','.')
        title('Step-size');grid on
    end
    

end


