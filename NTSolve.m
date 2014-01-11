function [ X, Z, mu, X_sens, X_sens_tau ] = NTSolve(Q, C, lambda, A, a, tau, tol, X, Z , mu, Participation, display)
    % Solve the NT system up to tolerance tol


    
    Nvar   = size(X,1);
    Nconst = size(A,3);
    n = Nvar*(Nvar+1)/2;
    
    norm_residual = 1e6;

    Record.tau      = [];
    Record.res_norm = [];
    Record.alpha    = [];


    while (norm_residual > tol)

        %Form NT system
        [ rhs, NTMat, norm_residual ] = FormNTSystem(  X, Z, mu, lambda, Q, A, a, C, Participation, tau  );

        % Compute NT direction
        NTstep = NTMat\rhs;

        dmu = NTstep(1:Nconst);
        dX  = smat(NTstep(Nconst+1:Nconst+n));
        dZ  = smat(NTstep(Nconst+n+1:end));

        % Compute update
        alpha = 1;

        tolcond = min([eig(X);
                       eig(Z)]);

        cond = min([eig(X + alpha*dX);
                    eig(Z + alpha*dZ)]);

        while (cond < 0.9*tolcond) 
            alpha = 0.9*alpha;
            cond = min([eig(X + alpha*dX);
                        eig(Z + alpha*dZ)]);
            if (alpha < 1e-6)
                alpha
            end
        end

        %display(['NT Line-search, tolcond: ',num2str(tolcond),' | cond: ',num2str(cond),' | alpha = ',num2str(alpha)])
        
        X  = X  + alpha*dX;
        Z  = Z  + alpha*dZ;

%         if (min(eig(X)) < 1e-9) || (min(eig(Z)) < 1e-9)
%             display(['Eig(X) = ',num2str(eig(X).')])
%             display(['Eig(Z) = ',num2str(eig(Z).')])
%             keyboard
%         end
        mu = mu + alpha*dmu;

        Record.res_norm = [Record.res_norm;norm_residual];
        Record.alpha    = [Record.alpha;alpha];
        Record.tau      = [Record.tau;tau];
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

    %Compute Gradient of X w.r.t lambda
    for k = 1:length(Participation)
        rhs_sens(:,k) = [zeros(size(a));
                         svec(C(:,:,Participation(k)));
                         zeros(n,1)];
    end

    rhs_sens(:,k+1) = [zeros(size(a));
                       zeros(n,1);
                       svec(eye(Nvar))];
    
    Sol_sens = NTMat\rhs_sens;

    for k = 1:length(Participation)
        X_sens(:,:,Participation(k)) = smat(Sol_sens(Nconst+1:Nconst+n,k));
    end

    X_sens_tau = smat(Sol_sens(Nconst+1:Nconst+n,k+1));
    
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


