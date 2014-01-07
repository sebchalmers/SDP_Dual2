function [ X, Z, mu, X_sens ] = NTSolve(Q, C, lambda, A, a, tau, tol, X, Z , mu, Participation, display)
    % Solve the NT system up to tolerance tol


    
    Nvar   = size(X,1);
    Nconst = size(A,3);
    
    tolIP = 1;

    Record.tau   = [];
    Record.tolIP = [];
    Record.alpha = [];


    while (tolIP > tol)

        %Form NT system
        [ W, G, IG ] = computeWandG( X, Z );
        P = IG;IP = G;

        n = Nvar*(Nvar+1)/2;
        AMat = zeros(Nconst,n);
        for k = 1:Nconst
            AMat(k,:) = svec(A(:,:,k)).';
        end
        I = eye(n);
        E = skron(P,IP.'*Z);
        F = skron(P*X,IP.');

        NTMat = [zeros(Nconst,Nconst)   AMat        zeros(Nconst,n);
                      AMat.'           zeros(n,n)        I         ;
                 zeros(n,Nconst)         E               F          ];

        rp = a - AMat*svec(X);

        %Dualize constraints
        DC = lambdaC( lambda, C, Participation );
        
        Rd = Q + DC - Z;
        for k = 1:Nconst
            Rd = Rd - mu(k)*A(:,:,k);
        end

        Rc = tau*eye(Nvar) - HP(P,IP,X*Z);

        rc = svec(Rc);
        rd = svec(Rd);
        
        rhs = [rp;
               rd;
               rc];


        tolIP = norm(rhs);

        % Compute NT step
        NTstep = NTMat\rhs;

        dmu = NTstep(1:Nconst);
        dX  = smat(NTstep(Nconst+1:Nconst+n));
        dZ  = smat(NTstep(Nconst+n+1:end));

        % Compute update

        % Step-size (maintain positivity)
         
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
                display(['Tiny step-size: alpha < 1e-6']);
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

        Record.tolIP = [Record.tolIP;tolIP];
        Record.alpha = [Record.alpha;alpha];
        Record.tau   = [Record.tau;tau];

    end
    
    %Compute Gradient of X w.r.t lambda
    for k = 1:length(Participation)
        rhs_sens(:,k) = [zeros(size(rp));
                         svec(C(:,:,Participation(k)));
                         zeros(size(rc))];
    end
            
    Sol_sens = NTMat\rhs_sens;

    for k = 1:length(Participation)
        X_sens(:,:,Participation(k)) = smat(Sol_sens(Nconst+1:Nconst+n,k));
    end
    %X*Z
    if (display == 1)
        figure(1);clf
        subplot(1,2,1)
        semilogy(Record.tolIP,'linestyle','none','marker','.')
        line([1 length(Record.tolIP)],[tol tol],'color','r')
        title('Tol');grid on;axis tight
        subplot(1,2,2)
        plot(Record.alpha,'linestyle','none','marker','.')
        title('Step-size');grid on
    end

end

