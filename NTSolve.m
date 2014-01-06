function [ X, Z, mu, X_sens ] = NTSolve(Q, C, lambda, A, a, tau, tol, X, Z , mu, Participation)
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

        sig = 1;%trace(X*Z)/Nvar;
        Rc = tau*sig*eye(Nvar) - HP(P,IP,X*Z);

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

        while (cond < 0.01*tolcond) 
            alpha = 0.9*alpha;
            cond = min([eig(X + alpha*dX);
                        eig(Z + alpha*dZ)]);

        end


        X  = X  + alpha*dX;
        Z  = Z  + alpha*dZ;

        mu = mu + alpha*dmu;

        Record.tolIP = [Record.tolIP;tolIP];
        Record.alpha = [Record.alpha;alpha];
        Record.tau   = [Record.tau;tau];

    end
    
    %Compute Gradient of X w.r.t lambda
    for k = 1:size(C,3)
        rhs_sens(:,k) = [zeros(size(rp));
                         svec(C(:,:,k));
                         zeros(size(rc))];
    end
            
    Sol_sens = NTMat\rhs_sens;
    for k = 1:size(Sol_sens,2)
        X_sens(:,:,k) = smat(Sol_sens(Nconst+1:Nconst+n,k));
    end
    %X*Z
% 
%     figure(1);clf
%     subplot(1,2,1)
%     semilogy(Record.tolIP,'linestyle','none','marker','.')
%     line([1 length(Record.tolIP)],[tol tol],'color','r')
%     title('Tol');grid on;axis tight
%     subplot(1,2,2)
%     plot(Record.alpha,'linestyle','none','marker','.')
%     title('Step-size');grid on


end

