function [ rhs, NTMat, norm_residual ] = FormNTSystem( X, Z, mu, lambda, Q, A, a, C, Participation, tau )
% Form Nesterov-Todd system (residual and KKT matrix)

    Nvar   = size(X,1);
    Nconst = size(A,3);
    
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


    norm_residual = norm(rhs);

end

