function [ rhs, norm_residual ] = FormResidualLinear( X, Z, mu, lambda, Q, A, a, C, Participation, tau, dX, dZ, dmu, dlambda, dtau )
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

    rp = a - AMat*svec(X+dX);

    %Dualize constraints
    DC = lambdaC( lambda+dlambda, C, Participation );

    Rd = Q + DC - (Z+dZ);
    for k = 1:Nconst
        Rd = Rd - (mu(k)+dmu(k))*A(:,:,k);
    end

    Rc = (tau+dtau)*eye(Nvar) - HP(P,IP,X*Z + X*dZ + dX*Z); 

    rc = svec(Rc);
    rd = svec(Rd);

    rhs = [rp;
           rd;
           rc];


    norm_residual = norm(rhs);

end

