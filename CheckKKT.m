function [ KKT_res_PD, KKT_res_Barrier ] = CheckKKT( X, Z, mu, lambda, Q, C, A, a, tau, Participation )
%Check KKTs for a local problem
    
    Nvar = size(Q,1);

    %Primal-dual (8)

    Rd_PD = Q - Z + lambdaC(lambda,C,Participation);
    for k = 1:size(A,3)
        Rd_PD = Rd_PD - mu(k)*A(:,:,k);
        Rp_PD(k) = trace(A(:,:,k)*X) - a(k);
    end
    Rc_PD = Z*X - tau*eye(Nvar);


    KKT_res_PD = norm([svec(Rd_PD);svec(Rp_PD);svec(Rc_PD)],inf);


    %Primal Barrier (7)

    Rd_Barrier = Q - tau*inv(X) + lambdaC(lambda,C,Participation);
    for k = 1:size(A,3)
        Rd_Barrier = Rd_Barrier - mu(k)*A(:,:,k);
        Rp_Barrier(k) = trace(A(:,:,k)*X) - a(k);
    end
    


    KKT_res_Barrier = norm([svec(Rd_Barrier);svec(Rp_Barrier)],inf);
 
end

