function [ X, Z, mu, lambda, iter ] = Center( Q, WeightNuclear, C, A, a, X, Z, mu, lambda, P, Cset, tau, tol )
%Center the local solutions and Lagrange multipliers for a fixed, given tau

    %Some dimensions
    Nagent = length(X);
    Nvar   = size(X{1},1);
    Nconst = size(A{1},3);
    n = Nvar*(Nvar+1)/2;

    
    delta_func = tol+eps;iter = 0;
    while delta_func > tol
        iter = iter + 1;
        
        for agent = 1:Nagent    
            [a1{agent}, a2{agent}, a3{agent}, a4{agent}, a5{agent}, a6{agent}, a7{agent}, a8{agent}, a9{agent}, a10(agent)] = NTSolve(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau, delta_func^2, X{agent}, Z{agent}, mu{agent}, P{agent}, 0);                     
        
        end

        %for using parfor...
        X      = a1;
        Z      = a2;
        mu     = a3;
        
        %Check local residuals

        X_sens = a4;
        X_sens_tau = a5;
        Z_sens = a6;
        Z_sens_tau = a7;
        mu_sens = a8;
        mu_sens_tau = a9;

        
        
        %Construct dual residuals
        res = 0; %Initiate res for summation
        for agent = 1:Nagent
            %Dual value
            DC = lambdaC( lambda, C{agent}, P{agent} );
           
            %Dual Residual 
            for const = 1:length(Cset) %walk through the C-sets
                res(const,1) = 0; %res for constraint 'const'
                for index = 1:length(Cset{const}) %sum over all the agents participating in constraint 'const'
                    participating_agent = Cset{const}(index);
                    res(const,1) = res(const,1) + trace(C{participating_agent}(:,:,const)*X{participating_agent});
                end
            end
        end


        %Compute Hessian
        DH = zeros(length(Cset),length(Cset));
        for m = 1:length(Cset)
            for l = 1:length(Cset)
                for i = 1:length(Cset{l})
                    k = Cset{l}(i);
                    if max(m == P{k})                      
                        DH(m,l) = DH(m,l) + trace(C{k}(:,:,l)*X_sens{k}(:,:,m));                            
                    end
                end
            end
        end

        delta_func = sqrt(-res.'*(DH\res));
        display(['d = ',num2str(delta_func),' | |res| = ',num2str(norm(res))]);

        %Lagrange update
        dlambda = - DH\res;

        %Update dual variables with predictor
        lambda = lambda + dlambda;

        %Local predictors
        for agent = 1:Nagent
            dmu = 0;dX  = 0;dZ  = 0;

            for k = 1:length(P{agent})
                dmu = dmu + mu_sens{agent}(:,P{agent}(k))*dlambda(P{agent}(k)) ;
                dX  = dX  + X_sens{agent}(:,:,P{agent}(k))*dlambda(P{agent}(k));
                dZ  = dZ  + Z_sens{agent}(:,:,P{agent}(k))*dlambda(P{agent}(k));
            end

            alphaX = alphaPos( X{agent}, dX );
            alphaZ = alphaPos( Z{agent}, dZ );

            mu{agent} = mu{agent} + dmu;
            X{agent}  = X{agent}  + alphaX*dX;
            Z{agent}  = Z{agent}  + alphaZ*dZ;

            PosCheck( X{agent}, Z{agent} );

        end
    end
    


    
end

