clear all
close all
clc

% run /Users/sebastien/Desktop/cvx/cvx_setup

% Test on 2 agents

%First or Second order Dual update
Order = 'First';
Lipschitz = 10;

Nvar   = 8;
Nconst = 2;
Nagent = 2;
Nlink  = 1;

WeightNuclear = 10;

for agent = 1:Nagent
    Q(:,:,agent) = random('norm',0,1,Nvar,Nvar);
    Q(:,:,agent) = Q(:,:,agent) + Q(:,:,agent).';
    while (min(real(eig(Q(:,:,agent)))) < 0)
        Q(:,:,agent) = Q(:,:,agent) + (1-min(real(eig(Q(:,:,agent)))))*eye(Nvar);
    end
    
    C(:,:,agent) = random('norm',0,1,Nvar,Nvar);
    C(:,:,agent) = C(:,:,agent) + C(:,:,agent).';

    for k = 1:Nconst
        A(:,:,k,agent) = random('norm',0,1,Nvar,Nvar);
        A(:,:,k,agent) = A(:,:,k,agent) + A(:,:,k,agent).';
        a(k,agent)     = random('norm',0,1);
    end
    
end

lambda = 0;

for agent = 1:Nagent
    cvx_begin

        cvx_precision( 1e-1 );

        variable Xk(Nvar,Nvar);% symmetric;

        dual variables y{1+Nconst}  

        minimize( trace((Q(:,:,agent) + WeightNuclear*eye(Nvar) + lambda*C(:,:,agent))*Xk) );
        subject to
            Xk == semidefinite(Nvar) : y{1};
            for k = 1:Nconst
                trace(A(:,:,k,agent)*Xk) == a(k,agent) : y{k+1} ;
            end

        cvx_problem
    cvx_end

    Z0(:,:,agent) = y{1};
    for k = 1:Nconst
        mu(k,agent) = y{k+1};
    end
    X0(:,:,agent) = Xk;
end

tol = 1e-6;
tau_table = logspace(1,-3,20);

X = X0; Z = Z0;

for tau_index = 1:length(tau_table)
    tau = tau_table(tau_index);

    iter = 1;


    res = 1;DualIncrease = 1;
    while (abs(res) > 1e-3) && (DualIncrease == 1)
        D = 0;res = 0;
        for agent = 1:Nagent
            [ X(:,:,agent), Z(:,:,agent), mu(:,agent), X_sens(:,:,agent) ] = NTSolve(Q(:,:,agent) + WeightNuclear*eye(Nvar), C(:,:,agent), lambda, A(:,:,:,agent), a(:,agent), tau, tol, X(:,:,agent), Z(:,:,agent) , mu(:,agent));

            %Dual value
            D = D + trace((Q(:,:,agent)  + WeightNuclear*eye(Nvar)  + lambda*C(:,:,agent))*X(:,:,agent));

            %Residual
            res = res + trace(C(:,:,agent)*X(:,:,agent));
        end

        %Compute dual Hessian
        DH = 0;
        for agent = 1:Nagent
            DH = DH + trace(C(:,:,agent)*X_sens(:,:,agent));
        end
        
        %Store
        lambda_store{tau_index}(iter) = lambda;
        D_store{tau_index}(iter)      = D;
        res_store{tau_index}(iter)    = norm(res);
        DH_store{tau_index}(iter)     = norm(DH);
        
        %%%%%%%%%%%%%%%%
        %  Dual Update %
        %%%%%%%%%%%%%%%%
        
        switch Order
            case 'First'
                lambda = lambda + res/Lipschitz;
            case 'Second'     
                lambda = lambda - DH\res;
        end

%         if iter > 1
%             %Check Dual Increase
%             if (D_store{tau_index}(iter) < D_store{tau_index}(iter-1))
%                 DualIncrease = 0;
%             end
%         end
        
        %Next iterate
        iter = iter + 1;

%         figure(2);clf
%         subplot(1,2,1)
%         plot(lambda_store{tau_index},D_store{tau_index},'linestyle','none','marker','.');hold on
%         plot(lambda_store{tau_index}(end),D_store{tau_index}(end),'linestyle','none','marker','.','color','r');
%         title('Dual Value')
%         subplot(1,2,2)
%         semilogy(res_store{tau_index});hold on
%         title('Residual')
%         grid on
    end
    iter_store(tau_index) = iter;
%     pause
    Central_path(tau_index,:) = [lambda_store{tau_index}(end) D_store{tau_index}(end)];

    figure(3);
    plot(lambda_store{tau_index},D_store{tau_index},'linestyle',':','marker','.');hold on
    plot(lambda_store{tau_index}(end),D_store{tau_index}(end),'linestyle','none','marker','.','color','r');
    [Dmax,IndexDmax] = max(D_store{tau_index});
    plot(lambda_store{tau_index}(IndexDmax),Dmax,'linestyle','none','marker','o','color','c');
    plot(Central_path(:,1),Central_path(:,2))
    title('Central path & with Dual values over local updates')
    xlabel('Dual variable');ylabel('Dual Value')
    
    figure(4)
    subplot(2,1,1)
    plot(DH_store{tau_index});hold on;grid on
    title('Dual curvature at each barrier value')
    xlabel('Dual iteration');ylabel('||Dual Hessian||')
    subplot(2,1,2)
    bar(-log10(tau_table(1:tau_index)),iter_store)
    xlabel('-log10(Barrier)');ylabel('#Dual iter')
    
    title('# of iteration at each barrier value')
    
end