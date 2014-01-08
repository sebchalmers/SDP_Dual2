clear all
close all
clc

List_of_Methods = {'FGM','Second'};
%List_of_Methods = {'Second'};
Marker = containers.Map;
MarkerList = {'o','x'};

for k = 1:length(List_of_Methods)
    Marker(List_of_Methods{k}) = MarkerList{k};
end

% run /Users/sebastien/Desktop/cvx/cvx_setup

% Test on 2 agents

display_subproblems = 0;

Nvar   = 8;
Nconst = 2;

load Data1
% % Four agents
% Nagent = 4;
% 
% %Define coupling constraints (C-set)
% Cset = {[1 2],
%         [2 3],
%         [3 4]};

% % Three agents
% Nagent = 3;
% 
% %Define coupling constraints (C-set)
% Cset = {[1 2],
%         [2 3]};
% 
% % res = [trace(C{1}(:,:,1)*X{1} + C{2}(:,:,1)*X{2});
% %        trace(C{2}(:,:,2)*X{1} + C{3}(:,:,2)*X{2})]
% 
% % %% Two agents
% % Nagent = 2;
% % 
% % %Define coupling constraints (C-set)
% % Cset = {[1 2]};
% 
% %% res = trace(C{1}*X{1} + C{2}*X{2})
% 
% 
% 
% %Participating set
% for k = 1:Nagent
%     P{k} = []; %create set of empty sets
% end
% for const = 1:length(Cset)
%     for k = 1:length(Cset{const})
%         P{Cset{const}(k)} = [P{Cset{const}(k)} const];
%     end
% end

% WeightNuclear = 10;
% 
% for agent = 1:Nagent
%     Q{agent} = random('norm',0,1,Nvar,Nvar);
%     Q{agent} = Q{agent} + Q{agent}.';
%     while (min(real(eig(Q{agent}))) < 0)
%         Q{agent} = Q{agent} + (1-min(real(eig(Q{agent}))))*eye(Nvar);
%     end
%     
% 
%     % C{agent}(:,:,index of coupling constraint)
%     for const = 1:length(P{agent})
%         C{agent}(:,:,P{agent}(const)) = random('norm',0,1,Nvar,Nvar);
%         C{agent}(:,:,P{agent}(const)) = C{agent}(:,:,P{agent}(const)) + C{agent}(:,:,P{agent}(const)).';
% %         while (min(real(eig(C{agent}(:,:,P{agent}(const))))) < 0)
% %             C{agent}(:,:,P{agent}(const)) = C{agent}(:,:,P{agent}(const)) + (1-min(real(eig(C{agent}(:,:,P{agent}(const))))))*eye(Nvar);
% %         end
%     end
% 
% 
%     for k = 1:Nconst
%         A{agent}(:,:,k) = random('norm',0,1,Nvar,Nvar);
%         A{agent}(:,:,k) = A{agent}(:,:,k) + A{agent}(:,:,k).';        
%         a{agent}(k,1)   = random('norm',0,1);
%     end
% 
% end

for method_number = 1:length(List_of_Methods)
    Method = List_of_Methods{method_number}

    lambda    = zeros(length(Cset),1);
    xFGM      = lambda; %for Fast Gradient Method
    xFGM_prev = xFGM;


    tauEnd = -6; %(order)
    tau_table = logspace(0,tauEnd,40);

    %Naive guess
    for agent = 1:Nagent
        X{agent} = eye(Nvar); Z{agent} = eye(Nvar); 
        mu{agent} = zeros(Nconst,1);
    end

    %initiate storage
    all_res_store = [];
    res_store     = [];
    DH_store      = [];
    D_store       = [];
    lambda_store  = [];        
            
    iter_total = 1; %total number of dual iteration of the algorithm
    tau_vs_iter_total = [];
    for tau_index = 1:length(tau_table)
        tau     = tau_table(tau_index);
        tolDual = tau;
        tolIP   = 1e-6;

        tau_vs_iter_total = [tau_vs_iter_total;
                             tau     iter_total-1];
        
        iter = 1;   %number of dual iterations for the current value of the barrier parameter


        res = 1e6;
        while (norm(res) > tolDual) 
            D = 0;res = 0;
            for agent = 1:Nagent
                [ X{agent}, Z{agent}, mu{agent}, X_sens{agent} ] = NTSolve(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau, tolIP, X{agent}, Z{agent}, mu{agent}, P{agent}, display_subproblems);

                %Dual value
                DC = lambdaC( lambda, C{agent}, P{agent} );
                D = D + trace((Q{agent}  + WeightNuclear*eye(Nvar)  + DC)*X{agent}) - tau*log(det(X{agent}));

                %Residual (to be checked)
                for const = 1:length(Cset) %walk through the C-sets
                    res(const,1) = 0; %res for constraint 'const'
                    for index = 1:length(Cset{const}) %sum over all the agents participating in constraint 'const'
                        participating_agent = Cset{const}(index);
                        res(const,1) = res(const,1) + trace(C{participating_agent}(:,:,const)*X{participating_agent});
                        %display(['Agent = ',num2str(participating_agent),' / Const = ',num2str(const)]);
                    end
                end



            end

            %Compute dual Hessian
            DH = zeros(length(Cset),length(Cset));
            for m = 1:length(Cset)
                for l = 1:length(Cset)
                    for i = 1:length(Cset{l})
                        k = Cset{l}(i);
                        if max(m == P{k})                      
                            DH(m,l) = DH(m,l) + trace(C{k}(:,:,l)*X_sens{k}(:,:,m))   ;                     
                        end
                    end
                end
            end

            Lipschitz = max(abs(eig(DH)));

            %Store
            lambda_store{tau_index}(iter,:) = lambda;
            D_store{tau_index}(iter)        = D;
            res_store{tau_index}(iter)      = norm(res);
            all_res_store(iter_total)       = norm(res);
            DH_store{tau_index}(iter)       = Lipschitz;

            %%%%%%%%%%%%%%%%
            %  Dual Update %
            %%%%%%%%%%%%%%%%

            switch Method
                case 'GA' %Gradient ascent
                    lambda      = lambda + res/Lipschitz;
                    
                case 'FGM' %Fast Gradient with restart
                    restart = 0;
                    if (iter > 1)
                        if (res_store{tau_index}(iter) > res_store{tau_index}(iter-1))
                            display('restart')
                            restart   = 1;
                            xFGM_prev = xFGM; % restart inertia
                            lambda    = xFGM; % update lambda according to gradient ascent on the previous lambda
                        else                            
                            xFGM    = lambda + res/Lipschitz;
                        end
                    else
                        xFGM        = lambda + res/Lipschitz;
                    end
                    
                    if (restart == 0)
                        lambda    = xFGM + (iter-1)/(iter+2)*(xFGM - xFGM_prev);
                        xFGM_prev = xFGM;
                    end
                    
                case 'Second'   %Second order update   
                    lambda    = lambda - DH\res;
            end

            %Next iterate
            iter       = iter + 1;
            iter_total = iter_total + 1;
            
%             figure(1);clf
%             semilogy(res_store{tau_index});hold on
%             title('Residual')
%             grid on
%             ylim([tolDual max(10*tolDual,res_store{tau_index}(1))])

        end
        iter_store(tau_index,method_number) = iter-1;
        
        tau_vs_iter_total = [tau_vs_iter_total;
                             tau     iter_total-1];
        figure(1);clf
        semilogy(res_store{tau_index});hold on
        title('Residual')
        grid on
        ylim([tolDual max(10*tolDual,res_store{tau_index}(1))])
        %pause
            
    end


    figure(2);
    semilogy(all_res_store,'linestyle','none','marker',Marker(Method),'color','k');hold on
    plot(tau_vs_iter_total(1:end,2),tau_vs_iter_total(1:end,1),'linestyle','-','color','k','linewidth',2);hold on
    title('||Dual residual||')
    grid on
    
    figure(3)
    barh(log10(tau_table(1:tau_index)),iter_store)
    ylabel('log10(Barrier)');
    grid on    
    title('# of iteration at each barrier value')
    sum(iter_store(:,1))
    return
end