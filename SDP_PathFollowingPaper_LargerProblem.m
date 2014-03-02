clear all
close all
clc

% run /Users/sebastien/Desktop/cvx/cvx_setup

%List_of_Methods = {'Second_All_Predictor','Adaptive'};


Method = 'Second_Dual_Predictor';
Marker = {'o','*','x'};


label = '_Large2';

saveFigs = 1;
display_subproblems = 0;

Nvar   = 20;
Nconst = 5;

Nproblem = 100; %Number of problems to be solved

%works but a bit slow...
ratio_max = 10;
gamma_min = 0.25;%0.4;

gamma_damping = 0.0;

p_adaptive = [ratio_max 1;1 1]\[1;gamma_min];

SDPstd = 1e1;

for problem_index = 1:Nproblem
    cvx_status = 'Infeasible';
    while not(strcmp(cvx_status,'Solved'))
        Data = 100; %Number of local problems

        if isstr(Data)
            %Load predefined scenario
            load(Data)
        else
            %Create random scenario
            Nagent = Data;

            %Define coupling constraints (C-set)
            switch Data

                case 4 %Four agents
                    Cset = {[1 2],
                            [2 3],
                            [3 4]};

                case 3 %Three agents
                    Cset = {[1 2],
                            [2 3]};

                case 2 %Two agents
                    Cset = {[1 2]};
            end

            if (Data > 4)

                for k = 1:Nagent-1
                    Cset{k} = [k  k+1];
                end
            end

            %Participating set
            for k = 1:Nagent
                P{k} = []; %create set of empty sets
            end
            for const = 1:length(Cset)
                for k = 1:length(Cset{const})
                    P{Cset{const}(k)} = [P{Cset{const}(k)} const];
                end
            end

            WeightNuclear = 10;

            for agent = 1:Nagent
                Q{agent} = random('norm',0,SDPstd,Nvar,Nvar);
                Q{agent} = Q{agent} + Q{agent}.';
                while (min(real(eig(Q{agent}))) < 0)
                    Q{agent} = Q{agent} + (1-min(real(eig(Q{agent}))))*eye(Nvar);
                end


                % C{agent}(:,:,index of coupling constraint)
                for const = 1:length(P{agent})
                    C{agent}(:,:,P{agent}(const)) = random('norm',0,SDPstd,Nvar,Nvar);
                    C{agent}(:,:,P{agent}(const)) = C{agent}(:,:,P{agent}(const)) + C{agent}(:,:,P{agent}(const)).';
                end


                for k = 1:Nconst
                    A{agent}(:,:,k) = random('norm',0,SDPstd,Nvar,Nvar);
                    A{agent}(:,:,k) = A{agent}(:,:,k) + A{agent}(:,:,k).';        
                    a{agent}(k,1)   = random('norm',0,SDPstd);
                end

            end
        end
        Ncouple = length(Cset);

%         %% Check feasibility using cplex %%
%         cvx_begin
%             cvx_precision( 1e-6 );
% 
%             dual variables lambda{length(Cset)};
% 
%             Cost = 0;
%             for agent = 1:Nagent
%                 %Matrices
%                 eval(['variable X',num2str(agent),'(Nvar,Nvar);']);
% 
%                 %Local dual variables
%                 eval(['dual variables mu',num2str(agent),'{Nconst} ;']);
%                 eval(['dual variables Z',num2str(agent),' ;']);
% 
%                 eval(['Cost = Cost + trace( (Q{agent}+ WeightNuclear*eye(Nvar))*X',num2str(agent),'  ) ;'])
%                 subject to
%                     for k = 1:Nconst
%                         eval(['trace(A{agent}(:,:,k)*X',num2str(agent),') == a{agent}(k) : mu',num2str(agent),'{k}  ;'])
%                     end
% 
%                     eval(['X',num2str(agent),' == semidefinite(Nvar): Z',num2str(agent),';'])
%             end
% 
% 
%             for coupling = 1:length(Cset)
%                 eval(['Const',num2str(coupling),' = 0;'])
%                 for k = 1:length(Cset{coupling})
%                     agent = Cset{coupling}(k);
%                     eval(['Const',num2str(coupling),' = Const',num2str(coupling),' + trace(C{agent}(:,:,coupling)*X',num2str(agent),');'])
%                 end
%                 eval(['Const',num2str(coupling),' == 0 : lambda{coupling};'])
%             end
% 
% 
%             minimize(Cost);
% 
%             cvx_problem
%         cvx_end
        cvx_status = 'Solved';
    end


%     for agent = 1:Nagent
%         eval(['X_cvx{agent} = X',num2str(agent),';']);
%         for k = 1:Nconst
%             eval(['mu_cvx{agent}(k) = mu',num2str(agent),'{k};']);
%         end
%         eval(['Z_cvx{agent} = Z',num2str(agent),';']);
%     end
%     for k = 1:length(Cset)
%         eval(['lambda_cvx(k) = lambda{k};']);
%     end

    %%%%%   Check solution   %%%%%%
    %
    % for agent = 1:Nagent
    %     %Coupling Residual 
    %     for const = 1:length(Cset) %walk through the C-sets
    %         resC{agent}(const,1) = 0; %res for constraint 'const'
    %         for index = 1:length(Cset{const}) %sum over all the agents participating in constraint 'const'
    %             participating_agent = Cset{const}(index);
    %             resC{agent}(const,1) = resC{agent}(const,1) + trace(C{participating_agent}(:,:,const)*X_cvx{participating_agent});
    %             %display(['Agent = ',num2str(participating_agent),' / Const = ',num2str(const)]);
    %         end
    %     end
    % 
    %     %Primal residual
    %     for k = 1:Nconst
    %         resP{agent}(k) = trace(A{agent}(:,:,k)*X_cvx{agent}) - a{agent}(k)  ;
    %     end        
    %     
    %     %Dual residual 
    %     DC = lambdaC( lambda_cvx, C{agent}, P{agent} );
    % 
    %     resD{agent} = Q{agent} + WeightNuclear*eye(Nvar) - DC - Z_cvx{agent};
    %     for k = 1:Nconst
    %         resD{agent} = resD{agent} - mu_cvx{agent}(k)*A{agent}(:,:,k);
    %     end
    % 
    % end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    save DataLarge
    
    %%%% Start Distributed Algorithm %%%%
    display('Algorithm 2')

        lambda    = zeros(length(Cset),1);
        xFGM      = lambda; %for Fast Gradient Method
        xFGM_prev = xFGM;

        toltau    = 1e-6;

        %Naive guess
        for agent = 1:Nagent
            X{agent} = eye(Nvar); Z{agent} = eye(Nvar); 
            mu{agent} = zeros(Nconst,1);
        end    



        gamma = gamma_min;

        dlambda_dtau = 0;
        tau = 1;
        res = 1;
        tau_index = 1;

        tol = tau/2;%sqrt(toltau);
                
        %%%% Initialize on the central path (step 1)
        display(['Initialization phase (computable in parallel, but performed serially here)'])
        
        [ X, Z, mu, lambda, iter(problem_index) ] = Center( Q, WeightNuclear, C, A, a, X, Z, mu, lambda, P, Cset, tau, tol );

%         display('-------------------------------------------')
%         Q{1}
%         WeightNuclear*eye(Nvar)
%         Q{1}+WeightNuclear*eye(Nvar)
%         return
%         for agent = 1:Nagent
%             [ rhs{agent}, NTMat{agent}, norm_residual(agent) ] = FormNTSystem(  X{agent}, Z{agent}, mu{agent}, lambda, Q{agent}+WeightNuclear*eye(Nvar), A{agent}, a{agent}, C{agent}, P{agent}, tau  );
%             
%         end
%         norm_residual
%         return
%         for agent = 1:Nagent
%             [X2{agent}, Z2{agent}, mu2{agent}] = NTStep(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau,  X{agent}, Z{agent}, mu{agent}, P{agent});
%         end
        
        
        %return
        display(['Dual iterations']);h = waitbar(0,'completed');
        while tau > toltau || norm(res) > toltau

            Store{problem_index}.tau(tau_index) = tau;

            %%%% Step NT system (step 2)

            display('Local NT steps (computable in parallel, but performed serially here)')

            for agent = 1:Nagent
                [a1{agent}, a2{agent}, a3{agent}, a4{agent}, a5{agent}, a6{agent}, a7{agent}, a8{agent}, a9{agent}] = NTStep(Q{agent} + WeightNuclear*eye(Nvar), C{agent}, lambda, A{agent}, a{agent}, tau,  X{agent}, Z{agent}, mu{agent}, P{agent});
                waitbar(agent/Nagent)
            end

            
            %for using parfor...
            X      = a1;
            Z      = a2;
            mu     = a3;
            X_sens = a4;
            X_sens_tau = a5;
            Z_sens = a6;
            Z_sens_tau = a7;
            mu_sens = a8;
            mu_sens_tau = a9;
            
            %Initiate res and D for summation
            res = 0;D = 0;
            for agent = 1:Nagent
                %Dual value
                DC = lambdaC( lambda, C{agent}, P{agent} );
                D = D + trace((Q{agent}  + WeightNuclear*eye(Nvar)  + DC)*X{agent}) - tau*log(det(X{agent}));

                %Dual Residual 
                for const = 1:length(Cset) %walk through the C-sets
                    res(const,1) = 0; %res for constraint 'const'
                    for index = 1:length(Cset{const}) %sum over all the agents participating in constraint 'const'
                        participating_agent = Cset{const}(index);
                        res(const,1) = res(const,1) + trace(C{participating_agent}(:,:,const)*X{participating_agent});
                    end
                end
            end


            Store{problem_index}.res(tau_index) = norm(res);
            Store{problem_index}.lambda(:,tau_index) = lambda;

            %%%% Dual gradient and Hessian (step 3)

            %Compute dual Hessian
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
            
            %Delta function 
            delta_func{problem_index}(tau_index) = sqrt(-res.'*(DH\res));%/sqrt(tau);
            
            %Compute gradient of D w.r.t barrier tau
            Dtau = zeros(length(Cset),1);
            for l = 1:length(Cset)
                for i = 1:length(Cset{l})
                    k = Cset{l}(i);                     
                    Dtau(l) = Dtau(l) + trace(C{k}(:,:,l)*X_sens_tau{k});                            
                end
            end

            %%%% Update barrier tau (step 4)          


            %Update tau
            Store{problem_index}.gamma(tau_index) = gamma;

            %Adaptive heuristic
            %gamma_new = max(gamma_min,min(1,p_adaptive(1)*(delta_func{problem_index}(tau_index)/tau) + p_adaptive(2)));
            
            if (tau < 1)
                p(1) = (1-gamma_min)/(tau^(-0.5) - 1);p(2) = gamma_min - p(1);
                gamma_new = p(1)*delta_func{problem_index}(tau_index)/tau + p(2); 
                gamma_new = max(gamma_min,min(1,p_adaptive(1)*(delta_func{problem_index}(tau_index)/tau) + p_adaptive(2))); 
            else
                gamma_new = gamma;
           
            end
            
            
            
            %Damping of the gamma reduction
            if gamma_new < gamma
                gamma_new = (1-gamma_damping)*gamma_new + gamma_damping*gamma;
            end
            gamma = gamma_new;
            
            gamma = 0.5;
            
            if (tau <= toltau)
                gamma = 1;
            end
            
            display(['Iter: ',num2str(tau_index),'  | Dual Residual: ',num2str(norm(res)),'  | tau: ',num2str(tau),'  | d: ',num2str(delta_func{problem_index}(tau_index)),'  | gamma: ',num2str(gamma)])
            
            
            tau_new = gamma*tau;
            


            dtau = tau_new - tau;
            tau = tau_new;



            %%%%  Dual Update  (step 5)  

            switch Method
                case 'Second_Classic_Predictor'
                    dlambda = - DH\res;
                    %dlambda = 0*dlambda;
                case 'Second_Dual_Predictor'
                    dlambda = - DH\(res + Dtau*dtau);

            end

            


            %Update dual variables with predictor
            lambda = lambda + dlambda;

            %%%% Predictor on the local solutions (step 6)
            switch Method

                case 'Second_Classic_Predictor'
                     for agent = 1:Nagent
                        dmu = mu_sens_tau{agent}*dtau;
                        dX  = X_sens_tau{agent}*dtau;
                        dZ  = Z_sens_tau{agent}*dtau;

                        PosCheck( X{agent}, Z{agent} );

                    end

                case 'Second_Dual_Predictor'

                    for agent = 1:Nagent
                        dmu = mu_sens_tau{agent}*dtau;
                        dX  = X_sens_tau{agent}*dtau;
                        dZ  = Z_sens_tau{agent}*dtau;

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

            tau_index = tau_index + 1;


        end
        Store{problem_index}.iter = tau_index-1;

        close(h)
    

    %%%%%   Check solution   %%%%%%
    % for agent = 1:Nagent
    %     %Coupling Residual 
    %     for const = 1:length(Cset) %walk through the C-sets
    %         resC2{agent}(const,1) = 0; %res for constraint 'const'
    %         for index = 1:length(Cset{const}) %sum over all the agents participating in constraint 'const'
    %             participating_agent = Cset{const}(index);
    %             resC2{agent}(const,1) = resC2{agent}(const,1) + trace(C{participating_agent}(:,:,const)*X{participating_agent});
    %             %display(['Agent = ',num2str(participating_agent),' / Const = ',num2str(const)]);
    %         end
    %     end
    % 
    %     %Primal residual
    %     for k = 1:Nconst
    %         resP2{agent} = trace(A{agent}(:,:,k)*X{agent}) - a{agent}(k)  ;
    %     end        
    % 
    % 
    %     %Dual residual 
    %     DC = lambdaC( lambda, C{agent}, P{agent} );
    % 
    %     resD2{agent} = Q{agent} + WeightNuclear*eye(Nvar) + DC - Z{agent};
    %     for k = 1:Nconst
    %         resD2{agent} = resD2{agent} - mu{agent}(k)*A{agent}(:,:,k);
    %     end
    % end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

LS = '-';FS = 16;



%%%% Compare all methods, version 2
fig = figure(1);clf

for problem_index = 1:Nproblem
    semilogy(Store{problem_index}.res,'linestyle',LS,'marker',Marker{1},'color','k');hold on
    semilogy(Store{problem_index}.tau,'linestyle','-','color','k');hold on
   % semilogy(delta_func{problem_index},'linestyle','--','marker',Marker{1},'color','r');hold on
    
    if (problem_index == 1)
        legend({'Algorithm 2\qquad','$$\tau$$'},'interpreter','latex','fontsize',FS,'location','southwest')
    end
    line([1 Store{problem_index}.iter],[toltau toltau],'linestyle','--','color','k')

end



ylabel('$$\|\nabla_\lambda D\|$$ and $$\tau$$','interpreter','latex','fontsize',FS)

%title(List_titles{problem_index})
grid on
axis tight

xlabel('Iteration','fontsize',FS)


title('NDSDP Path-Following, barrier and dual residual','fontsize',FS)
set(gca,'fontsize',FS);



if saveFigs
    PapPos = get(gcf,'PaperPosition');%PapPos(4) = 1.5*PapPos(4);
    set(gcf,'PaperPosition',PapPos)

    FileName = ['/Users/sebastien/Desktop/OPTICON/Publications/CDC2014/GP/SDP_Decomposition/Figures/PathFollowing_tau_and_res',label];
    exportfig(fig, FileName,'color','cmyk')
end

%save Results_Large1


figure(10);clf 
for problem_index = 1:length(delta_func)
    semilogy(delta_func{problem_index}./Store{problem_index}.tau(1:length(delta_func{problem_index})));hold on
end
grid on
