function [V,tau] = backwardHJB(TT,M,VT,gma,nu,Cmax,eps,K,mu0)

%gma,nu are col. vectors of resp. growth- and deathrates, NOT yet
%multiplied with tau
%gma*tau = g, nu*tau = mu

T = length(TT);
V = zeros(T,length(VT));
V(end,:) = VT;
tau = zeros(length(TT),length(VT));

fprintf('Backward Hamilton Jacobi Bellman: current t = %.f / %.f\n',TT(end),TT(end));

for i = T-1:-1:1
    
    tau(i,:) = findTypeIIOpt(V(i+1,:),gma(i+1,:),nu(i+1,:),M,Cmax,eps,K);
    
    gmatau = gma(i,:).*tau(i,:);
    g = (eps*gmatau./(gmatau./Cmax + 1)) - K;
    g = g.*(g>0); %necessary to ensure real positives numbers (sometimes -1<<g<0)
    %starv = (g<0); remember that starving should now be impossible
    %g = g - g.*starv;
    mu = nu(i,:).*tau(i,:) +mu0; %+ lambda*starv;
    
    gss = g./[diff(M), Inf]; %state-specific growth
    
    
    G = gmuGenerator(gss,mu);
    
    V(i,:) = (eye(length(VT))-G*(TT(i+1)-TT(i)))\V(i+1,:)';
    
    if i<T-100 && mod(i-1,100) == 0
        for j = 1:4+strlength(num2str(TT(i+100),'%.f'))+strlength(num2str(TT(end),'%.f'))
            fprintf('\b');
        end
        fprintf('%.f / %.f\n',TT(i),TT(end));
    end
end
tau(end,:) = findTypeIIOpt(V(end,:),gma(end,:),nu(end,:),M,Cmax,eps,K);