function [N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,tau)

%gma,nu are col. vectors of resp. growth- and deathrates, NOT yet
%multiplied with tau
%gma*tau = g, nu*tau = mu
%gma here is the SPECIFIC ingesstion rwa. contrary to the model showing by
%default the ablosute ingestion rwa.

N = zeros(length(T),length(M));
N(1,1) = omega(1);
nu = zeros(length(T)+1,length(M));
gma = zeros(length(T)+1,length(M));

fprintf('Forward transport: current t = 1 / %.f\n',T(end));

for i = 2:length(T)
    
    gma(i,:) = growthVector(N(i-1,:),M,rho,gma0,tau(i-1,:),p,offset); 
    nu(i,:) = deathVector(N(i-1,:),rho,nu0,tau(i-1,:),p,offset); 
    
   
    gmatau = gma(i,:).*tau(i,:);
    g = (eps*gmatau./(gmatau./Cmax + 1)) - K;
    g = g.*(g>0);
    %starv = (g<0);
    %g = g - g.*starv;
    mu = nu(i,:).*tau(i,:) +mu0;%+ lambda*starv;
    
    gss = g./[diff(M), Inf]; %state-specific growth
    


    G = gmuGenerator(gss,mu);
    
    N(i,:) = (eye(length(M))-G'*(T(i)-T(i-1)))\N(i-1,:)'; %implicit forward step
    
    N(i,1) = N(i,1) + omega(i); %oviposition
    
    
    if mod(i-1,100) == 0 || (i == length(T) && length(T)> 100)
        for j = 1:4+strlength(num2str(T(i-100),'%.f'))+strlength(num2str(T(end),'%.f'))
            fprintf('\b');
        end
        fprintf('%.f / %.f\n',T(i),T(end));
    end
end

gma(end,:) = growthVector(N(end,:),M,rho,gma0,tau(end,:),p,offset);
nu(end,:) = deathVector(N(end,:),rho,nu0,tau(end,:),p,offset);
gma(1,:) = [];
nu(1,:) = [];