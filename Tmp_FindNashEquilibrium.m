LSD = 1/17; %log-sate-difference: The amount of indices between m = exp(1) and m = exp(2)
M = exp((1:200)*LSD)/exp(LSD); %Creating all indices of the mass states

dt = .03;
T = 1:dt:50; %Creating all indices of the time states

offset = round(log(400)/LSD); %The peak of the preference function (bete)
p = normpdf(-(offset-1)*LSD:LSD:(offset-1)*LSD); %creating the preference function (with alpha^2 = exp(1))
omega = zeros(length(T),1);
sigma = 20/4.5; %sigma of egg hatching
Nomega = normpdf(-4:dt/sigma:4); %the peak will be at 4.5 times sigma
omega(1:length(Nomega)) = Nomega/sum(Nomega); %normalise
rho = 1e5; %clearance rate
eps = .1; %efficiency of conversion
C0 = 2/eps; %coefficient for metablolic loss
Cmax = C0*M.^(3/4); %metabolic loss
K = Cmax*eps/10; % the denominator is phi_0, basic feeding level
mu0 = 0.05*M.^(-.25); %basi mortality rate
nu0 = ones(1,length(M))*6/(400/exp(1)*eps).*M.^(-.25); %background mortality rwa. (b is the first number)
gma0 = ones(1,length(M))*6.*M.^(.75); %background ingestion rwa. (b is the first number)

lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));% least viable feeding level (given only background)
taubar = repmat(lambdab,length(T),1); %taubar at lambdab

VT = ones(1,length(M))./(1+exp(-4*(log(M)-log(1000)))); %Value function at T_end
%%
[taunash,msqe,recordedTaubars,StepScalar] = searchBrouwers(T,M,VT,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar,200,1e-9,[],@bisectingH,8);
%The algorithm that seaches for a fixpoint by underrelaxed best response
%procedure
%The last 5 arguments are repsectively:
%1 maximum allowed main loop iterations
%2 tolerance for stopping
%3 vector of indices of iterations where the computer should save current
%process (slow)
%4 choice of algortihm (don't change)
%maximum number of bisections

%The outputs are:
%tauopt(t,y), the optimal strategy (T x M array)
%mean squared error(#iterations), the mean squared error after each main
%loop iteration
%recordedTaubars, the specified save points (from 3rd last input)
%StepScalar(#iterations): the constant (j) that was used at iteration each #
%%
[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taunash); 
[V,~] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0); %These two needs to be run again to get all results