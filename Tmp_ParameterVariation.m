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

currtau = taubar;
%% 
Example = cell(length(R),1); %Here I search for different clearance rates, as an example
R = exp(0:1:50)-1; %vector of parameter choices

for i = 1:length(R)
    rho = R(i); %setting the current parameter
    tic
    [Example{i},~,~,~] = searchBrouwers(T,M,VT,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar,20,1e-8,[],@bisectingH,8); %finding Nash for the parameter (very slow)
    fprintf(2,'Iteration %.f, with b = %.4f complete!\n It took %.4f seconds\n',i,R(i),toc);
    currtau = Example{i}; %Using current Nash for next initial guess
end
%%
F = zeros(length(R),1); %total Fitness
for i = 1:length(R)
    N = forwardTransport(T,M,omega,p,offset,R(i),Cmax,eps,K,mu0,gma0,nu0,Example{i}); %finding population dynamic given the Nash strategy
    F(i) = N(end,:)*VT'; %total fitness
end