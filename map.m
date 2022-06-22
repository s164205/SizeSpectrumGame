function tau = map(T,M,VT,N0,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar)

[~,gma,nu] = forwardTransport(T,M,N0,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
[~,tau] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0);