function [newtau,tauopt,URconstant] = underRelaxedIBR(delta,taubar,msqe,T,M,VT,N0,p,offset,rho,Cmax,eps,K,lambda,mu0,gma0,nu0,URconstant)

newtau = taubar + delta*URconstant;
tauopt = [];