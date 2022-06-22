function [newtau,tauopt,scalar] = bisectingH(delta,taubar,msqe0,T,M,VT,N0,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,maxIt)
bisector = 0;
disp('Starting the bisecting houristic to current optimum')
for j = 1:maxIt
    try
        scalar = 2^(-bisector);    
        newtau = taubar + delta*scalar;
        tauopt = map(T,M,VT,N0,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,newtau);
        msqe = sum((newtau-tauopt).^2,'all')/(length(T)*length(M));
        if msqe <= msqe0
            disp(['Best response underrelaxed by factor ' num2str(2^bisector) ' was a success!'])
            return
        end
        disp(['Best response underrelaxed by factor ' num2str(2^bisector) ' had a higher error of ' num2str(msqe)])
    catch
        disp(['Best response underrelaxed by factor ' num2str(2^bisector) ' resulted in an error'])
    end
    bisector = bisector + 1;
    disp(['Now trying to underrelax by factor ' num2str(2^bisector)])
end
disp('bisecting houristic ran out of iterations')