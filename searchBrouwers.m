function [tau,msqe,recordedTaubars,stepScalar] = searchBrouwers(T,M,VT,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar,maxIt,tol,savePoints,steppingFun,varargin)
try
    msqe = zeros(maxIt,1);
    stepScalar = zeros(maxIt,1);
    recordedTaubars = [];
    tauopt = [];
    for i = 1:maxIt
        if i == 1 || msqe(i-1) > tol
            if isempty(tauopt)
                tauopt = map(T,M,VT,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
            end
            delta = tauopt-taubar;
            msqe(i) = sum(delta.^2,'all')/(length(T)*length(M));
            disp(['Mean square error = ' num2str(msqe(i))])
            if ~isempty(savePoints == i)
                recordedTaubars = cat(3,recordedTaubars,taubar);
            end
            if msqe(i) > tol
                [taubar,tauopt,stepScalar(i)] = feval(steppingFun,delta,taubar,msqe(i),T,M,VT,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,varargin{:});
                disp(['Iteration no. ' num2str(i) ' of main loop complete'])
            end
        else
            msqe = msqe(1:i-1);
            stepScalar = stepScalar(1:i-1);
            tau = taubar;
            disp('Finished to tolerence!')
            return
        end
    end
    tau = taubar;
    disp('Max iterations reached')
catch er
    tau = taubar;
    msqe = msqe(1:i-1);
    stepScalar = stepScalar(1:i-1);
    disp(['There was an error on iteration ' num2str(i) '!'])
    fprintf(2,'The message was:\n%s \n',er.message);
end