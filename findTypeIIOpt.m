function tauopt = findTypeIIOpt(V,gma,nu,M,Cmax,eps,K)

lambda = Cmax.*K./(gma.*(eps*Cmax-K)); %taumin s.t. eps*C = K
lambda(end) = 0;
lambda(isnan(lambda)) = 0;

dVdm = [0, diff(V)]./[1,  diff(M)];
dVdm(1) = dVdm(2);
tauopt = zeros(1,length(V));
for i = 1:length(V)-1
    if dVdm(i) > 0 
        tauopt(i) = min(1,max(lambda(i),((sqrt(dVdm(i)*V(i)*nu(i)*gma(i)*eps) - nu(i)*V(i))*Cmax(i))/(gma(i)*nu(i)*V(i))));
    else
        tauopt(i) = lambda(i); %If dVdm negative, then minimise
        tauopt(V <= 0) = 1; %If V = 0 (dVdm is then probably = 0), then maximise 
    end
end