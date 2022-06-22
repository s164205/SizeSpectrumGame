function nu = deathVector(N,rho,nu0,tau,p,offset)

%p = row vector (odd amount of entries) of the kernel around offset that
%can eat "me"
%offset = the amount of mass partitions bigger a fish has to be (center of
%kernel) to eat "me"
%nu0 = vector with length(N) entries

K = length(p);
k = K/2+.5;
NN = [N,zeros(1,offset+k-1)];
tautau = [tau,zeros(1,offset+k-1)];
%mu1 = 0; Back from when N(1) had special meaning
% for i = 2:offset+k %death of mu(1), that takes all indices < 1 "on its shoulders"
%      mu1 = mu1 + sum(NN(i:i+K-1).*tautau(i:i+K-1).*p)*tau(1);
% end
nu = zeros(1,length(N));
%nu(1) = mu1;
for i = 1:length(N)
    nu(i) = rho*sum(NN(offset-k+i+1:offset-k+i+K).*tautau(offset-k+i+1:offset-k+i+K).*p) + nu0(i);
end
    