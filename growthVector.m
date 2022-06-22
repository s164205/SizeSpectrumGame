function gma = growthVector(N,M,rho,gma0,tau,p,offset)

%p = row vector (odd amount of entries) of the kernel around offset that
%can be eaten
%offset = the amount of mass partitions smaller a fish has to be (center of
%kernel) to be eaten
%rho is the relavtive change of mass at absorbtion, i.e. g = rho*mu

K = length(p);
k = K/2+.5;
NN = [zeros(1,offset+k-1),N];
MM = [zeros(1,offset+k-1),M];
tautau = [ones(1,offset+k-1)*tau(1),tau];
gma = zeros(1,length(N));
for i = 1:length(N)-1
    gma(i) = rho*sum(NN(i:i+K-1).*MM(i:i+K-1).*tautau(i:i+K-1).*p) + gma0(i);
end