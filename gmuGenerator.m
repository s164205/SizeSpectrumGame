function G = gmuGenerator(g,mu)
%g,mu are col. vectors of resp. growth- and deathrates

n = length(g);
G = spdiags([zeros(1,n); -mu-g; [0,g(1:end-1)]]',-1:1,n,n);

