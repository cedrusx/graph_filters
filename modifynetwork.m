function Am = modifynetwork(A,p,METHOD,Loc,radius)

N = size(A,1);
if METHOD == 1 % radius increase by p
    Am = (loc2dis(Loc) < radius * (1+p)) - eye(N);
elseif METHOD == 2 % p link failures
    Am = A .* (random('unif',0,1,N,N) > p);
    Am = triu(Am) + triu(Am)' - diag(diag(Am));
elseif METHOD == 3 % for each sensor, real_radius = radius * unif(1-p,1+p) 
    dis = loc2dis(Loc,ones(N,N));
    Am = dis <= radius * random('unif',1-p,1+p,N,N);
    Am = triu(Am) + triu(Am)' - 2 * diag(diag(Am));
end