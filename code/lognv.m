function vector=lognv(n,min,max,mean,sd)
mu = log((mean^2)/sqrt(sd+mean^2));
sigma = sqrt(log(sd/(mean^2)+1));
dist = makedist('LogNormal','mu',mu,'sigma',sigma);
trdist=truncate(dist,min,max);
vector=random(trdist,[n,1]);
end


