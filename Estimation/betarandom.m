%generate random numbers from beta distribution

function rnd = betarandom(s,m)
%s = a+b,m=a/s

if (nargin ~= 2)
    error('betarnd: you must give two arguments');
	end
	
	if (~isscalar(s) || ~isscalar(m))
	  error('betarnd: A and B should be scalar parameters');
	end
	
	if (s <= 0 || s == Inf || m <= 0 || m == Inf)
	  rnd = NaN;
	else
	  x = gamrnd(s*m, 1);
	  rnd = x/(x+gamrnd(s-s*m, 1));
	end


