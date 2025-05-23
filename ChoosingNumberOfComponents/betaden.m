function y = betaden(x,s,m)

if nargin < 3
   error('stats:betapdf:TooFewInputs','Requires three input arguments.');
end

% Return NaN for out of range parameters.
s(s<=0) = NaN;
m(m<=0) = NaN;

% Out of range x could create a spurious NaN*i part to y, prevent that.
% These entries will get set to zero later.
xOutOfRange = (x<0) | (x>1);
x(xOutOfRange) = .5;

logkerna = ((s.*m)-1).*log(x);   
logkernb = ((s-s.*m)-1).*log(1-x); 
    y = exp(logkerna+logkernb - betaln(s.*m,s-s.*m));