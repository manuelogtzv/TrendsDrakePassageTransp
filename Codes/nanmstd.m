function [m, s, n] = nanmstd(x, varargin);
% function [m, s, n] = nanmstd(x, dim, squeeze);
% Mean, Standard Deviation, and Number of good points
% in an array with NaN-filled gaps.  By default calculations are
% done column-wise or along the first non-singleton dimension, as with mean, std. 
%
%    dim: optional argument; dimension along which to calculate
%
%    squeeze: optional string argument; if 'squeeze', then
%           the squeeze operator is applied to the outputs.
%
% YLF 2012/06, based on mstdgap by EF
% but can calculate along dimension N even with Mx1 input
% also uses unbiased stdev estimator (normalized by N-1)
% keeps things as arrays rather than working only with matrices so may be slower
% YLF 2012/09 - added if statement lines 24-28 to return scalar

xdim = size(x);

squeeze_flag = 0;
if nargin==1; 
   %dimension not specified, average along 1st non-singleton dimension
   %dim = find(xdim>1); dim = dim(1); %commented out 09-2012
   dim = find(xdim>1);
   if length(dim)==0
      dim = 1;
   end
   dim = dim(1);

elseif nargin>1
   for ii = 1:length(varargin)
      arg = varargin{ii};
      if strcmp(arg, 'squeeze'),
         squeeze_flag = 1;
      elseif arg > 0
	 dim = arg;
      else
         error('invalid optional argument')
      end
   end
end

if length(xdim)<dim
   %nothing to average over along this dimension
   m = x;
   s = zeros(xdim);
   n = s+1;

else

   xmsk = isnan(abs(x));
   n = sum(~xmsk, dim);
   nmsk = (n==0);
   ns = n-1;
   nsmsk = (ns<=0);
   badval = NaN; if ~isreal(x); badval = NaN+i*NaN; end

   %mean
   x(xmsk) = 0;
   m = sum(x, dim)./n;
   m(nmsk) = NaN;
   if ~isreal(x)
      m(nmsk) = m(nmsk) + i*NaN;
   end
      
   %stdev
   mrdim = ones(1, length(xdim)); mrdim(dim) = xdim(dim);
   xp = x - repmat(m, mrdim);
   xp(xmsk) = 0;
   s = sqrt(sum(xp.*conj(xp), dim)./ns);
   s(nsmsk) = NaN;

end

%optional squeeze
if squeeze_flag == 1,
   m = squeeze(m);
   s = squeeze(s);
   n = squeeze(n);
end
