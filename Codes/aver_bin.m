function [xb, yb] = aver_bin(x, y, gsize);
%function [xb, yb] = line_bin(x, y, gsize);
%
%averages y in gsize-length bins along a line
%x is 1xN, y is MxN, gsize is 1x1
%
%will lose from the end if gsize doesn't fit evenly into range(x)
% se changed code to match YLF's, 06-05-12, to use YLF's nanmstd.m for statistics


xb = [min(x)+gsize/2:gsize:max(x)+gsize/2];
yb = NaN*ones(length(xb), 1);
for no = 1:length(xb)
   ii = find(x>=xb(no)-gsize/2 & x<xb(no)+gsize/2);
   %yb(:, no) = nanmean(y(:, ii), 2); % ylf
   %yb(:,no) =mstdgap(y(:, ii), 2);    % codas
   yb(no, 1) = nanmstd(y(ii));
end

