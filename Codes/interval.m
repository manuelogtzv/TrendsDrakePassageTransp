% I = interval(A, MINSTEP) returns the groups of running indices between
% minimun step separation MINSTEP
% for instance, if A is a list of dates and MINSTEP is 5, then I is a cell
% arrayof strings containing the indices of each series of days separated
% by less than 5 days.
% I =  interval(A, MINSTEP,'n') returns a cell array of the indices of each 
% continuous series. 

function I = interval(A, minstep, varargin)

op = 's';
if size(varargin)>0
	op = varargin{1};
end

B = find(diff(A)>minstep);
I =  [repmat('[',length(B)+1,1),num2str([1;B(:)+1]),repmat(':',length(B)+1,1),num2str([B(:);length(A)]),...
     repmat(']',length(B)+1,1)];
I = cellstr(I);

if op == 'n'
	A  = cell(1,length(I));
	for n = 1:length(I);
		A{n} = eval(I{n});
	end
	I = A;
end

