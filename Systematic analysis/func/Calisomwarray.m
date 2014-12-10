function [chowildtypecomps,glycanmwarray] = Calisomwarray(glycancomp)
% Calglycanmwarray calculate all the isotopic mass based on the composition
%     string.
%
% Syntax:
%       [chowildtypecomps,glycanmwarray] = Calglycanmwarray(glycancomp)
%
% Input:
%
%
% Output:
%
%
% Example:
%
%
%
%

% Author: Yusen Zhou
%
if(~iscellstr(glycancomp))
    error('WRONG INPUT TYPE')
end
chowildtypecomps = glycancomp;
glycanstringarray     = cellfun(@gly1charformat,chowildtypecomps,...
    'UniformOutput', false);
glycanformulaarray    = cellfun(@glycanFormula,glycanstringarray);
glycanmwarray         = arrayfun(@(x)isotopicdist(x,'SHOWPLOT',false),...
    glycanformulaarray,'UniformOutput', false);
end