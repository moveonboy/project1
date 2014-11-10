function varargout= readCellNGlycanFromExcel(i,nglycanname)
%
%
%
%
%
%
%
%
%
%
%
%
%
% Author: Yusen Zhou & Gang Liu 
% Date Lastly Updated: 10/27/14
nargoutchk(1,5)
[peak,glycancomp,rawdata]   = xlsread(nglycanname);
nglycancomplist                       = glycancomp(2:end,:);
nglycanmasslist                        = peak(:,1);
chowildtype                              = peak(:,i);
index                                          = isfinite(chowildtype);
chowildtypecomps                  = nglycancomplist(index);
chowildtypepeaks                    = nglycanmasslist(index);
% glycanstring = gly3to1letter();

glycanstringarray     = cellfun(@gly1charformat,chowildtypecomps,...
    'UniformOutput', false);
glycanformulaarray = cellfun(@glycanFormula,glycanstringarray);
glycanmwarray        = arrayfun(@(x)isotopicdist(x,'SHOWPLOT',false),...
    glycanformulaarray,'UniformOutput', false);

if(nargout==1)
    varargout{1} = chowildtypecomps;
elseif(nargout==2)
    varargout{1} = chowildtypecomps;
    varargout{2} = glycanmwarray;
elseif(nargout==3)
    varargout{1} = chowildtypecomps;
    varargout{2} = glycanmwarray;
    varargout{3} = chowildtypepeaks;
elseif(nargout==4)
    varargout{1} = chowildtypecomps;
    varargout{2} = glycanmwarray;
    varargout{3} = chowildtypepeaks;
    varargout{4} = glycanstringarray;
elseif(nargout==5)
    varargout{1} = chowildtypecomps;
    varargout{2} = glycanmwarray;
    varargout{3} = chowildtypepeaks;
    varargout{4} = glycanstringarray;
    varargout{5} = glycanformulaarray;
end
% disp('end');
end


