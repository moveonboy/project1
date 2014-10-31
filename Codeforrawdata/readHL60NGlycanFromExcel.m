function [chowildtypecomps,chowildtypepeaks,glycanstringarray,...
    glycanformulaarray,glycanmwarray] = readHL60NGlycanFromExcel(i)
nglycanname                            = 'HL60WT_glycanspeceis.xls';
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

% disp('end');
end


