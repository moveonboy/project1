function varargout = msrelativeabundance(MSdata,nglycanname,abundanceExcelname)
%msrelativeabundance: Calculate the relative abundance for each potential
% glycans based on the MSdata and potential glycans, and write the result
% into a Excel file.
%
%glycanlist = msrelativeabundance(MSdata,nglycanname,abundanceExcelname)
% Output the relative abundance Excel file and glycanlist of MSdata
%
%[glycanlist abundancelibrary] = msrelativeabundance(MSdata,nglycanname,abundanceExcelname)
% Output the relative abundance Excel file, potential glycanlist of MSdata
% and containers.Map for glycan relative abundance.
%
% For example:
%         [glycanlist isopeakareagroup,abundancelibrary] = msrelativeabundance('MSdata.mat',..
%                                   'HL60WT_glycanspeceis.xls','HL60WTN_relativeabundance.xls');
%
%
%
% Author: Yusen Zhou 
% Date Lastly Updated: 10/27/14
nargoutchk(1,2);


load(MSdata);
peaklist = MSdata.peaklist;   
pfwhh    = MSdata.pfwhh;
[chowildtypecomps,glycanmwarray] = readCellNGlycanFromExcel(1,nglycanname);
[glycanlist,abundancelibrary]= msfraction(peaklist,pfwhh,chowildtypecomps,glycanmwarray,abundanceExcelname);

if(nargout==1)
    varargout{1} = glycanlist;
elseif(nargout==2)
    varargout{1} = glycanlist;
    varargout{2} = abundancelibrary;
    save('HL60WTlibrary.mat','HL60library');
end

end