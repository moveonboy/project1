function blockdatabase = createblocks()
%createblocks create the building blocks database by using the
% containers.Map
%
% blockdatabase = createblocks() using the containers.Map to build the database, 
% the key is the name of block and the mapObj is the structure of the block.
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 9/23/14

blockdatabase = containers.Map;
blockdatabase('LacNAcI')       = glycanMLread('LacNAcI.glycoct_xml');
blockdatabase('SLacNAcI')      = glycanMLread('SLacNAcI.glycoct_xml');
blockdatabase('a2_6SLacNAcI')  = glycanMLread('a2_6SLacNAcI.glycoct_xml');
blockdatabase('LacNAcII')      = glycanMLread('LacNAcII.glycoct_xml');
blockdatabase('SLacNAcII')     = glycanMLread('SLacNAcII.glycoct_xml');
blockdatabase('a2_6SLacNAcII') = glycanMLread('a2_6SLacNAcII.glycoct_xml');
blockdatabase('LacNAcIII')     = glycanMLread('LacNAcIII.glycoct_xml');
blockdatabase('LacNAcIV')      = glycanMLread('LacNAcIV.glycoct_xml');
blockdatabase('SLeX')          = glycanMLread('SLeX.glycoct_xml');
blockdatabase('LeX')           = glycanMLread('LeX.glycoct_xml');
blockdatabase('a2_6SSLeX')     = glycanMLread('a2_6SSLeX.glycoct_xml'); %a2/6 similar sialyl lewis X
blockdatabase('coreFuc')       = glycanMLread('coreFuc.glycoct_xml');
blockdatabase('bisecting')     = glycanMLread('bisecting.glycoct_xml');
blockdatabase('core1')         = glycanMLread('core1.glycoct_xml');
blockdatabase('core2')         = glycanMLread('core2.glycoct_xml');
blockdatabase('core3')         = glycanMLread('core3.glycoct_xml');
blockdatabase('core4')         = glycanMLread('core4.glycoct_xml');
blockdatabase('core5')         = glycanMLread('core5.glycoct_xml');
blockdatabase('core6')         = glycanMLread('core6.glycoct_xml');
blockdatabase('core7')         = glycanMLread('core7.glycoct_xml');
blockdatabase('core8')         = glycanMLread('core8.glycoct_xml');
end