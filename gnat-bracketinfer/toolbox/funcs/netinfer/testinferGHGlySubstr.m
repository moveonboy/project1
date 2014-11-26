             mani  = GHEnz.loadmat('mani.mat');
             m5species  = glycanMLread('M5.glycoct_xml');
             m6species = inferGHGlySubstr(m5species,mani,0);
             options  = displayset('showmass',true,'showLinkage',true,...
                              'showRedEnd',true);
             for i = 1: m6species.length
                 glycanViewer(m6species.get(i).glycanStruct,options);
             end
 
             mani       = GHEnz.loadmat('mani.mat');
             m5species  = glycanMLread('M5.glycoct_xml');
             m6species  = inferGHGlySubstr(m5species,mani,1);
             options    = displayset('showmass',true,'showLinkage',true,...
                              'showRedEnd',true);
             for i = 1: m6species.length
                 glycanViewer(m6species.get(i).glycanStruct,options);
             end
 
             mani       = GHEnz.loadmat('mani.mat');
             m5species  = glycanMLread('M5.glycoct_xml');
             [nsubstr,m6species,m5rxns] = inferGHGlySubstr(m5species,mani,1);
             options  = displayset('showmass',true,'showLinkage',true,...
                              'showRedEnd',true);
             for i = 1: nsubstr
                glycanRxnViewer(m5rxns.get(i));
             end
 
             mani       = GHEnz.loadmat('mani.mat');
             m5species  = glycanMLread('M5.glycoct_xml');
             [nsubstr, m7species,m6rxns,m6pathway] = inferGHGlySubstr(m5species,mani,1);
             glycanPathViewer(m6pathway);
 
             mani  = GHEnz.loadmat('mani.mat');
             m6species  = glycanMLread('m6bracket.glycoct_xml');
             [nsubstr, m7species,m7rxns,m7pathway] = inferGHGlySubstr(m6species,mani,1);
             glycanPathViewer(m6pathway);
             
             mani  = GHEnz.loadmat('mani.mat');
             m6species  = glycanMLread('m7bracket.glycoct_xml');
             [nsubstr, m8species,m8rxns,m8pathway] = inferGHGlySubstr(m7species,mani,1);
             glycanPathViewer(m6pathway);
             
             mani  = GHEnz.loadmat('mani.mat');
             m6species  = glycanMLread('m8bracket.glycoct_xml');
             [nsubstr, m9species,m9rxns,m9pathway] = inferGHGlySubstr(m8species,mani,1);
             glycanPathViewer(m6pathway);
             