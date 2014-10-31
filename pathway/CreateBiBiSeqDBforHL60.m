function bibiSeqKineticsDB = CreateBiBiSeqDBforHL60()


%% load enzdatabase
enzdbmatfilename   = 'HL60enzDB.mat';
enzdb              = enzdbmatLoad(enzdbmatfilename);

%% load enzyme
mgat1     = enzdb('mgat1');
mgat2     = enzdb('mgat2');
mgat3     = enzdb('mgat3');
mgat4     = enzdb('mgat4');
mgat5     = enzdb('mgat5');
mani      = enzdb('mani');
b4GalI    = enzdb('b4GalI');
manii     = enzdb('manii');
Fut8      = enzdb('Fut8');
ignt      = enzdb('ignt');
Fut7      = enzdb('Fut7');
ST3GalIII = enzdb('ST3GalIII');
ST3GalIV  = enzdb('ST3GalIV');

%% setup EnzKinetics 
rt1                       =  struct('Vm',450,'Kmd_G1E',0.65,'Kmd_SNE',0,'Keq',1);
enzkineticsobj            =  BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt1);
mgat1.enzkinetics         = enzkineticsobj;
mgat1.enzkinetics.enz     = mgat1;

rt2                       = struct('Vm',140,'Kmd_G1E',0.475,'Kmd_SNE',2.4,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt2);
mgat2.enzkinetics         = enzkineticsobj;
mgat2.enzkinetics.enz     = mgat2;

rt3                       = struct('Vm',4000,'Kmd_G1E',0.475,'Kmd_SNE',7.75,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt3);
mgat3.enzkinetics         = enzkineticsobj;
mgat3.enzkinetics.enz     = mgat3;

rt4                       = struct('Vm',10,'Kmd_G1E',8.5,'Kmd_SNE',20.75,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt4);
mgat4.enzkinetics         = enzkineticsobj;
mgat4.enzkinetics.enz     = mgat4;

rt5                       = struct('Vm',10,'Kmd_G1E',0.325,'Kmd_SNE',8.75,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt5);
mgat5.enzkinetics         = enzkineticsobj;
mgat5.enzkinetics.enz     = mgat5;

rt6                       = struct('Vm',300,'Km',0.25);
enzkineticsobj6           = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt6);
mani.enzkinetics          = enzkineticsobj6;
mani.enzkinetics.enz      = mani;

rt7                       = struct('Vm',580,'Kmd_G1E',0.325,'Kmd_SNE',0,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt7);
b4GalI.enzkinetics        = enzkineticsobj;
b4GalI.enzkinetics.enz    = b4GalI;

rt8                       = struct('Vm',300,'Km',0.25);
enzkineticsobj            = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt8);
manii.enzkinetics         = enzkineticsobj;
manii.enzkinetics.enz     = manii;

rt9                       = struct('Vm',38,'Kmd_G1E',0.0625,'Kmd_SNE',0.115,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt9);
Fut8.enzkinetics          = enzkineticsobj;
Fut8.enzkinetics.enz      = Fut8;

rt10                      = struct('Vm',38,'Kmd_G1E',0.0625,'Kmd_SNE',0.115,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt10);
ignt.enzkinetics          = enzkineticsobj;
ignt.enzkinetics.enz      = ignt;

rt11                      = struct('Vm',38,'Kmd_G1E',7.7,'Kmd_SNE',0.041,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt11);
Fut7.enzkinetics          = enzkineticsobj;
Fut7.enzkinetics.enz      = Fut7;

rt12                      = struct('Vm',70,'Kmd_G1E',7.5,'Kmd_SNE',0.075,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt12);
ST3GalIII.enzkinetics     = enzkineticsobj;
ST3GalIII.enzkinetics.enz = ST3GalIII;

rt13                      = struct('Vm',70,'Kmd_G1E',0.55,'Kmd_SNE',0.075,'Keq',1);
enzkineticsobj            = BiBiKinetics(TypeStringConst.bibiSeqOrderwSubstInhib,rt13);
ST3GalIV.enzkinetics      = enzkineticsobj;
ST3GalIV.enzkinetics.enz  = ST3GalIV;

bibiSeqKineticsDB = containers.Map;
bibiSeqKineticsDB('mgat1')     = mgat1.enzkinetics;
bibiSeqKineticsDB('mgat2')     = mgat2.enzkinetics;
bibiSeqKineticsDB('mgat3')     = mgat3.enzkinetics;
bibiSeqKineticsDB('mgat4')     = mgat4.enzkinetics;
bibiSeqKineticsDB('mgat5')     = mgat5.enzkinetics;
bibiSeqKineticsDB('mani')      = mani.enzkinetics;
bibiSeqKineticsDB('manii')     = manii.enzkinetics;
bibiSeqKineticsDB('b4GalI')    = b4GalI.enzkinetics;
bibiSeqKineticsDB('Fut8')      = Fut8.enzkinetics;
bibiSeqKineticsDB('ignt')      = ignt.enzkinetics;
bibiSeqKineticsDB('Fut7')      = Fut7.enzkinetics;
bibiSeqKineticsDB('ST3GalIII') = ST3GalIII.enzkinetics;
bibiSeqKineticsDB('ST3GalIV')  = ST3GalIV.enzkinetics;
