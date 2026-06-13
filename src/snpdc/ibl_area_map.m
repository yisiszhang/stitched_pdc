function areaMap = ibl_area_map()
%IBL_AREA_MAP Allen CCF remapping used in the IBL spontaneous analysis.

pairs = {
    'DG-sg','DG'; 'DG-mo','DG'; 'DG-po','DG'; 'DG','DG'
    'CA1','CA1'; 'CA2','CA2'; 'CA3','CA3'; 'SUB','SUB'; 'ProS','ProS'; 'POST','POST'; 'PAR','PAR'; 'PRE','PRE'; 'HATA','HATA'; 'FC','FC'
    'VISa6a','VISa'; 'VISa6b','VISa'; 'VISa5','VISa'; 'VISa4','VISa'; 'VISa2/3','VISa'; 'VISa1','VISa'
    'VISp6a','VISp'; 'VISp6b','VISp'; 'VISp5','VISp'; 'VISp4','VISp'; 'VISp2/3','VISp'; 'VISp1','VISp'
    'VISam1','VISam'; 'VISam6a','VISam'; 'VISam6b','VISam'; 'VISam5','VISam'; 'VISam4','VISam'; 'VISam2/3','VISam'
    'VISpm1','VISpm'; 'VISpm6a','VISpm'; 'VISpm6b','VISpm'; 'VISpm5','VISpm'; 'VISpm4','VISpm'; 'VISpm2/3','VISpm'
    'VISl6a','VISl'; 'VISl6b','VISl'; 'VISl5','VISl'; 'VISl2/3','VISl'
    'VISal1','VISal'; 'VISal6a','VISal'; 'VISal2/3','VISal'; 'VISli5','VISli'; 'VISpl6a','VISpl'
    'VISpor5','VISpor'; 'VISpor6a','VISpor'; 'VISrl4','VISrl'; 'VISrl5','VISrl'; 'VISrl6b','VISrl'
    'RSPv1','RSPv'; 'RSPv2/3','RSPv'; 'RSPv5','RSPv'; 'RSPv6a','RSPv'; 'RSPv6b','RSPv'
    'RSPd1','RSPd'; 'RSPd6a','RSPd'; 'RSPd6b','RSPd'; 'RSPd5','RSPd'; 'RSPd2/3','RSPd'
    'RSPagl1','RSPagl'; 'RSPagl6a','RSPagl'; 'RSPagl6b','RSPagl'; 'RSPagl5','RSPagl'; 'RSPagl2/3','RSPagl'
    'SCig','SCm'; 'SCiw','SCm'; 'SCdg','SCm'; 'SCdw','SCm'; 'SCop','SCs'; 'SCsg','SCs'; 'SCzo','SCs'
    'VPM','VPM'; 'VPMpc','VPMpc'; 'VPL','VPL'; 'VPLpc','VPLpc'; 'LP','LP'; 'PO','PO'; 'POL','POL'
    'LGd-co','LGd'; 'LGd-ip','LGd'; 'LGd-sh','LGd'; 'LGv','LGv'
    'MGv','MG'; 'MGd','MG'; 'MGm','MG'; 'MD','MD'; 'VM','VM'; 'VL','VL'; 'CM','CM'; 'CL','CL'; 'CUN','CUN'
    'LD','LD'; 'LH','LH'; 'PF','PF'; 'PVT','PVT'; 'PCN','PCN'; 'SMT','SMT'; 'TH','TH'; 'RT','RT'
    'PoT','PoT'; 'SGN','SGN'; 'VAL','VAL'; 'AMd','AMd'; 'AMv','AMv'; 'AV','AV'; 'AD','AD'; 'IAD','IAD'; 'PIL','PIL'; 'PP','PP'; 'MH','MH'; 'IntG','IntG'; 'RH','RH'; 'PR','PR'; 'TRN','TRN'
    'ZI','ZI'; 'Eth','Eth'; 'SPFp','SPFp'
    'MRN','MRN'; 'APN','APN'; 'PAG','PAG'; 'RN','RN'; 'SNr','SNr'; 'SNc','SNc'; 'MB','MB'; 'NOT','NOT'
    'PPT','PPT'; 'PPN','PPN'; 'NPC','NPC'; 'RPF','RPF'; 'OP','OP'; 'CLI','CLI'
    'NB','NB'; 'DR','DR'; 'RR','RR'; 'INC','INC'; 'ND','ND'; 'DT','DT'; 'PBG','PBG'; 'LDT','LDT'; 'SOCl','SOCl'; 'VTN','VTN'; 'EW','EW'; 'SCO','SCO'; 'Su3','Su3'
    'MOp6a','MOp'; 'MOp6b','MOp'; 'MOp5','MOp'; 'MOp2/3','MOp'; 'MOp1','MOp'
    'MOs6a','MOs'; 'MOs6b','MOs'; 'MOs5','MOs'; 'MOs2/3','MOs'; 'MOs1','MOs'
    'SSp-bfd6a','SSp'; 'SSp-bfd6b','SSp'; 'SSp-bfd5','SSp'; 'SSp-bfd4','SSp'; 'SSp-bfd2/3','SSp'
    'SSp-tr6a','SSp'; 'SSp-tr6b','SSp'; 'SSp-tr5','SSp'; 'SSp-tr4','SSp'; 'SSp-tr2/3','SSp'
    'SSp-m6a','SSp'; 'SSp-m6b','SSp'; 'SSp-m5','SSp'; 'SSp-m4','SSp'
    'SSp-ul6a','SSp'; 'SSp-ul5','SSp'; 'SSp-ul4','SSp'; 'SSp-ul2/3','SSp'; 'SSp-ll6a','SSp'; 'SSp-ll6b','SSp'; 'SSp-ll5','SSp'; 'SSp-ll4','SSp'; 'SSp-n6a','SSp'; 'SSp-n5','SSp'; 'SSp-un6a','SSp'
    'SSs6a','SSs'; 'SSs6b','SSs'; 'SSs5','SSs'; 'SSs4','SSs'; 'SSs2/3','SSs'
    'ACAd6a','ACA'; 'ACAd6b','ACA'; 'ACAd5','ACA'; 'ACAd2/3','ACA'; 'ACAv6a','ACA'; 'ACAv6b','ACA'; 'ACAv5','ACA'; 'ACAv2/3','ACA'
    'PL6a','PL'; 'PL6b','PL'; 'PL5','PL'; 'PL2/3','PL'
    'ILA6a','ILA'; 'ILA5','ILA'; 'ILA2/3','ILA'; 'ILA1','ILA'
    'ORBm6a','ORB'; 'ORBm5','ORB'; 'ORBl5','ORB'; 'ORBl2/3','ORB'; 'ORBl1','ORB'; 'ORBvl6a','ORB'; 'ORBvl5','ORB'; 'ORBvl2/3','ORB'; 'FRP6a','FRP'; 'FRP5','FRP'
    'CP','CP'; 'ACB','ACB'; 'STR','STR'; 'GPe','GPe'; 'GPi','GPi'; 'SI','SI'; 'LSc','LSc'; 'PAL','PAL'
    'AUDd5','AUD'; 'AUDd2/3','AUD'; 'AUDpo2/3','AUD'; 'AUDp4','AUDp'; 'AUDp5','AUDp'; 'AUDp6a','AUDp'; 'AUDpo4','AUDpo'; 'AUDpo6a','AUDpo'; 'AUDpo6b','AUDpo'; 'AUDv5','AUDv'; 'AUDv6a','AUDv'
    'IRN','IRN'; 'GRN','GRN'; 'MV','MV'; 'NTS','NTS'; 'NLL','NLL'
    'ICe','IC'; 'PB','PB'; 'PSV','PSV'; 'MY','MY'; 'P','P'
    'CB','CB'; 'IP','IP'; 'CENT2','CENT2'; 'CENT3','CENT3'; 'ANcr1','ANcr1'; 'ANcr2','ANcr2'; 'NOD','NOD'; 'FL','FL'; 'COPY','COPY'; 'VeCB','VeCB'; 'SIM','SIM'; 'PRM','PRM'; 'PFL','PFL'
    'fiber tracts','fiber tracts'; 'fp','fp'; 'dhc','dhc'; 'cpd','cpd'; 'ccs','ccs'; 'or','or'; 'alv','alv'; 'ec','ec'; 'ml','ml'; 'int','int'; 'arb','arb'
    'mcp','mcp'; 'dscp','dscp'; 'sctv','sctv'; 'scwm','scwm'; 'cing','cing'; 'fa','fa'; 'ff','ff'; 'FF','FF'; 'icp','icp'; 'cst','cst'
    'ccb','ccb'; 'V','V'; 'll','ll'; 'em','em'; 'scp','scp'; 'nst','nst'
    'BLAa','BLA'; 'BLAp','BLA'; 'BMAp','BMA'; 'BMAa','BMA'; 'CEAc','CEA'; 'CEAm','CEA'; 'CEAl','CEA'
    'LA','LA'; 'AAA','AAA'; 'COApm','COA'; 'COApl','COA'; 'COAa','COA'; 'MEA','MEA'; 'PAA','OLF'; 'PA','PA'
    'AON','AON'; 'OLF','OLF'; 'PIR','PIR'; 'TR','TR'; 'OT','OT'; 'IA','IA'; 'TTd','OLF'; 'DP','OLF'
    'AId6a','AI'; 'AId6b','AI'; 'AId5','AI'; 'AIv6a','AI'; 'AIv5','AI'; 'AIv2/3','AI'; 'AIp6a','AI'; 'AIp6b','AI'; 'GU6a','GU'; 'ECT5','ECT'; 'ECT2/3','ECT'; 'TEa6a','TEa'; 'TEa5','TEa'; 'PRC','PRC'
    'HY','HY'; 'LHA','LHA'; 'SF','SF'; 'BST','BST'; 'MPO','MPO'; 'MS','MS'; 'NDB','NDB'; 'SAG','SAG'; 'PSTN','PSTN'; 'STN','STN'; 'PH','PH'; 'TU','TU'; 'PeF','PeF'; 'AVP','AVP'; 'DMH','DMH'; 'ADP','ADP'
    'CTXsp','CTXsp'
    'ENTl6a','ENTl'; 'ENTl5','ENTl'; 'ENTl3','ENTl'; 'ENTl2','ENTl'
    'ENTm3','ENTm'; 'ENTm2','ENTm'; 'ENTm6','ENTm'; 'ENTm1','ENTm'; 'ENTm5','ENTm'
    'SPVI','SPVI'; 'TRS','TRS'; 'SUT','SUT'; 'MARN','MARN'; 'PCG','PCG'; 'III','III'; 'ACVII','ACVII'; 'VII','VII'; 'LAV','LAV'; 'P5','P5'; 'SPVC','SPVC'; 'ICB','ICB'; 'ICc','ICc'; 'IMD','IMD'; 'SPIV','SPIV'; 'SUM','SUM'; 'y','y'; 'APr','APr'
    'CLA','CLA'; 'EPd','EP'; 'EPv','EP'; 'LSr','LS'; 'LSv','LS'; 'HPF','HPF'; 'DCO','DCO'; 'DN','DN'; 'FS','FS'; 'CS','CS'
    'VISpm','VISpm'; 'VISa','VISa'; 'VISp','VISp'; 'root','root'
    };

excludeAreas = {'fiber tracts', 'fp', 'dhc', 'cpd', 'ccs', 'or', 'alv', 'ec', 'ml', 'int', 'arb', ...
    'mcp', 'dscp', 'sctv', 'scwm', 'cing', 'fa', 'ff', 'FF', 'icp', 'cst', 'ccb', 'V', 'll', 'em', 'scp', 'nst', ...
    'ccg', 'bsc', 'fi', 'tspc', 'mlf', 'rust', 'bic', 'ee', 'mtg', 'pc', 'sm', 'sptV', 'aco', 'ar', 'opt', 'uf', 'fr', 'csc', 'tb', 'act', 'py', 'st', 'IVn', ...
    'void', 'chpl', 'V3', 'V4', 'HPF', 'root'};

areaMap = struct();
areaMap.map = containers.Map(pairs(:,1), pairs(:,2));
areaMap.exclude = excludeAreas;
areaMap.allowed_outputs = unique(string(pairs(:,2)));
end
