

%LOC='/usr2/home/matlab/OPA-DA/V0-static-data/';
%LOC='/usr2/home/matlab/OPA-DA/V0_4.9-static-data/';
%LOC='/usr2/home/matlab/OPA-DA/V1_4.9-static-data/';
%LOC='/home/matlab/OPA-DA/V4_4.9-static-data/CCI_CLIM_CMEMS/';
LOC='/home/matlab/OPA-DA/V5_4.9-static-data/CCI1km_CMEMScoast/';

% IN LOC VANNO I RISULTATI (3D_VAR e MISFIT)

% bisogna avere le due griglie, fornite da fortran90
% per SUB e SUB_ALL

%% varSat
%verificare che non sia già esistente
%dirdata='/storage/matlab/Archive_satellite/CCI_CNR/';
%dirdata='/storage/matlab/Archive_satellite/CCI_CNR_CMEMS/CHECKED/';
%dirdatasat='/storage/matlab/Archive_satellite/NEW_CCI_10days';
dirdatasat='/storage/matlab/Archive_satellite/CCI1km_Interp24_10days';
limstd=2;  % Nella funzione usati solo
dayF=10;   % per nomi dei file
%W=SWIFS_mask16;
%W=V4_mask;
W=V5_mask;
%climfile='/storage/matlab/MOVED/Archive_Climatologies/Satellite/CCI02/CLIMATOLOGY_15.mat';
%varsat_on_origmesh='/home/matlab/OPA-DA/VAR_SAT/CCI_CMEMS/';
%dirvarsat='/home/matlab/OPA-DA/VAR_SAT/CCI_CMEMS/';
dirvarsat='/home/matlab/OPA-DA/VAR_SAT/CCI1km_Interp24/';
%var_sat_CCI_10gg(limstd,dayF,W,dirdata,climfile,varsat_on_origmesh)
var_sat_CCI_10gg(limstd,dayF,W,dirdatasat,dirvarsat)
%questa funzione al momento è qui
%/home/matlab/OPA-DA/VAR_SAT

%% Increase VarSat in summer
% the results should NOT be used for evaluation of varmod

VV=V5_mask;
%dirVS = '/home/matlab/OPA-DA/VAR_SAT/CCI_CMEMS';
dirVS = '/home/matlab/OPA-DA/VAR_SAT/CCI1km_Interp24';
%outDIR   = '/home/matlab/OPA-DA/VAR_SAT/CCI_CMEMS/SUMMER_INCREASED50_100/'; 
outDIR   = '/home/matlab/OPA-DA/VAR_SAT/CCI1km_Interp24/SUMMER_INCREASED50_100/'; 
% increasing factor hardcoded!!! <--

increaseSummerSatVar(dirVS,outDIR,VV) %4 mappe mensili
% la funzione è qui 
% /home/matlab/OPA-DA/VAR_SAT


%% varErr
VV=V5_mask;
dirVE='/home/matlab/OPA-DA/VAR_SAT_MOD/V5/VARERR_CCI1km_CMEMScoast'; %the output
dirMod='/myo1/CMEMScoast/WEEKLY/CHL_SUP24';
dirSat='/storage/matlab/Archive_satellite/CCI1km_Interp24_weekly'; %calcolata altrove
% è da ricalcolare sui checked/
%diff_CCINEW_RACMEMS(dirVE,dirMod,dirSat,VV)
diff_CCI1km_CMEMScoast(dirVE,dirMod,dirSat,VV)
%diff_CCI_RA3(dirOutput,dirMod,dirSat,VV)
% questa funzione al momento è qui
% /home/matlab/OPA-DA/VAR_SAT_MOD/V5

%% creo i var2D.12.nc per le eof SUB_ALL
limGib=-4.9;
VV=V5_mask;
minMod=0.5;
%minModC=0.5; %coast v1
minModC=3/4; %coast v2
dirVS = '/home/matlab/OPA-DA/VAR_SAT/CCI1km_Interp24';
outDIR   = [LOC 'VAR_MOD/CCI1km_CMEMScoast']; 

varMod_CMEMScoast_CCI1km(dirVE,dirVS,limGib,VV,minMod,minModC,outDIR) %12 mappe mensili
% questa funzione al momento è qui
% /home/matlab/OPA-DA/VAR_SAT_MOD/V5

%varMod_RA3_CCI(dirVE, DIR_INT,limGib,VV,MaskVersion,minErrMod,outDIR) %12 mappe mensili
% questa funzione al momento è qui
% /home/matlab/annared/TESTVARMOD

%% creo gli eof.9.

% con KB
gridFile = [LOC '3D_VAR/GRID/SUB/BFM_grid.nc']; 
outDIR   = [LOC '3D_VAR/EOF/' ];
readRUN_doEOFs_GRID(gridFile,outDIR); % scrive in SUB/eof.9.01.nc

%%
%con rea OPEC
LOCV1 = '/usr2/home/matlab/OPA-DA/V1_4.9-static-data/';
gridVar = [LOCV1 '3D_VAR/GRID/SUB/BFM_grid.nc']; 
outDir   = [LOC '3D_VAR/EOF/SUB/RA03'];
inDir    = '/myo1/OPEC/OPA_872_RA-R03/AVE';
readRUN_doEOFs(gridVar,inDir,outDir)

%%
% con rea CMEMScoast 
gridVar = [LOC '3D_VAR/GRID/SUB15/BFM_grid15.nc'];%griglia a 1/16 
outDir  = [LOC '3D_VAR/EOF'];               %15 sottobacini perché UN Adriatico
inDir   = '/myo1/CMEMScoast/WEEKLY/AVE';    %non ci sarebbero abbastanza punin adr1
var = 'P_l'                           %per calcolo EOF
readRUN_doEOFs(gridVar,inDir,outDir,var)

%%
% con rea CMEMS per float (EOF costanti per 15 subbasin)
gridVar = [LOC '3D_VAR/GRID/SUB15_16/BFM_grid.nc'];
outDir  = [LOC '3D_VAR/EOF/TEST_WEEKLY/'];
inDir   = '/myo1/CMEM/WEEKLY/AVE';
var     = ['P_i';'N3n';'O2o'];
readRUN_doEOFs(gridVar,inDir,outDir,var)


%% eof puntuali con rea CMEMS

%gridVar = [LOC '3D_VAR/GRID/ALL/BFM_grid.nc'];
%outDir  = [LOC '3D_VAR/EOF/TEST_WEEKLY'];
%inDir   = '/myo1/CMEMS/WEEKLY/AVE';
%readRUN_doEOFs(gridVar,inDir,outDir)

%% generazione degli EOF SUB-ALL
%coast on fluxus

 inDIRvar = [LOC 'VAR_MOD/CCI1km_CMEMScoast/'     ];
 inDIReof = [LOC '3D_VAR/EOF/SUB24/'  ];
% inDIReof = [LOC '3D_VAR/EOF/SUB/KB/'     ];
%outDIReof = [LOC '3D_VAR/EOF/SUB_ALL/' ];
outDIReof = ['/myo1/Archive_opech_V5/EOF/SUB_ALL/']
% outDIReof = [LOC '3D_VAR/EOF/SUB_ALL/KB' ];

 coastv = 'v02';
 VV='V5_mask';
 
generateSUB_ALL_Eofs(inDIRvar,inDIReof, outDIReof,VV,coastv); 


% ATTENZIONE A:
%uso interno agli script di V4 (in readRUN_doEOFs.m ad esempio)
%numero di livelli imposto negli script (in generateSUB_ALL_Eofs.m ad esempio)
%quante EOF vogliamo?


