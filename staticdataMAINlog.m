

LOC='/usr2/home/matlab/OPA-DA/V0-static-data/';

% bisogna avere le due griglie, fornite da fortran90
% per SUB e SUB_ALL




%% creo i var2D.12.nc per il misfit
gridFile = [LOC '3D_VAR/GRID/SUB_ALL/BFM_grid.nc']; 
outDIR   = [LOC 'MISFIT_LOG/']; 

eofSUB_VAR_LOG(gridFile,outDIR) ; 


%% creo gli var2D.


gridFile = [LOC '3D_VAR/GRID/SUB/BFM_grid.nc']; 
outDIR   = [LOC '3D_VAR/EOF/' ];
readRUN_doEOFs_GRID_LOG(gridFile,outDIR); % scrive in SUB_LOG/eof.9.01.nc


%% generazione degli EOF SUB-ALL

 inDIRvar = [LOC 'MISFIT_LOG/'             ];
 inDIReof = [LOC '3D_VAR/EOF/SUB_LOG/'     ];
outDIReof = [LOC '3D_VAR/EOF/SUB_ALL_LOG/' ];

generateSUB_ALL_Eofs(inDIRvar,inDIReof, outDIReof,'V0'); 

