% readRUN_doEOFs_GRID_LOG(gridFile,outDir)
% read RUN BFM on grid43 and prepare EOFs on grid 72 on logaritmic value



function readRUN_doEOFs_GRID(gridFile,outDir)


% --------- i dati sono su mesh KB ----
dirdata='/myo1/Archive_DA/CH3_R29_KB' ;
V = KB_mask;
%-------------------------------------

% ----------------- nearest interpolation data --------------
filemask1=V.maskfile;
M1=ncread(filemask1,'nav_lon','nav_lat','nav_lev','tmask','gdept');% gdept altezza del centro cella

M2 = ncread(gridFile);

nx=size(M2.lon,2);
ny=size(M2.lon,1);

tmask1 = logical(squeeze(M1.tmask(1,:,:)));
tmask2 = logical(squeeze(M2.tmsk( 1,:,:)));

IND = getIndForFastinterp2(M1.nav_lon, M1.nav_lat, tmask1, M2.lon, M2.lat, tmask2);

z200_M1 =getDepthIndex(M1.nav_lev,200) + 1 ; % profondita fino a 208.4metri nella griglia 43
z200_M2 =getDepthIndex(M2.dep    ,200) + 1 ;  % il 26 livello e' a quota 209.5

%%% considerare solo i punti oltre 200 metri 14^lev (208.448)
%%% fare eof fino a 200 m
% profondita' delle 2 griglie fino a 200 metri
depth_1 = M1.gdept(1:z200_M1);
depth_2 = M2.dep(  1:z200_M2);

nlev_2 = numel(M2.dep);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% verifico quante regioni ci sono
qualiregs = unique(M2.regs(:));
NREGIONI  = length(qualiregs);
disp(['number of regions =' num2str(NREGIONI)]); %disp('regions: ');disp(qualiregs);


ndmese=[31 28 31 30 31 30 31 30 31 31 30 31];


% conto quanti profili ci sono per profondita' fino a 209.5 metri (26lev)
A01   = M2.tmsk(1:z200_M2,:,:); %
SUP   = squeeze(sum(A01,1)); % numero di celle da 0 a z200_M2

REGSnoterra = M2.regs;
REGSnoterra(SUP<z200_M2) = 0;
vREGSnoterra       = REGSnoterra(:); % vettore con le regioni solo acqua oltre 200 metri

nTimes=1000;
%%%%%%%LETTURA FILE RUN AND DO EOFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for mm=1:12
    ii         =  0;
    profiliREG = zeros(z200_M1, nTimes, NREGIONI) ;

    for yyyy=1999:2005
        for dd=1:ndmese(mm)
            yyyytxt  = num2str(yyyy);
            mmtxt    = num2str(mm,'%02d');
            ddtxt    = num2str(dd,'%02d');
            filename = [dirdata '/ch3.' yyyytxt mmtxt ddtxt '.nc'];
            if exist(filename,'file')
                disp(['reading ...' filename(end-14:end)]);
                STR   = ncread(filename);
                chl3D = STR.ch3; chl3D(chl3D >= 9.999E19)=NaN;
                chl3D=log(chl3D);

                % trasferisco la map43 su map72

                Mv=zeros(z200_M1,ny*nx);
                for ip=1:z200_M1
                    map_1=squeeze(chl3D(ip,:,:));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % map72 = grid43togrid72(map43); % ricostruisco il campo 3D su grid72
                    map_2 = interp2_WetPoint_fast(map_1,tmask1,tmask2,IND);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Mv(ip,:)=map_2(:);

                end %for ip (levs)
                % attenzione mantengo i livelli di grid43 ... interpolero'
                % poi i profili degli autovettori
                %                Mv(1:z200_M1, 1:nx72*ny72)=reshape(M,[z200_M1,nx72*ny72]);% vettorizzo la matrice
                ii=ii+1 ;
               % inan=0;

                % ------------- ciclo sulle regioni, DEVE essere veloce -----------------------------------------
                for ireg=1:NREGIONI,
                    quali = vREGSnoterra==ireg;
                    if any(quali),
                        Mv0 = Mv(:,quali); 
                        %%%% ATTENZIONE USO il nanmean perche' nel trasferire map43 in map72 alcuni punti
                        %%%% hanno profili diversamente lunghi

                        TheProfile=nanmean(Mv0,2) ;
                        if ~isnan(TheProfile)
                            profiliREG(1:z200_M1,ii,ireg)=TheProfile; % profilo medio su tutta la regione
                        else
                            profiliREG(1:z200_M1,ii,ireg)=0;
                            % inan=inan+1;
                        end %if isnan
                    end
                end %for ireg
                % ------------------------------------------------------------------------------------------------
               % disp(inan);
            end %if exist
        end %for dd
    end %for yyyy


profiliREG = profiliREG(:,1:ii,:); 

    % tante EOFs quanti sono i profili
    disp('calcolo anomalie')
    
    E_2nc = zeros(NREGIONI, nlev_2 , z200_M1 ) ; % ma sul secondo indice lo riempio fino a z200_M2
    Taunc = zeros(NREGIONI,          z200_M1 ) ;
    
    for ireg=1:NREGIONI
        disp(['profili ' num2str(ireg)])
        profili=squeeze(profiliREG(:,:,ireg));
        [n m]=size(profili); % n e' uguale a z200_M1

        Mprofili   = mean(profili,2); % profilo medio per tutti  casi di un mese
        ANOprofili = profili-repmat(Mprofili,1,m); % anomalie % NON SONO normalizzate


        [E psev V]=svd(ANOprofili,0);
        [n nE]=size(E); % nE number of EOFs
        
        E_2= interp1(depth_1, E, depth_2,'linear','extrap');
        E_2nc(ireg,1:z200_M2,:)=E_2;
        
        Taunc(ireg,:)=diag(psev/sqrt(m)); % eigenvalue diagional of Tau [np x np]
       
    end

    %% Save date in eof.NREGIONI.mm.nc file for 3dvar
   
    caseNREG=NREGIONI ;
    if NREGIONI>1000, caseNREG=100; end

    switch caseNREG
        case   1, eofFile=[outDir 'MED_LOG/eof.',num2str(NREGIONI),'.',num2str(mm,'%02d'),'.nc'];
        case   2, eofFile=[outDir 'ESO_LOG/eof.',num2str(NREGIONI),'.',num2str(mm,'%02d'),'.nc'];
        case   9, eofFile=[outDir 'SUB_LOG/eof.',num2str(NREGIONI),'.',num2str(mm,'%02d'),'.nc'];
        case 100, eofFile=[outDir 'ALL_LOG/eof.',num2str(NREGIONI),'.',num2str(mm,'%02d'),'.nc'];
        otherwise
            disp('Region type should be 1, 2 or 3')
    end

    disp(eofFile)
    % --------- scrittura su file -------------------------------------
    cdfid = mexcdf('CREATE',eofFile, 'CLOBBER');

    % Define dimensions.
    nregdimid  = mexcdf('DIMDEF', cdfid, 'nreg', NREGIONI);
    nlevdimid  = mexcdf('DIMDEF', cdfid, 'nlev', nlev_2  );
    neofdimid  = mexcdf('DIMDEF', cdfid, 'neof', nE      );

    % Define variables.
    evaid     = mexcdf('VARDEF', cdfid, 'eva',   'FLOAT', 2, [neofdimid           nregdimid] );
    evcid     = mexcdf('VARDEF', cdfid, 'evc',   'FLOAT', 3, [neofdimid nlevdimid nregdimid] );
    mexcdf('ENDEF', cdfid);

    % Store data
    mexcdf('VARPUT', cdfid,  evaid, [0 0  ], [nE         NREGIONI], Taunc);
    mexcdf('VARPUT', cdfid,  evcid, [0 0 0], [nE nlev_2  NREGIONI], E_2nc);
    mexcdf('SYNC', cdfid);

    % Close output  file.
    mexcdf('CLOSE', cdfid);

end  

