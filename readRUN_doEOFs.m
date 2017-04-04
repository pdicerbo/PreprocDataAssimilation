% readRUN_doEOFs
% read RUN BFM and prepare EOFs



function readRUN_doEOFs(gridVar,inDir,outDir,var)

%% EOF parameter
LISTneof = [26];
aa=size(LISTneof);
NNeof = aa(1);

%%
% dati from RA CMEMScoast
V = V4_mask;

%% Griglia modello
filemask1=V.maskfile;
M=ncread(filemask1,'nav_lon','nav_lat','nav_lev','tmask');

% considerare solo i punti con profondità maggiore di 200m
% fare eof fino a 200 m

% un livello in più rispetto ai 200m
iz200 = getDepthIndex(M.nav_lev,200) + 1 ;
iz450 = getDepthIndex(M.nav_lev,450) + 1 ;
nlev = size(M.nav_lev,1);

%% Griglia 3DVAR
MV = ncread(gridVar);

% numero regioni
qualiregs = unique(MV.regs(:));
NREGIONI  = length(qualiregs);
disp(['number of regions =' num2str(NREGIONI)]); %disp('regions: ');disp(qualiregs);

%% ???

% conto quanti profili ci sono per profondita' fino a 209.5 metri (26lev)
A01   = MV.tmsk(1:iz200,:,:); %
SUP   = squeeze(sum(A01,1)); % numero di celle da 0 a z200_M2

REGSnoterra = MV.regs;
REGSnoterra(SUP<iz200) = 0;

%% LETTURA FILE RUN AND DO EOFs
disp('Eof for variable(s):')
numvar = size(var,1);
for iv=1:numvar
    disp(var(iv,:))
end
infile = dir([inDir '/ave.*nc']);
nTimes = numel(infile);

%%
%disp(['Number of eof ' num2str(neof)])
for im=1:12
    E_2nc = zeros(NNeof,NREGIONI,iz200,max(LISTneof),numvar);
    Taunc = zeros(NNeof,NREGIONI,max(LISTneof),numvar);
    profiliREG = zeros(iz200, nTimes, NREGIONI,numvar);
    ii=0;
    for it=1:nTimes
        filename = infile(it).name;        
        imf = str2num(filename(9:10));
        if imf==im
            ii=ii+1;
            disp(['reading ...' filename]);
            %STR   = ncread([inDir '/' filename],'P1i','P2i','P3i','P4i');
            STR = ncread([inDir '/' filename],'P_l');
            %chl3D = STR.P_i+STR.P2i+STR.P3i+STR.P4i;
            for iv=1:numvar
                vv = var(iv,:);
                chl3D = STR.(vv);
                % ciclo sulle regioni
                for ireg=1:NREGIONI
                    quali = REGSnoterra==ireg;
                    if any(quali(:))
                       chlv_M = chl3D(:,quali);
                       TheProfile=nanmean(chlv_M,2) ;
                       if ~all(isnan(TheProfile))
                          profiliREG(1:iz200,ii,ireg,iv)=TheProfile(1:iz200); % profilo medio su tutta la regione
                       else
                          profiliREG(1:iz200,ii,ireg,iv)=0;
                       end %if isnan
                    end %if any
                end % for ireg
            end % ivar
        end %if im
    end %for it
    profiliREG = profiliREG(:,1:ii,:,:);

    for ireg=1:NREGIONI
        for iv=1:numvar
            disp(['Variable ' var(iv,:)])
            disp(['profili ' num2str(ireg)])
            profili=squeeze(profiliREG(:,:,ireg,iv));
            [n,m]=size(profili); % n e' uguale a iz200
            
            Mprofili   = mean(profili,2); % profilo medio per tutti  casi di un mese
            ANOprofili = profili-repmat(Mprofili,1,m); % anomalie % NON SONO normalizzate
            
            
            [E,psev,V]=svd(ANOprofili,0);
            [n,nE]=size(E); % nE number of EOFs
            disp('nE')
            disp(nE)
            
            eva = diag(psev/sqrt(m-1));
            evc = E(:,:);
            
            
            for iin=1:NNeof
                neof = LISTneof(iin);
                if neof~=nE
                    %calibrate variance
                    lambda = eva.^2;
                    TOTvariance = sum(lambda);
                    NEOFvariance = sum(lambda(1:neof));
                    alpha = TOTvariance/NEOFvariance;
                    evanew = eva(1:neof) * alpha^.5;
                else
                    evanew = eva;
                end
                
                E_2nc(iin,ireg,1:iz200,1:neof,iv)=evc(:,1:neof);
                Taunc(iin,ireg,1:neof,iv)=evanew; % eigenvalue diagional of Tau [np x np]
            end
        end
        
    end
    disp(sum(squeeze(E_2nc(1,1,1,:,1).^2)).*(squeeze(Taunc(1,1,:,1).^2)))
    any(E_2nc(:)>1.e+19)
    any(Taunc(:)>1.e+19)

    %% Save date in eof.I.mm.nc file for 3dvar

    for iv=1:numvar
        vv = var(iv,:);
        for iin=1:NNeof
            neof = LISTneof(iin);
            caseNREG=NREGIONI ;
            if NREGIONI>1000, caseNREG=100; end
            
            switch caseNREG
                case   1, eofFile=[outDir '/MED/eof.',num2str(NREGIONI),'.',num2str(im,'%02d'),'.nc'];
                case   2, eofFile=[outDir '/ESO/eof.',num2str(NREGIONI),'.',num2str(im,'%02d'),'.nc'];
                case   9, eofFile=[outDir '/SUB/eof.',num2str(NREGIONI),'.',num2str(im,'%02d'),'.nc'];
                %for float case  15, eofFile=[outDir '/SUB15_450/' vv '/NEOF',num2str(neof),'/eof.',num2str(NREGIONI),'.',num2str(im,'%02d'),'.nc'];
                case  15, eofFile=[outDir '/SUB/eof.',num2str(NREGIONI),'.',num2str(im,'%02d'),'.nc']; 
                case 100, eofFile=[outDir '/ALL/eof.',num2str(NREGIONI),'.',num2str(im,'%02d'),'.nc'];
                otherwise
                    disp('Region type should be 1, 2 or 3')
            end
            
            disp(eofFile)
            % --------- scrittura su file -------------------------------------
            cdfid = mexcdf('CREATE',eofFile, 'CLOBBER');
            
            % Define dimensions.
            nregdimid  = mexcdf('DIMDEF', cdfid, 'nreg', NREGIONI);
            nlevdimid  = mexcdf('DIMDEF', cdfid, 'nlev', iz200  );
            neofdimid  = mexcdf('DIMDEF', cdfid, 'neof', neof      );
            
            % Define variables.
            evaid     = mexcdf('VARDEF', cdfid, 'eva',   'FLOAT', 2, [neofdimid           nregdimid] );
            evcid     = mexcdf('VARDEF', cdfid, 'evc',   'FLOAT', 3, [neofdimid nlevdimid nregdimid] );
            mexcdf('ENDEF', cdfid);
            
            % Store data
            mexcdf('VARPUT', cdfid,  evaid, [0 0  ], [neof       NREGIONI], Taunc(iin,:,1:neof));
            mexcdf('VARPUT', cdfid,  evcid, [0 0 0], [neof iz200 NREGIONI], E_2nc(iin,:,:,1:neof));
            mexcdf('SYNC', cdfid);
            
            % Close output  file.
            mexcdf('CLOSE', cdfid);
        end
    end

end %im

