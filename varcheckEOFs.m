% readRUN_doEOFs
% read RUN BFM and prepare EOFs



function varcheckEOFs(gridVar,inDir,var)

%% EOF parameter
LISTneof = [26];
aa=size(LISTneof);
NNeof = aa(1);

%%
% dati from RA CARBO CMEMS
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
for im=1:1%2
    E_2nc = zeros(NNeof,NREGIONI,iz450,max(LISTneof),numvar);
    Taunc = zeros(NNeof,NREGIONI,max(LISTneof),numvar);
    profiliREG = zeros(iz450, nTimes, NREGIONI,numvar);
    ii=0;
    for it=1:nTimes
        filename = infile(it).name;        
        imf = str2num(filename(9:10));
        if imf==im
            ii=ii+1;
            disp(['reading ...' filename]);
            %STR   = ncread([inDir '/' filename],'P1i','P2i','P3i','P4i');
            STR = ncread([inDir '/' filename],'P_i','N3n','O2o');
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
                          profiliREG(1:iz450,ii,ireg,iv)=TheProfile(1:iz450); % profilo medio su tutta la regione
                       else
                          profiliREG(1:iz450,ii,ireg,iv)=0;
                       end %if isnan
                    end %if any
                end % for ireg
            end % ivar
        end %if im
    end %for it
    profiliREG = profiliREG(:,1:ii,:,:);
    disp(im)
    disp(size(profiliREG))
    TT = (squeeze(profiliREG(:,:,1,:)));
    disp(size(TT))
    VT = std(TT').^2
    disp(size(VT))


end %im

