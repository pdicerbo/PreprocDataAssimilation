%function diff_CCI_dec(dirOutput,dirMod,dirSat,VV)
% evaluate difference variance of CCI and decadal model results
% dirMod,dirSat are directories of input files
function diff_CCINEW_RACMEMS(dirVE,dirMod,dirSat,VV)

%% Reading mask

jpiN=VV.jpi;
jpjN=VV.jpj;

MS=ncread(VV.maskfile,'nav_lon','nav_lat');


%% Reading data and evaluating variance
sattxt='_d-OC_CNR-L3-CHL-MedOC4AD4_SAM_1KM-MED-REP-v02.nc';

LM=dir([dirMod '/chl.*.nc']);
numF=numel(LM);
MM=ncread([dirMod '/' LM(1).name]);
jpi = size(MM.lchlm,2);
jpj = size(MM.lchlm,1);

S.DIMS.lon=jpiN;
S.DIMS.lat=jpjN;

S.lon.value=squeeze(MS.nav_lon(1,:));
S.lat.value=squeeze(MS.nav_lat(:,1));

%lonNew = squeeze(MS.nav_lon(1,:));
%latNew = squeeze(MS.nav_lat(:,1));

varDiff=zeros(jpj,jpi,12);

for im=1:12
    disp(['month ' num2str(im)])
    nump=zeros(jpj,jpi);
    sumDiff2=zeros(jpj,jpi);
    for ii=1:numF
        date=LM(ii).name(5:12);
        if strcmp(num2str(im,'%02d'),date(5:6))
            disp(date)
            dateS=datestr(datenum(date,'yyyymmdd')-3,'yyyymmdd');
            fileSat=[dirSat '/' dateS sattxt];
            if exist(fileSat,'file')                
                disp(' ... reading sat and model file')
                SS=ncread(fileSat);
                MM=ncread([dirMod '/' LM(ii).name],'lchlm');
                diffCHL=SS.CHL-MM.lchlm;
                diffCHL(isnan(diffCHL))=0;
                nump(diffCHL~=0)=nump(diffCHL~=0)+1;
                sumDiff2=sumDiff2+diffCHL.*diffCHL;
            end
        end
    end
    varDiff(:,:,im)=sumDiff2./nump;

    fileOut=[dirVE '/varErr.' num2str(im,'%02d') '.nc'];
    S.variance.value=squeeze(varDiff(:,:,im));
    disp(['writing ' fileOut])
    ncwrite(fileOut,S)
end

%%


