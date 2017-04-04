% function eofSUB_VAR(gridFile,outDir)

% calculate surface variance (here from a KB archive)
% * modulated with the surface variance (for each grid point)


% and interpolates it on gridFile

function eofSUB_VAR(gridFile,outDir)
V = KB_mask;

jpi=V.jpi;
jpj=V.jpj;

% ----------------- nearest interpolation data --------------
filemask1=V.maskfile;
M1=ncread(filemask1,'nav_lon','nav_lat','nav_lev','tmask');


filemask2=gridFile;
M2=ncread(filemask2,'lon','lat','tmsk');

tmask1 = logical(squeeze(M1.tmask(1,:,:))); 
tmask2 = logical(squeeze(M2.tmsk( 1,:,:))); 

IND = getIndForFastinterp2(M1.nav_lon, M1.nav_lat, tmask1, M2.lon, M2.lat, tmask2); 
% ----------------------------------------------------------




dirdata='/myo1/Archive_DA/CH3_R29_KB';
ndmese=[31 28 31 30 31 30 31 30 31 31 30 31];

for mm=1:12
    chl2DT =zeros(jpj,jpi);
    chl2D2T=zeros(jpj,jpi);
    np=0;
    for yyyy=1999:2005
        for dd=1:ndmese(mm)
            yyyytxt= num2str(yyyy);
            mmtxt  = num2str(mm,'%02d'); 
            ddtxt  = num2str(dd,'%02d'); 
            filename=[dirdata '/ch3.' yyyytxt mmtxt ddtxt '.nc'];
           
            if exist(filename,'file')
                np=np+1;
                disp(['Reading ' filename]);
                CHL=ncread(filename);
                chl2D=squeeze(CHL.ch3(1,:,:));
                chl2D(chl2D>=9.999E19)=NaN;
                chl2DT = chl2DT+chl2D;
                chl2D2T= chl2D2T+chl2D.*chl2D;
            end
        end
    end
    chl2Dm  = chl2DT/np;
    chl2D2m = chl2D2T/np;
    var2D   = chl2D2m-chl2Dm.*chl2Dm;
    

    % from KB to V0_DA mesh ---------------------------------
    var2D72 = interp2_WetPoint_fast(var2D,tmask1,tmask2,IND);
    % -------------------------------------------------------
    
%     %from 843 to 872 mesh
%     var2D72=grid43togrid72(var2D);
%     var2D72(var2D72==0)=NaN;


    filevar=[outDir 'var2D.' num2str(mm,'%02d') '.nc'];
    S.variance.value=var2D72;
    S.DIMS.lon=362;
    S.DIMS.lat=128;
    ncwrite(filevar,S)
end


