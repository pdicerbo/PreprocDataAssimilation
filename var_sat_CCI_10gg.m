% function var_sat_CCI_10gg(limstd,dayF,W,dirdata
% calculate surface variance for satellite 
% on 10 day mean (calculated extra script on the V4 model grid)
function var_sat_CCI_10gg(limstd,dayF,W,dirdata,dirvarsat)

%load(climfile)


jpi=W.jpi;
jpj=W.jpj;

%M=ncread(W.maskfile);
%lon=squeeze(M.nav_lon(1,:));
%lat=squeeze(M.nav_lat(:,1));
lon=W.lon24;
lat=W.lat24;

LL=dir([dirdata '/*.nc' ]);
numd=size(LL,1);


valCax=[.4 .4 .9 1 .6 .2 .1 .1 .1 .3 .5 .6];
Month{1}='January';
Month{2}='February';
Month{3}='March';
Month{4}='April';
Month{5}='May';
Month{6}='June';
Month{7}='July';
Month{8}='August';
Month{9}='September';
Month{10}='October';
Month{11}='November';
Month{12}='December';


%%
figure(1)
set(1,'Visible','off','Renderer','zbuffer')
chl2DT=zeros(jpj,jpi,12);
chl2D2T=zeros(jpj,jpi,12);
np=zeros(jpj,jpi,12);
countMonth=zeros(12,1);


for ii=1:numd
    dataMM=LL(ii).name(5:6);
    disp(['Month ' dataMM ' --- day index ' num2str(ii)])
    for mm=1:12
        if strcmp(num2str(mm,'%02d'),dataMM)
            countMonth(mm)=countMonth(mm)+1;
            npp=zeros(jpj,jpi);
            

                dataC=LL(ii).name(1:8);
                disp(['Reading data ' dataC ]);
                MM=ncread([dirdata '/' LL(ii).name]);
                chl2D=MM.CHL;
                disp(size(chl2D))
                chl2D(isnan(chl2D))=0.;
            

            app=zeros(jpj,jpi);
            app(chl2D~=0)=1;
            np(:,:,mm)=np(:,:,mm)+app;
            chl2DT(:,:,mm)=chl2DT(:,:,mm)+chl2D;
            chl2D2T(:,:,mm)=chl2D2T(:,:,mm)+chl2D.*chl2D;
        end
    end
end

np(np==0)=NaN;
chl2Dm=chl2DT./np;
chl2D2m=chl2D2T./np;
var2D=chl2D2m-chl2Dm.*chl2Dm;

%% Save varsat


for mm=1:12
    filevar=[dirvarsat 'var2Dsat.CCI.F' num2str(dayF) '.' ...
        num2str(limstd) '.'  num2str(mm,'%02d') '.nc'];
    S.variance.value=squeeze(var2D(:,:,mm));
    S.lon.value=lon;
    S.lat.value=lat;
    S.DIMS.lon=jpi;
    S.DIMS.lat=jpj;
    ncwrite(filevar,S)

    filenp=[dirvarsat 'np.CCI.F' num2str(dayF) '.' ...
        num2str(limstd) '.' num2str(mm,'%02d') '.nc'];
    T.np.value=squeeze(np(:,:,mm));
    T.DIMS.lon=jpi;
    T.DIMS.lat=jpj;
    ncwrite(filenp,T)
    
    varFig=(squeeze(var2D(:,:,mm)).^0.5);
    varFig(varFig==0)=NaN;
    pcolor(varFig);shading flat; colorbar
    caxis([0 valCax(mm)])
    title(Month(mm))
    nomeFig=[dirvarsat, 'devst.sat.CCI.F' num2str(dayF) '.' ...
        num2str(limstd) '.' num2str(mm,'%02d') '.jpg'];
    print(gcf,'-djpeg',nomeFig)
end


