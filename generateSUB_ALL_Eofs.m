function generateSUB_ALL_Eofs(inDIRvar,inDIReof,outDir,VV,coastv)

maskfunc=str2func([VV]);

submaskfile=['/myo1/Archive_opech_' VV(1:2) '/masks/sub16mask24.nc'];

MaskVersion=maskfunc(); 
%submaskfile=MaskVersion.submaskfiledat;

jpi = MaskVersion.jpi;
jpj = MaskVersion.jpj;
jpk = MaskVersion.jpk;
%M=ncread(MaskVersion.maskfile,'tmask');

if isempty(coastv)
    varfiletxt = 'var2D'; 
else
    varfiletxt = 'var2Dcoast';
end

eoffiletxt = 'eof';


%% setup submask
submask = true(jpj,jpi,16);
subtot  = zeros(jpj, jpi);

P = Polygonal_Med_SubBasin_Def_V2;
jjm = 1;
for jmask=1:16
    sb = Med_subBasin_switch(jmask,P);
    disp(sb)
    %disp(jjm)
    SS = ncread(submaskfile, sb);
    submask(:,:,jmask)=logical( SS.(sb) );
    subtot=subtot+submask(:,:,jmask)*jjm;
    if strcmp(sb,'adr1')
       continue
    else
       jjm = jjm+1;
    end
end

%% 


VAR=ncread([inDIRvar varfiletxt '.01' coastv '.nc']);
points = ~isnan(VAR.variance);
%points = (subtot~=0);

[I,J]=find(points'); %column arrays
np=numel(I);


month = 1;
eofFile  =[inDIReof 'eof.15.' num2str(month,'%02d') '.nc'];
N   = ncread(  eofFile);
eva = N.eva; evc = N.evc ; 
neof = size(eva,1);
nlevEOF = size(evc,2);

%%

for month=1:12

    eofFile  =[inDIReof 'eof.15.' num2str(month,'%02d') '.nc'];
    var2DFile=[inDIRvar varfiletxt '.' num2str(month,'%02d') coastv '.nc'] ;
    disp(eofFile)
    disp(var2DFile)

    N   = ncread(  eofFile);
    eva = N.eva ; evc = N.evc;
    VAR = ncread(var2DFile);  var2D=VAR.variance;

    evaNEW = zeros(neof,np);
    evcNEW = zeros(neof,jpk,np);


    for ip=1:np
        ii = I(ip);
        jj = J(ip);
        sub    = subtot(jj, ii);
        lambda = eva(:,  sub).^2;
        e      = evc(:,1,sub);
        sigma2 = lambda(1:neof)' * (e(1:neof).^2) ;
        alpha  = var2D(jj,ii)/sigma2;

        %    lambdaN = lambda* alpha^.5 ;
        %    eN      = e     * alpha^.25;

        evaNEW(:,  ip)   = eva(1:neof,  sub) * alpha^.25;
        evcNEW(:,1:nlevEOF,ip)   = evc(1:neof,:,sub) * alpha^.25;
        %evaNEW(:,  ip)   = eva(1:4,  sub) * alpha^.25;
        %evcNEW(:,:,ip)   = evc(1:4,:,sub) * alpha^.25;
        %if ii==605 && jj==44
        %    appl = lambda;
        %    appe = e;
        %    sigapp(neof) = sigma2;
        %    disp(sigma2)np
        %    disp(alpha)
        %    disp(var2D(jj,ii))
        %    appnl = evaNEW(:,ip).^2;
        %    appne = evcNEW(:,1,ip);
        %end

    end

%%end
    
    fileout =  [outDir eoffiletxt '.' num2str(np) '.' num2str(month,'%02d') '.' coastv '.nc'];
    disp(fileout)
    
    S.DIMS.nreg  = np ;
    S.DIMS.nlev  = jpk ;
    S.DIMS.neof  = neof ;

    S.eva.value = evaNEW ;
    S.evc.value = evcNEW ;
    S.eva.type = 'FLOAT';
    S.evc.type = 'FLOAT';

    ncwrite(fileout,S);
    clear('S')

end


