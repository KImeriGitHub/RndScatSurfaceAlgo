function RMat_new = getInflation(r_inc_Data, rExt, ZetaExt, NFourier, k)

eg=0.57721566490153286060651209008240243104215933593992;
twopi=2*pi;
fourpi=4*pi;

r_init = r_inc_Data(1);

diffeps=0.00001; % < 1e-6 might yield machine precision errors.

DomExt = shape.Ellipse(rExt, rExt, NFourier)+ZetaExt;
[~,~,~,R0zetazeta, ~]= getNeumannFct(DomExt, [0;0], [0;0], k);

Dom_Null = shape.Ellipse(r_init, r_init, NFourier)+[0;0];
theta_ip=linspace(0,twopi,NFourier+1);
theta_ip(end)=[];
tPts_ip = linspace(-pi,pi,NFourier+1);
tPts_ip = tPts_ip(1:end-1);


%PRE -Computation
sigmaMat_init=repmat(circshift(interp1(Dom_Null.theta,Dom_Null.sigma,theta_ip),NFourier/2,2).',1,NFourier);
CosMat=cos((tPts_ip)-(tPts_ip).');
sqrt2CosMat=sqrt(2-2*CosMat);
logMatp4pi=log(1-CosMat)/fourpi;
for j=1:size(logMatp4pi,1)
     logMatp4pi(j,j)=(log(1-cos(theta_ip(2)))-4)/4/pi;%(logMatp4pi(j,mod(j,NumDiscPts)+1)-4)/4/pi;
end % For testing: (logMatp4pi(1,:)*(sigmaMat_init(:,1)/r_init)) == -log(2)/2
logp4pi_FCoeffs = [-log(2)/2, -1./(1:NFourier-1)]./twopi;

%% Main Loop

% CASE IF r_init IS NOT VERY VERY SMALL
% [~,~,~,~, RMat_init]= getNeumannFct(Dom_init, pts_normed*r_init, [nan;nan], k);
% RMat_init=RMat_init(NFourier+1:end,:);
% RMat_init=circshift(RMat_init,NFourier/2);
% R_Mat=RMat_init+(-1i/4*besselh(0,k*r_init*sqrt2CosMat)-log(1-CosMat)/fourpi);
% for j=1:size(R_Mat,1)
%      R_Mat(j,j)=(RMat_init(j,j)+(2*eg-1i*pi+log(2)+2*log(k*r_init/2))/4/pi);
% end
% [R_CosConst, R_SinConst] = getFourierConst(R_Mat); %Length of Coeffs is NumDiscPts/2
% R_CosConst=[R_CosConst,zeros(NFourier,NFourier/2)];
% R_SinConst=[R_SinConst(:,1:NFourier/2),zeros(NFourier,NFourier/2)];

R_CosConst=[(log(k*r_init)/fourpi + eg/fourpi - 1i/8 + R0zetazeta/2)*ones(NFourier,1), zeros(NFourier,NFourier-1)];
R_SinConst=zeros(NFourier,NFourier);

for i=1:length(r_inc_Data)-1
    %disp(['i = ',num2str(i), ' of ', num2str(length(r_inc_Data)-1)]);
    r=r_inc_Data(i);
    R=r_inc_Data(i+1);
    
    thetaPts=[cos(linspace(-pi,pi,NFourier+1));sin(linspace(-pi,pi,NFourier+1))];
    thetaPts(:,end)=[];
    
    [~,~,~,R0zRxR, ~]= getNeumannFct(DomExt, thetaPts*R, [thetaPts*R,thetaPts*(R+diffeps)], k);
    delxR_R0zRxR=(R0zRxR(NFourier+1:end,:)-R0zRxR(1:NFourier,:))/diffeps;
    [delxR_R0zRxR_CosConst, delxR_R0zRxR_SinConst] = getFourierConst(delxR_R0zRxR.');

    [~,~,~,R0zRxr, ~]= getNeumannFct(DomExt, [thetaPts*R,thetaPts*(R+diffeps)], [thetaPts*r,thetaPts*(r+diffeps)], k);   
    delxr_R0zRxr=(R0zRxr(NFourier+1:end, 1:NFourier)-R0zRxr(1:NFourier,1:NFourier))/diffeps;
    delzRdelxr_R0zRxr=((R0zRxr(NFourier+1:end,NFourier+1:end)+R0zRxr(1:NFourier,1:NFourier)    )...
                      -(R0zRxr(1:NFourier,NFourier+1:end)    +R0zRxr(NFourier+1:end,1:NFourier)))/diffeps^2; % This thing is symmetric (good enough)
    %delxr_R0zRxr == delxr_R0zrxr+ (R-r)*delzRdelxr_R0zRxr ? Confirmed up to e-5  colormap('hot'); imagesc(real(delxr_R0zRxr -( delxr_R0zrxr+ (R-r)*delzRdelxr_R0zRxr))); colorbar;
    %Could use here a Taylor series to decrease number of computations
    
    R_Mat = R_CosConst*[ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]...
           +R_SinConst*[zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
       
    sqrtRrCosMat=sqrt(R^2+r^2-2*R*r*CosMat);
    
    logRrp2pi_FCoeffs = 2*[log(R)/2/pi-log(2*R*r)/4/pi, -(r/R).^(1:NFourier-1)./(1:NFourier-1)./(twopi)];
    logRrp4pi_FCoeffs = [log(R/r/2)/fourpi, -(r/R).^(1:NFourier-1)./(1:NFourier-1)./(twopi)];
      
    log_CosConst=cos(tPts_ip.'*(0:NFourier-1)).*repmat(logp4pi_FCoeffs,NFourier,1);
    logRr_CosConst=cos(tPts_ip.'*(0:NFourier-1)).*repmat(logRrp4pi_FCoeffs,NFourier,1);
    log_SinConst=sin(tPts_ip.'*(0:NFourier-1)).*repmat(logp4pi_FCoeffs,NFourier,1);
    logRr_SinConst=sin(tPts_ip.'*(0:NFourier-1)).*repmat(logRrp4pi_FCoeffs,NFourier,1);
    IntegConst=repmat([twopi,pi*ones(1,NFourier-1)],NFourier,1);
    
    DelRRGamma_NonSing=1i*k/4/sqrt(2)*(sqrt(1-CosMat).*besselh(1,k*R*sqrt2CosMat))-1/(2*R*twopi);
    DelrrGamma_NonSing=1i*k/4/sqrt(2)*(sqrt(1-CosMat).*besselh(1,k*r*sqrt2CosMat))-1/(2*r*twopi);
    DelRrGamma_NonSing=1i*k/4*(r-R*CosMat)./sqrtRrCosMat.*besselh(1,k*sqrtRrCosMat)-(r-R*CosMat)./(R^2+r^2-2*R*r*CosMat)./twopi;
    for j=1:size(DelRRGamma_NonSing,1)
        DelRRGamma_NonSing(j,j)=0;
        DelrrGamma_NonSing(j,j)=0;
    end
    [DelRRGns_CosConst, DelRRGns_SinConst] = getFourierConst(DelRRGamma_NonSing);
%     [DelrrGns_CosConst, DelrrGns_SinConst] = getFourierConst(DelrrGamma_NonSing);
%     [DelRrGns_CosConst, DelRrGns_SinConst] = getFourierConst(DelRrGamma_NonSing);
    
    tempNonSingGammaR=(-1i/4*besselh(0,k*R*sqrt2CosMat)-logMatp4pi);
    tempNonSingGammar=(-1i/4*besselh(0,k*r*sqrt2CosMat)-logMatp4pi);
    for j=1:size(tempNonSingGammaR,1)
        tempNonSingGammaR(j,j)=((2*eg-1i*pi+log(2)+2*log(k*R/2))/4/pi);
        tempNonSingGammar(j,j)=((2*eg-1i*pi+log(2)+2*log(k*r/2))/4/pi);
    end
    
    %% GETTING delyR_NrxRyR
    DDG2_FCoeffs=zeros(1,NFourier/2);
    DDG12_FCoeffs=zeros(1,NFourier/2);
    for j=1:NFourier
        DDG2_FCoeffs(j) = 2*0.5*integral(@(t)(-k^2).*(R-r*cos(t)).*(r-R*cos(t))./(twopi*(R^2+r^2-2*R*r*cos(t)))...
            .*log(R^2+r^2-2*R*r*cos(t)).*cos((j-1)*t),0,pi,'Waypoints',(1/2+(0:max(1,(j-1))))*pi/max(1,(j-1)));
        DDG12_FCoeffs(j) = 2*0.5*integral(@(t)(-k^2).*(-2*R*r+(R^2+r^2).*cos(t))./(twopi*(R^2+r^2-2*R*r*cos(t)))...
            .*log(k*sqrt(R^2+r^2-2*R*r*cos(t))).*cos((j-1)*t),0,pi,'Waypoints',(1/2+(0:max(1,(j-1))))*pi/max(1,(j-1)));
    end
            
    SingSingConv  = 2*1/4/pi*r/R*(r^2-R^2*CosMat)./(R^4+r^4-2*R^2*r^2*CosMat);
%      SingSingConv2 = (cos((tPts_ip).'*(0:NFourier-1))) * (repmat(DDG2_FCoeffs.'.*logRrp2pi_FCoeffs.',1,NFourier).*(cos((0:NFourier-1).'*(tPts_ip))))...
%         +(sin((tPts_ip).'*(0:NFourier-1))) * (repmat(DDG2_FCoeffs.'.*logRrp2pi_FCoeffs.',1,NFourier).*(sin((0:NFourier-1).'*(tPts_ip))));
    SingSingConv2_CosConst=cos(tPts_ip.'*(0:NFourier-1)).*repmat(DDG2_FCoeffs.*logRrp2pi_FCoeffs,NFourier,1);
    SingSingConv2_SinConst=sin(tPts_ip.'*(0:NFourier-1)).*repmat(DDG2_FCoeffs.*logRrp2pi_FCoeffs,NFourier,1);
    
%     SingSingConv3 = (cos((theta.'-pi)*(0:NFourier/2-1))) * (repmat(DDG12_FCoeffs.'.*logRrp2pi_FCoeffs.',1,NumDiscPts).*(cos((0:NFourier/2-1).'*(theta-pi))))...
%         +(sin((theta.'-pi)*(0:NFourier/2-1))) * (repmat(DDG12_FCoeffs.'.*logRrp2pi_FCoeffs.',1,NumDiscPts).*(sin((0:NFourier/2-1).'*(theta-pi))));
    SingSingConv3_CosConst=cos(tPts_ip.'*(0:NFourier-1)).*repmat(DDG12_FCoeffs.*logRrp2pi_FCoeffs,NFourier,1);
    SingSingConv3_SinConst=sin(tPts_ip.'*(0:NFourier-1)).*repmat(DDG12_FCoeffs.*logRrp2pi_FCoeffs,NFourier,1);
    
    
    smoothDelDelGammaTerm=-(1/4)*1i.*k.*(...
        k.*(-R + r.*CosMat).*(r- R.*CosMat)./(sqrtRrCosMat.^2).*(besselh(0,k.*sqrtRrCosMat)-2i./pi.*log(sqrtRrCosMat)) ...
    - ((-2.*r.*R + (r.^2 + R.^2).*CosMat)./(sqrtRrCosMat.^3).*(besselh(1,k.*sqrtRrCosMat)+2i./pi./(k.*sqrtRrCosMat)-(2i*k.*sqrtRrCosMat.*log(k.*sqrtRrCosMat))./twopi)));
     
    NSmoo_CosConst = 2*R_CosConst;
    NSmoo_SinConst = 2*R_SinConst;
    Inc3_CosConst=NSmoo_CosConst.*repmat([0,(r/R).^(0:NFourier-2).*(1:NFourier-1)/2/R^2],NFourier,1); %critical!
    Inc3_SinConst=NSmoo_SinConst.*repmat([0,(r/R).^(0:NFourier-2).*(1:NFourier-1)/2/R^2],NFourier,1); %critical!

    Inc4_CosConst=NSmoo_CosConst.*repmat(DDG2_FCoeffs,NFourier,1);
    Inc4_SinConst=NSmoo_SinConst.*repmat(DDG2_FCoeffs,NFourier,1);

    Inc5_CosConst=NSmoo_CosConst.*repmat(DDG12_FCoeffs,NFourier,1);
    Inc5_SinConst=NSmoo_SinConst.*repmat(DDG12_FCoeffs,NFourier,1);

    SmooNTerm=2*(R_Mat);
    Inc6  =(smoothDelDelGammaTerm.'*((SmooNTerm+2*(log(R^2+r^2-2*R*r*CosMat)/4/pi-log(2*R*r)/4/pi)).'.*(sigmaMat_init*r/r_init))).'; %This might blow up?!
    Inc6R =(delzRdelxr_R0zRxr.'*((SmooNTerm+2*(log(R^2+r^2-2*R*r*CosMat)/4/pi-log(2*R*r)/4/pi)).'.*(sigmaMat_init*r/r_init))).'; %This might blow up?!
    [Inc6_CosConst, Inc6_SinConst] = getFourierConst(-Inc6-Inc6R);
    
    DelInc_CosConst=-r*(SingSingConv2_CosConst+SingSingConv3_CosConst+Inc3_CosConst+Inc4_CosConst+Inc5_CosConst)...
                    +[Inc6_CosConst + DelRRGns_CosConst+delxR_R0zRxR_CosConst,zeros(NFourier,NFourier/2)]; % ZEROS extension
    DelInc_SinConst=-r*(SingSingConv2_SinConst+SingSingConv3_SinConst+Inc3_SinConst+Inc4_SinConst+Inc5_SinConst)...
                    +[Inc6_SinConst + DelRRGns_SinConst+delxR_R0zRxR_SinConst,zeros(NFourier,NFourier/2)]; % ZEROS extension
    
    %% GETTING NRxRyR
    %temp =-2*(DelIncNonSing.'  *  ((R_Mat+logMatp4pi).'.*(sigmaMat_init*R/r_init))).';
    temp=(-2*R)*(R_CosConst*(DelInc_CosConst.*IntegConst).'+R_SinConst*(DelInc_SinConst.*IntegConst).'+...
        +log_CosConst*(DelInc_CosConst.*IntegConst).'+log_SinConst*(DelInc_SinConst.*IntegConst).');
    [temp_CosConst,temp_SinConst]=getFourierConst(temp);
    
    temp2=-2*((1i*k/4*(r-R*CosMat)./sqrtRrCosMat.*besselh(1,k*sqrtRrCosMat)-(r-R*CosMat)./(R^2+r^2-2*R*r*CosMat)./twopi).'  *  ...
        ((R_Mat+log((R^2+r^2)/(2*R*r)-CosMat)./fourpi).'.*(sigmaMat_init*r/r_init))).';
    temp2R=-2*(delxr_R0zRxr.' * ((R_Mat+log((R^2+r^2)/(2*R*r)-CosMat)./fourpi).'.*(sigmaMat_init*r/r_init))).';
    [temp2_CosConst, temp2_SinConst] = getFourierConst(temp2+temp2R);
    [GammaNS_CosConst,  GammaNS_SinConst]  = getFourierConst(tempNonSingGammaR + R0zRxR(1:NFourier,:).');
    
    R_CosConst_new = 0.5*([GammaNS_CosConst,zeros(NFourier,NFourier/2)]...
        +R_CosConst.*repmat([0,1+(r/R).^(1:NFourier-1)-(r/R).^(2*(1:NFourier-1))],NFourier,1)...
        +(log(2)/fourpi).*[ones(NFourier,1),zeros(NFourier,NFourier-1)]...
        +[temp_CosConst+temp2_CosConst,zeros(NFourier,NFourier/2)]); % ZEROS extension

    R_SinConst_new = 0.5*([GammaNS_SinConst,zeros(NFourier,NFourier/2)]...
        +R_SinConst.*repmat([0,1+(r/R).^(1:NFourier-1)-(r/R).^(2*(1:NFourier-1))],NFourier,1)...
        +[temp_SinConst + temp2_SinConst,zeros(NFourier,NFourier/2)]);  % ZEROS extension
    
    RMat_new =  R_CosConst_new * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]...
              + R_SinConst_new * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
    
    %% Post Loop Computations
    assert(norm(RMat_new)<1e10)
          
    R_CosConst=R_CosConst_new;
    R_SinConst=R_SinConst_new;
end
end

