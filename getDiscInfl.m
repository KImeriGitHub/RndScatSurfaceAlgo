function [R_Cos, R_Sin]= getDiscInfl(R_Cos, R_Sin, RTens, r, R, k, dPts, DDG2_FCoeffs, DDG12_FCoeffs)
%PRE: R_Cos, R_Sin are each complex NFourier x NFourier matrices
%     RTens is a complex 5 x 5 x 5 x 5 x 7 tensor
%     CenterIdx is a (2 x NCenter) matrix of pos Integer
%     r, R, k are positive real numbers, R>r
%POST:
%DESC: Via the input from getOptCircles

NFourier  = size(R_Cos,1);
NZGridRow = size(R_Cos,1);
eg=0.57721566490153286060651209008240243104215933593992;
twopi=2*pi;
fourpi=4*pi;

tPts_ip = linspace(-pi,pi,NFourier+1);
tPts_ip = tPts_ip(1:end-1);
thetaPts = [cos(linspace(-pi,pi,NFourier+1));sin(linspace(-pi,pi,NFourier+1))];
thetaPts(:,end) = [];
theta_ip = linspace(0,twopi,NFourier+1);
theta_ip(end) = [];

%% PRE -Allocation
Inc6_CosConst=zeros(NFourier,NFourier);
Inc6_SinConst=zeros(NFourier,NFourier);

%% PRE -Computation
sigmaMat_r=r*twopi/NFourier;

CosMat=cos((tPts_ip)-(tPts_ip).');
sqrt2CosMat=sqrt(2-2*CosMat);
logMatp4pi=log(1-CosMat)/fourpi;
for j=1:size(logMatp4pi,1)
     logMatp4pi(j,j)=(log(1-cos(theta_ip(2)))-4)/4/pi;
end 
logp4pi_FCoeffs = [-log(2)/2, -1./(1:NFourier-1)]./twopi;

%nD interpolation
Rzz=RTens(:, :, :, :, 1); 
nRzz1=RTens(:, :, :, :, 2);  
nRzz2=RTens(:, :, :, :, 3); 
nnRz1z1=RTens(:, :, :, :, 4); 
nnRz1z2=RTens(:, :, :, :, 5); 
nnRz2z1=RTens(:, :, :, :, 6);  
nnRz2z2=RTens(:, :, :, :, 7); 
zR=R/dPts*[cos(tPts_ip);sin(tPts_ip)];
xr=r/dPts*[cos(tPts_ip);sin(tPts_ip)];
nor=[cos(tPts_ip);sin(tPts_ip)];
nor1Mat=repmat(nor(1,:).',1,NFourier);
nor2Mat=repmat(nor(2,:).',1,NFourier);

% RzRxr = zeros(NFourier, NFourier);
% nRzRxr1 = zeros(NFourier, NFourier);
% nRzRxr2 = zeros(NFourier, NFourier);
% nnRz1Rxr1 = zeros(NFourier, NFourier);
% nnRz1Rxr2 = zeros(NFourier, NFourier);
% nnRz2Rxr1 = zeros(NFourier, NFourier);
% nnRzR2xr2 = zeros(NFourier, NFourier);
% for i = 1: NFourier
%     zR1iO=zR(1,i)*ones(1,NFourier);
%     zR2iO=zR(2,i)*ones(1,NFourier);
%     RzRxr(:,i) = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1, Rzz        ,zR1iO,zR2iO,xr(1,:),xr(2,:),'cubic');
%     nRzRxr1(:,i) = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1, nRzz1    ,zR1iO,zR2iO,xr(1,:),xr(2,:),'cubic');
%     nRzRxr2(:,i) = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1, nRzz2    ,zR1iO,zR2iO,xr(1,:),xr(2,:),'cubic');
%     nnRz1Rxr1(:,i) = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1, nnRz1z1,zR1iO,zR2iO,xr(1,:),xr(2,:),'cubic');
%     nnRz1Rxr2(:,i) = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1, nnRz1z2,zR1iO,zR2iO,xr(1,:),xr(2,:),'cubic');
%     nnRz2Rxr1(:,i) = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1, nnRz2z1,zR1iO,zR2iO,xr(1,:),xr(2,:),'cubic');
%     nnRzR2xr2(:,i) = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1, nnRz2z2,zR1iO,zR2iO,xr(1,:),xr(2,:),'cubic');
% end

zR1Mat=repmat(zR(1,:).',1,NFourier);
zR2Mat=repmat(zR(2,:).',1,NFourier);
xr1Mat=repmat(xr(1,:),NFourier,1);
xr2Mat=repmat(xr(2,:),NFourier,1);
RzRxr = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],   Rzz      ,zR1Mat,zR2Mat,xr1Mat,xr2Mat,'makima');
nRzRxr1 = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], nRzz1    ,zR1Mat,zR2Mat,xr1Mat,xr2Mat,'makima');
nRzRxr2 = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], nRzz2    ,zR1Mat,zR2Mat,xr1Mat,xr2Mat,'makima');
nRzRxR1 = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], nRzz1    ,zR1Mat,zR2Mat,zR1Mat.',zR2Mat.','makima');
nRzRxR2 = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], nRzz2    ,zR1Mat,zR2Mat,zR1Mat.',zR2Mat.','makima');
nnRzR1xr1 = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], nnRz1z1,zR1Mat,zR2Mat,xr1Mat,xr2Mat,'makima');
nnRzR1xr2 = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], nnRz1z2,zR1Mat,zR2Mat,xr1Mat,xr2Mat,'makima');
nnRzR2xr1 = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], nnRz2z1,zR1Mat,zR2Mat,xr1Mat,xr2Mat,'makima');
nnRzR2xr2 = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], nnRz2z2,zR1Mat,zR2Mat,xr1Mat,xr2Mat,'makima');

%We could make this even faster by combining all interpolation in one
%interpolation with something like
% comb = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], 1:9, RTens, kron(...),kron(...),kron(...),kron(...),kron(...),'makima');


delxr_RzRxr = nRzRxr1.*nor1Mat + nRzRxr2.*nor2Mat;
delzRdelxr_RzRxr = nnRzR1xr1.*nor1Mat.'.*nor1Mat + nnRzR1xr2.*nor1Mat.'.*nor2Mat +...
                   nnRzR2xr1.*nor2Mat.'.*nor1Mat + nnRzR2xr2.*nor2Mat.'.*nor2Mat;
delxR_RzRxR = nRzRxR1.*nor1Mat + nRzRxR2.*nor2Mat;
RzRxR = RzRxr + (R-r)*delxr_RzRxr;

R_Mat = R_Cos*[ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]...
       +R_Sin*[zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];

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
%[DelrrGns_CosConst, DelrrGns_SinConst] = getFourierConst(DelrrGamma_NonSing);
%[DelRrGns_CosConst, DelRrGns_SinConst] = getFourierConst(DelRrGamma_NonSing);

tempNonSingGammaR=(-1i/4*besselh(0,k*R*sqrt2CosMat)-logMatp4pi);
tempNonSingGammar=(-1i/4*besselh(0,k*r*sqrt2CosMat)-logMatp4pi);
for j=1:size(tempNonSingGammaR,1)
    tempNonSingGammaR(j,j)=((2*eg-1i*pi+log(2)+2*log(k*R/2))/4/pi);
    tempNonSingGammar(j,j)=((2*eg-1i*pi+log(2)+2*log(k*r/2))/4/pi);
end

%% GETTING delyR_NrxRyR
% DDG2_FCoeffs=zeros(1,NFourier/2);
% DDG12_FCoeffs=zeros(1,NFourier/2);
% for j=1:NFourier
%     DDG2_FCoeffs(j) = 2*0.5*integral(@(t)(-k^2).*(R-r*cos(t)).*(r-R*cos(t))./(twopi*(R^2+r^2-2*R*r*cos(t)))...
%         .*log(R^2+r^2-2*R*r*cos(t)).*cos((j-1)*t),0,pi,'Waypoints',(1/2+(0:max(1,(j-1))))*pi/max(1,(j-1)));
%     DDG12_FCoeffs(j) = 2*0.5*integral(@(t)(-k^2).*(-2*R*r+(R^2+r^2).*cos(t))./(twopi*(R^2+r^2-2*R*r*cos(t)))...
%         .*log(k*sqrt(R^2+r^2-2*R*r*cos(t))).*cos((j-1)*t),0,pi,'Waypoints',(1/2+(0:max(1,(j-1))))*pi/max(1,(j-1)));
% end


%SingSingConv  = 2*1/4/pi*r/R*(r^2-R^2*CosMat)./(R^4+r^4-2*R^2*r^2*CosMat);
%      SingSingConv2 = (cos((tPts_ip).'*(0:NFourier-1))) * (repmat(DDG2_FCoeffs.'.*logRrp2pi_FCoeffs.',1,NFourier).*(cos((0:NFourier-1).'*(tPts_ip))))...
%         +(sin((tPts_ip).'*(0:NFourier-1))) * (repmat(DDG2_FCoeffs.'.*logRrp2pi_FCoeffs.',1,NFourier).*(sin((0:NFourier-1).'*(tPts_ip))));
SingSingConv2_CosConst=cos(tPts_ip.'*(0:NFourier-1)).*repmat(DDG2_FCoeffs.*logRrp2pi_FCoeffs,NFourier,1);
SingSingConv2_SinConst=sin(tPts_ip.'*(0:NFourier-1)).*repmat(DDG2_FCoeffs.*logRrp2pi_FCoeffs,NFourier,1);

%     SingSingConv3 = (cos((theta.'-pi)*(0:NFourier/2-1))) * (repmat(DDG12_FCoeffs.'.*logRrp2pi_FCoeffs.',1,NumDiscPts).*(cos((0:NFourier/2-1).'*(theta-pi))))...
%         +(sin((theta.'-pi)*(0:NFourier/2-1))) * (repmat(DDG12_FCoeffs.'.*logRrp2pi_FCoeffs.',1,NumDiscPts).*(sin((0:NFourier/2-1).'*(theta-pi))));
SingSingConv3_CosConst=cos(tPts_ip.'*(0:NFourier-1)).*repmat(DDG12_FCoeffs.*logRrp2pi_FCoeffs,NFourier,1);
SingSingConv3_SinConst=sin(tPts_ip.'*(0:NFourier-1)).*repmat(DDG12_FCoeffs.*logRrp2pi_FCoeffs,NFourier,1);


smoothDelDelGammaTerm = -(1/4)*1i.*k.*(...
    k.*(-R + r.*CosMat).*(r- R.*CosMat)./(sqrtRrCosMat.^2).*(besselh(0,k.*sqrtRrCosMat)-2i./pi.*log(sqrtRrCosMat)) ...
    - ((-2.*r.*R + (r.^2 + R.^2).*CosMat)./(sqrtRrCosMat.^3).*(besselh(1,k.*sqrtRrCosMat)+2i./pi./(k.*sqrtRrCosMat)-(2i*k.*sqrtRrCosMat.*log(k.*sqrtRrCosMat))./twopi)));

NSmoo_CosConst = 2*R_Cos;
NSmoo_SinConst = 2*R_Sin;
Inc3_CosConst=NSmoo_CosConst.*repmat([0,(r/R).^(0:NFourier-2).*(1:NFourier-1)/2/R^2],NFourier,1); %critical!
Inc3_SinConst=NSmoo_SinConst.*repmat([0,(r/R).^(0:NFourier-2).*(1:NFourier-1)/2/R^2],NFourier,1); %critical!

Inc4_CosConst=NSmoo_CosConst.*repmat(DDG2_FCoeffs,NFourier,1);
Inc4_SinConst=NSmoo_SinConst.*repmat(DDG2_FCoeffs,NFourier,1);

Inc5_CosConst=NSmoo_CosConst.*repmat(DDG12_FCoeffs,NFourier,1);
Inc5_SinConst=NSmoo_SinConst.*repmat(DDG12_FCoeffs,NFourier,1);

Inc6  =(smoothDelDelGammaTerm.'*(((2*(R_Mat))+2*(log(R^2+r^2-2*R*r*CosMat)/4/pi-log(2*R*r)/4/pi)).'.*(sigmaMat_r))).'; %This might blow up?!
Inc6R =(delzRdelxr_RzRxr.'*(((2*(R_Mat))+2*(log(R^2+r^2-2*R*r*CosMat)/4/pi-log(2*R*r)/4/pi)).'.*(sigmaMat_r))).'; %This might blow up?!
[Inc6_CosConst, Inc6_SinConst] = getFourierConst(-Inc6-Inc6R);

[delxR_RzRxR_CosConst,delxR_RzRxR_SinConst] = getFourierConst(delxR_RzRxR);
DelInc_CosConst=-r*(SingSingConv2_CosConst+SingSingConv3_CosConst+Inc3_CosConst+Inc4_CosConst+Inc5_CosConst)...
    +[Inc6_CosConst + DelRRGns_CosConst+delxR_RzRxR_CosConst,zeros(NFourier,NFourier/2)]; % ZEROS extension
DelInc_SinConst=-r*(SingSingConv2_SinConst+SingSingConv3_SinConst+Inc3_SinConst+Inc4_SinConst+Inc5_SinConst)...
    +[Inc6_SinConst + DelRRGns_SinConst+delxR_RzRxR_SinConst,zeros(NFourier,NFourier/2)]; % ZEROS extension

%% GETTING NRxRyR
%temp =-2*(DelIncNonSing.'  *  ((R_Mat+logMatp4pi).'.*(sigmaMat_init*R/r_init))).';
temp=(-2*R)*(R_Cos*(DelInc_CosConst.*IntegConst).'+R_Sin*(DelInc_SinConst.*IntegConst).'+...
    +log_CosConst*(DelInc_CosConst.*IntegConst).'+log_SinConst*(DelInc_SinConst.*IntegConst).');
[temp_CosConst,temp_SinConst]=getFourierConst(temp);

temp2=-2*((1i*k/4*(r-R*CosMat)./sqrtRrCosMat.*besselh(1,k*sqrtRrCosMat)-(r-R*CosMat)./(R^2+r^2-2*R*r*CosMat)./twopi).'  *  ...
    ((R_Mat+log((R^2+r^2)/(2*R*r)-CosMat)./fourpi).'.*(sigmaMat_r))).';
temp2R=-2*(delxr_RzRxr.' * ((R_Mat+log((R^2+r^2)/(2*R*r)-CosMat)./fourpi).'.*(sigmaMat_r))).';
[temp2_CosConst, temp2_SinConst] = getFourierConst(temp2+temp2R);
[GammaNS_CosConst,  GammaNS_SinConst]  = getFourierConst(tempNonSingGammaR + RzRxR(1:NFourier,:).');

R_Cos = 0.5*([GammaNS_CosConst,zeros(NFourier,NFourier/2)]...
    +R_Cos.*repmat([0,1+(r/R).^(1:NFourier-1)-(r/R).^(2*(1:NFourier-1))],NFourier,1)...
    +(log(2)/fourpi).*[ones(NFourier,1),zeros(NFourier,NFourier-1)]...
    +[temp_CosConst+temp2_CosConst,zeros(NFourier,NFourier/2)]); % ZEROS extension

R_Sin = 0.5*([GammaNS_SinConst,zeros(NFourier,NFourier/2)]...
    +R_Sin.*repmat([0,1+(r/R).^(1:NFourier-1)-(r/R).^(2*(1:NFourier-1))],NFourier,1)...
    +[temp_SinConst + temp2_SinConst,zeros(NFourier,NFourier/2)]);  % ZEROS extension

% RMat_new =  R_CosConst_new * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]...
%     + R_SinConst_new * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];

%% Post Loop Computations
%assert(norm(RMat_new)<1e10)

end



%%How to nD interpolate:
% V = zeros(3,3,3,3);
% V(:,:,1,1)=(-1:1:1).'+(-1:1:1) -1 -1;
% V(:,:,1,2)=(-1:1:1).'+(-1:1:1) -1 -0;
% V(:,:,1,3)=(-1:1:1).'+(-1:1:1) -1 +1;
% V(:,:,2,1)=(-1:1:1).'+(-1:1:1) -0 -1;
% V(:,:,2,2)=(-1:1:1).'+(-1:1:1) -0 +0;
% V(:,:,2,3)=(-1:1:1).'+(-1:1:1) -0 +1;
% V(:,:,3,1)=(-1:1:1).'+(-1:1:1) +1 -1;
% V(:,:,3,2)=(-1:1:1).'+(-1:1:1) +1 -0;
% V(:,:,3,3)=(-1:1:1).'+(-1:1:1) +1 +1;
% vq = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1,V,...
%     [0.5,0],[0.5,-0.75],[0.5,0.75],[0.5,0.42],'linear');
%Yields: vq == [2,0.42]

%Extrapolate: Do not extrapolate
% vq = interpn(-1:1:1,-1:1:1,-1:1:1,-1:1:1,V,...
%     [0.5,0],[0.5,-0.75],[0.5,0.75],[0.5,0.42],'linear','extrap');

