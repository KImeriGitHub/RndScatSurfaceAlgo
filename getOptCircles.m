function [NMat, CenterCoord, DiskRad, min_val_Data, BigN] = getOptCircles(U,Y, k, z, x, r_inc_Data, itermax, NZeta, NFourier)
%PRE:  U is a n x p complex matrix, n=size(z,2)
%      Y is a n x 10  {0,1} matrix, 10=size(x,2) 
%      k is a real pos number, not zero
%      z is a 2 x n array of real numbers
%      x is a 2 x 10 array of real numbers
%      r_inc_Data is row of pos inc real  numbers
%      itermax, NZeta pos Integer, >2 
%      NFourier == 2^(integer) pos Integer, 
%POST: 
%      NMat is a NFourier x NFourier complex matrix, N=size(SourcePts,2)
%      CenterCoord is a 2 x itermax array of real numbers
%      DiskRad is a 1 x itermax real array
%      min_val_Step is a length(r_inc_Data) x (itermax +1) pos real array
%      BigN is a itermax*NFourier x itermax*NFourier complex matrix
%Assert: z in (-1.15, -1.05) x (-1, 1)
%        x in (-1.15, -1.05) x (-1, 1)

%Desc: Trying to set disks in (-1,1) x (-1,1) at the locations given by 
%CenterCoord that those do not overlap and such that their 
%NMat matrix is close enough to || U*N - Y||.
%BigN is the Neumann function N(DiskBoundary,DiskBoundary).

%% Notes

%% Changeable configuration
r_max = r_inc_Data(end);
r_end = r_max;
r_init = r_inc_Data(1);
%Disk2DiskThresh = 2*(sqrt(1/4/itermax)-r_end); % Numerically stable if dis>r/4 (experience) 
Disk2DiskThresh = 0.01;
Pt2PtThresh = 0.01; % ideally dist>lambda/4 (but often not possible. We need N < 2/Pt2PtThresh^2)
%gradNum = 10;
%gradDist = 0.0001;
assert(itermax>2);

min_val_Data=zeros(1,itermax+1);
CenterCoord = zeros(2,itermax);
DiskRad = zeros(1,itermax);
ALL_DiscPts_x = zeros(1,itermax*NFourier);
ALL_DiscPts_y = zeros(1,itermax*NFourier);

%ZetaGrid is a row of NZeta points in (-1,1) x (-1,1) chosen randomly
% ZetaRow = [kron(ones(1,NZeta),linspace(-1,1,NZeta));
%            kron(linspace(-1,1,NZeta),ones(1,NZeta))];
ZetaRow = zeros(2,NZeta);
ZR_loop = 0;
whileCounter = 0;
while ZR_loop < NZeta
    RNGpt = (rand(2,1)-0.5)*2;
    if ~any(sum((ZetaRow(:,1:ZR_loop)-RNGpt).^2)<(Pt2PtThresh)^2) ...
      && ~any(sum((CenterCoord(:,1:0)-RNGpt).^2)<(r_max*(sqrt(2)+1)+Disk2DiskThresh)^2)
        ZetaRow(:,ZR_loop+1) = RNGpt;
        ZR_loop = ZR_loop + 1;
    end
    whileCounter = whileCounter+1;
    if whileCounter > 100*NZeta
        error('No initally location found. Check parameter.');
    end
end

% Extra variables
tPts_ip = linspace(-pi,pi,NFourier+1);
tPts_ip = tPts_ip(1:end-1);
nor=[cos(tPts_ip); sin(tPts_ip)];

eg = 0.57721566490153286060651209008240243104215933593992;
twopi=2*pi;
fourpi=4*pi;

CosM = cos((tPts_ip)-(tPts_ip).');
LMat = log(1-CosM)/fourpi;
for i = 1:NFourier
    LMat(i,i) = (log(tPts_ip(2)-tPts_ip(1))-2)/twopi;
end
% LMat(1,:)*ones(NFourier,1)*(tPts_ip(2)-tPts_ip(1)) should give -0.346574

%PRE ALLOC
IpVals_R=zeros(5,5,5,5,7);
BigN = zeros([NFourier*itermax, NFourier*itermax]); %TOO BIG

sigmaVal_nor = twopi/NFourier;

GzCirc = zeros(size(z,2),itermax*NFourier);%TOO BIG
Gz0Circ1 = zeros(size(z,2),itermax*NFourier);%TOO BIG
Gz0Circ2 = zeros(size(z,2),itermax*NFourier);%TOO BIG
GxCirc = zeros(size(x,2),itermax*NFourier);%TOO BIG
Gx0Circ1 = zeros(size(x,2),itermax*NFourier);%TOO BIG
Gx0Circ2 = zeros(size(x,2),itermax*NFourier);%TOO BIG
Gz0DelCirc = zeros(size(z,2),itermax*NFourier);%TOO BIG
Gx0DelCirc = zeros(size(x,2),itermax*NFourier);%TOO BIG

%% First Step
%Source - Receiver Connection
[Gzx, ~, ~, ~, ~, ~, ~, ~, ~, ~] = getGamma(z, x,k);
Gzx = Gzx.';
%Source - Zeta Connection
[GzZeta, Gz0Zeta1, Gz0Zeta2, ~, ~, ...
    ~, ~, ~, ~, ~] = getGamma(z, ZetaRow,k);
GzZeta=GzZeta.';Gz0Zeta1=Gz0Zeta1.';Gz0Zeta2=Gz0Zeta2.';
%Receiver - Zeta Connection
[GxZeta, Gx0Zeta1, Gx0Zeta2, ~, ~, ...
    ~, ~, ~, ~, ~] = getGamma(x, ZetaRow,k);
GxZeta=GxZeta.';Gx0Zeta1=Gx0Zeta1.';Gx0Zeta2=Gx0Zeta2.';

Zeta_Idx = 0;
CostMat = abs(U*Gzx).^2;
sigCost = 1./(1+exp(-CostMat));
min_val = sum((sigCost-Y).^2,'all');
min_val_Data(1,1)=min_val;
for i=1:NZeta
    NMat_add = r_init^2*pi*(k^2*GzZeta(:,i)*GxZeta(:,i).' ...
        - 2*Gz0Zeta1(:,i)*Gx0Zeta1(:,i).' - 2*Gz0Zeta2(:,i)*Gx0Zeta2(:,i).');
    CostMattemp = abs(U*(Gzx+NMat_add)).^2;
    sigCostTemp = 1./(1+exp(-CostMattemp));
    val_opt = sum((sigCostTemp-Y).^2,'all');
    
    if val_opt < min_val
        Zeta_Idx = i;
        min_val = val_opt;
    end
end

if Zeta_Idx == 0
    Zeta_Idx = randi(NZeta);
%     DiskRad   = nan;
%     CenterCoord = [nan;nan];
%     NMat=nan(size(z,2),size(x,2));
%     warning('FirstStep failed')
%     return;
end

CenterCoord(:,1) = ZetaRow(:,Zeta_Idx);
% ZetaX_idx = floor((Zeta_Idx-1)/size(ZetaX,1))+1;
% ZetaY_idx = mod((Zeta_Idx-1),size(ZetaX,1))+1;
% CenterIdx = [ZetaX_idx; ZetaY_idx];
% Zeta_Idx == (ZetaX_idx-1)*size(X,1)+ZetaY_idx;

NMat = Gzx + r_init^2*pi*(k^2*GzZeta(:,Zeta_Idx)*GxZeta(:,Zeta_Idx).' ...
        - 2*Gz0Zeta1(:,Zeta_Idx)*Gx0Zeta1(:,Zeta_Idx).' ...
        - 2*Gz0Zeta2(:,Zeta_Idx)*Gx0Zeta2(:,Zeta_Idx).');

%% PRE Computation
DDG2_FCoeffs_PREC=zeros(length(r_inc_Data),NFourier/2);
DDG12_FCoeffs_PREC=zeros(length(r_inc_Data),NFourier/2);
for j=1:NFourier
    for ii = 2:length(r_inc_Data)
        r=r_inc_Data(ii-1);
        R=r_inc_Data(ii);
        DDG2_FCoeffs_PREC(ii,j) = 2*0.5*integral(@(t)(-k^2).*(R-r*cos(t)).*(r-R*cos(t))./(twopi*(R^2+r^2-2*R*r*cos(t)))...
            .*log(R^2+r^2-2*R*r*cos(t)).*cos((j-1)*t),0,pi,'Waypoints',(1/2+(0:max(1,(j-1))))*pi/max(1,(j-1)));
        DDG12_FCoeffs_PREC(ii,j) = 2*0.5*integral(@(t)(-k^2).*(-2*R*r+(R^2+r^2).*cos(t))./(twopi*(R^2+r^2-2*R*r*cos(t)))...
            .*log(k*sqrt(R^2+r^2-2*R*r*cos(t))).*cos((j-1)*t),0,pi,'Waypoints',(1/2+(0:max(1,(j-1))))*pi/max(1,(j-1)));
    end
end

%% PreAlloc
% Circle Bdry Data
circleDATA_Sin=zeros([NFourier, NFourier, itermax]); %TOO BIG
circleDATA_Cos=zeros([NFourier, NFourier, itermax]); %TOO BIG

%PreAlloc Gamma(Zeta,Disc)
% GZetaDisc = zeros(NZeta,    itermax*NFourier);
% GZeta0DelDisc = zeros(NZeta,itermax*NFourier);
% GZeta1Disc0 = zeros(NZeta,  itermax*NFourier);
% GZeta2Disc0 = zeros(NZeta,  itermax*NFourier);
% GZeta1DelDisc = zeros(NZeta,itermax*NFourier);
% GZeta2DelDisc = zeros(NZeta,itermax*NFourier);

CC1 = CenterCoord(1,1);
CC2 = CenterCoord(2,1);
% surround_Coords = [CC1-r_end, CC1, CC1+r_end, CC1-r_end, CC1, CC1+r_end, CC1-r_end, CC1, CC1+r_end;
%                    CC2-r_end, CC2-r_end, CC2-r_end, CC2, CC2, CC2, CC2+r_end, CC2+r_end, CC2+r_end];
DiscPts = [r_init*cos(tPts_ip).'+ CC1,...
           r_init*sin(tPts_ip).'+ CC2].';

N0CircSurR  = zeros(NFourier, NFourier);
N0CircSurR1 = zeros(NFourier, NFourier);
N0CircSurR2 = zeros(NFourier, NFourier);

N0CircDelSurR = N0CircSurR1.*repmat(nor(1,:),NFourier,1)+N0CircSurR2.*repmat(nor(2,:),NFourier,1);

%% Main Loop
RCC = 0;
NZetaDim = NZeta;   %Used if we have not enough locations found.
for i=1:itermax  
    disp(['i = ',num2str(i),' of ',num2str(itermax)]);
    %% Initilization for new loop
    CirciIdx = (1:NFourier)+(i-1)*NFourier;
    circleDATA_Sin(:,:,i) = zeros([NFourier, NFourier]);
    circleDATA_Cos(:,:,i) = [(log(k*r_init)/4/pi + eg/4/pi - 1i/8 + RCC/2)...
            *ones(NFourier,1), zeros(NFourier,NFourier-1)];
    for ii = 2:length(r_inc_Data)
        %Info: 'n' means 'new' and refers to the update
        % 'R' is only for circle values and refers to Pts at circle with
        % radius R but domain with circle r
        
        % Actual increase
        r=r_inc_Data(ii-1);
        R=r_inc_Data(ii);
        [circleDATA_Cos(:,:,i), circleDATA_Sin(:,:,i)] = getDiscInfl(circleDATA_Cos(:,:,i), circleDATA_Sin(:,:,i),...
            IpVals_R, r, R, k, r_end, DDG2_FCoeffs_PREC(ii,:), DDG12_FCoeffs_PREC(ii,:));
        RMatn =  circleDATA_Cos(:,:,i) * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]... %'n' stand for new
               + circleDATA_Sin(:,:,i) * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
    end
    
    %% UPDATE N(y,w)
    BigN(CirciIdx,CirciIdx) = RMatn + LMat;
    %%% Update old circles
    for j=1:(i-1)
        temp = -0.5*(N0CircSurR((1:NFourier) + NFourier*(j-1),:))*(N0CircDelSurR((1:NFourier) + NFourier*(j-1),:).' *(sigmaVal_nor*r_end))...
            +0.5*(N0CircDelSurR((1:NFourier) + NFourier*(j-1),:)*(2*BigN((1:NFourier)+(i-1)*NFourier,(1:NFourier)+(i-1)*NFourier))*N0CircDelSurR((1:NFourier) + NFourier*(j-1),:).') *(sigmaVal_nor*r_end)^2;
        circleDATA_Cos(:,:,j) = circleDATA_Cos(:,:,j) + [getFourierConst(temp),zeros(NFourier,NFourier/2)];
        RMatn =  (circleDATA_Cos(:,:,j)) * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]... %'n' stand for new
                + circleDATA_Sin(:,:,j) * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
        BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(j-1)*NFourier) = RMatn+LMat;
    end
    
    %%% New CrossActions
    for j=1:(i-1)
        BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(i-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(i-1)*NFourier)...
            - (N0CircDelSurR((1:NFourier) + NFourier*(j-1),:))*(BigN((1:NFourier)+(i-1)*NFourier,(1:NFourier)+(i-1)*NFourier).' *(sigmaVal_nor*r_end));
        BigN((1:NFourier)+(i-1)*NFourier,(1:NFourier)+(j-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(i-1)*NFourier).';
    end
    
    %%% Update Old CrossActions
    for j=1:(i-1)
        for jj=1:(j-1)
            BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier)...
                - (N0CircDelSurR((1:NFourier) + NFourier*(j-1),:))*(BigN((1:NFourier)+(jj-1)*NFourier,(1:NFourier)+(i-1)*NFourier).'*(sigmaVal_nor*r_end));
            BigN((1:NFourier)+(jj-1)*NFourier,(1:NFourier)+(j-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier).';
        end
    end
    
    
    %% Randomly generating new possible locations
    ZR_loop = 0;
    whileCounter = 0;
    ZetaRow = zeros(2,NZeta);
    while ZR_loop<NZeta
        RNGpt = (rand(2,1)-0.5)*2;
        if ~any(sum((ZetaRow(:,1:ZR_loop)-RNGpt).^2)<(Pt2PtThresh)^2) ...
                && ~any(sum((CenterCoord(:,1:i)-RNGpt).^2)<(r_max*(sqrt(2)+1)+Disk2DiskThresh)^2)
            ZetaRow(:,ZR_loop+1) = RNGpt;
            ZR_loop = ZR_loop + 1;
        end
        whileCounter = whileCounter+1;
        if whileCounter>100*NZeta
            if ZR_loop>0
                warning('Not full NZeta locations found. Continuing with less.');
                NZetaDim = ZR_loop;
                ZetaRow = zeros(2,NZetaDim);
                disp(['NZetaDom = ',num2str(NZetaDim)]);
            else
                warning('No further rnd location found. Terminating.');
                BigN         = BigN(1:NFourier*i,1:NFourier*i);
                CenterCoord  = CenterCoord(:,1:i);
                DiskRad      = DiskRad(1,1:i);
                min_val_Data = min_val_Data(:,1:i+1);
                return;
            end
        end
    end
    
    %% Preparing N(z, zeta), N(x, zeta)
    DiscPts = [r_end*cos(tPts_ip).'+ CenterCoord(1,i),...
               r_end*sin(tPts_ip).'+ CenterCoord(2,i)].';
    ALL_DiscPts_x(CirciIdx) = DiscPts(1,:);
    ALL_DiscPts_y(CirciIdx) = DiscPts(2,:);
    
    % Source- Circle Conecction
    [GzCirc(:,CirciIdx), Gz0Circ1(:,CirciIdx), Gz0Circ2(:,CirciIdx),...
        ~, ~, ~, ~, ~, ~, ~] = getGamma(DiscPts, z,k);
    Gz0DelCirc(:,CirciIdx) = (Gz0Circ1(:,CirciIdx)*sigmaVal_nor*r_end).*repmat(nor(1,:),size(z,2),1)...
        +(Gz0Circ2(:,CirciIdx)*sigmaVal_nor*r_end).*repmat(nor(2,:),size(z,2),1);
    % Receiver- Circle Conecction
    [GxCirc(:,CirciIdx), Gx0Circ1(:,CirciIdx), Gx0Circ2(:,CirciIdx),...
        ~, ~, ~, ~, ~, ~, ~] = getGamma(DiscPts, x,k);
    Gx0DelCirc(:,CirciIdx) = (Gx0Circ1(:,CirciIdx)*sigmaVal_nor*r_end).*repmat(nor(1,:),size(x,2),1)...
        +(Gx0Circ2(:,CirciIdx)*sigmaVal_nor*r_end).*repmat(nor(2,:),size(x,2),1);
    %Source - Zeta Connection
    [GzZeta, Gz0Zeta1, Gz0Zeta2, ~, ~, ...
        ~, ~, ~, ~, ~] = getGamma(z, ZetaRow,k);
    GzZeta=GzZeta.';Gz0Zeta1=Gz0Zeta1.';Gz0Zeta2=Gz0Zeta2.';
    %Receiver - Zeta Connection
    [GxZeta, Gx0Zeta1, Gx0Zeta2, ~, ~, ...
        ~, ~, ~, ~, ~] = getGamma(x, ZetaRow,k);
    GxZeta=GxZeta.';Gx0Zeta1=Gx0Zeta1.';Gx0Zeta2=Gx0Zeta2.';
    
    [~, GZeta0Disc1Temp, GZeta0Disc2Temp, ~, ~, ...
        GZeta1Disc1Temp, GZeta2Disc1Temp, GZeta1Disc2Temp, GZeta2Disc2Temp, ~] ...
        = getGamma([ALL_DiscPts_x(1:NFourier*i);ALL_DiscPts_y(1:NFourier*i)], ZetaRow,k);
    
    %GZeta1Disc0(:,CirciIdx) = GZeta1Disc0Temp*(sigmaVal_nor*R);
    %GZeta2Disc0(:,CirciIdx) = GZeta2Disc0Temp*(sigmaVal_nor*R);
    
    DiskRad(i) = r_end;
    integConst = kron(sigmaVal_nor*DiskRad(1:i),ones(1,NFourier));
    GZeta0DelDisc =  (GZeta0Disc1Temp.*repmat(integConst,NZetaDim,1)).*repmat(nor(1,:),NZetaDim,i)...
                   + (GZeta0Disc2Temp.*repmat(integConst,NZetaDim,1)).*repmat(nor(2,:),NZetaDim,i);
    GZeta1DelDisc =  (GZeta1Disc1Temp.*repmat(integConst,NZetaDim,1)).*repmat(nor(1,:),NZetaDim,i)...
                   + (GZeta1Disc2Temp.*repmat(integConst,NZetaDim,1)).*repmat(nor(2,:),NZetaDim,i);
    GZeta2DelDisc =  (GZeta2Disc1Temp.*repmat(integConst,NZetaDim,1)).*repmat(nor(1,:),NZetaDim,i)...
                   + (GZeta2Disc2Temp.*repmat(integConst,NZetaDim,1)).*repmat(nor(2,:),NZetaDim,i);
    
    NMat = Gzx - GzCirc(:,1:i*NFourier)*Gx0DelCirc(:,1:i*NFourier).' + Gz0DelCirc(:,1:i*NFourier)*BigN(1:NFourier*i,1:NFourier*i)*Gx0DelCirc(:,1:i*NFourier).';
    NzZeta = GzZeta - GzCirc(:,1:i*NFourier)*GZeta0DelDisc(:,1:i*NFourier).' + Gz0DelCirc(:,1:i*NFourier)*BigN(1:NFourier*i,1:NFourier*i)*GZeta0DelDisc(:,1:i*NFourier).';
    NxZeta = GxZeta - GxCirc(:,1:i*NFourier)*GZeta0DelDisc(:,1:i*NFourier).' + Gx0DelCirc(:,1:i*NFourier)*BigN(1:NFourier*i,1:NFourier*i)*GZeta0DelDisc(:,1:i*NFourier).';
    NzZeta1 = Gz0Zeta1 - GzCirc(:,1:i*NFourier)*GZeta1DelDisc(:,1:i*NFourier).' + Gz0DelCirc(:,1:i*NFourier)*BigN(1:NFourier*i,1:NFourier*i)*GZeta1DelDisc(:,1:i*NFourier).';
    NxZeta1 = Gx0Zeta1 - GxCirc(:,1:i*NFourier)*GZeta1DelDisc(:,1:i*NFourier).' + Gx0DelCirc(:,1:i*NFourier)*BigN(1:NFourier*i,1:NFourier*i)*GZeta1DelDisc(:,1:i*NFourier).';
    NzZeta2 = Gz0Zeta2 - GzCirc(:,1:i*NFourier)*GZeta2DelDisc(:,1:i*NFourier).' + Gz0DelCirc(:,1:i*NFourier)*BigN(1:NFourier*i,1:NFourier*i)*GZeta2DelDisc(:,1:i*NFourier).';
    NxZeta2 = Gx0Zeta2 - GxCirc(:,1:i*NFourier)*GZeta2DelDisc(:,1:i*NFourier).' + Gx0DelCirc(:,1:i*NFourier)*BigN(1:NFourier*i,1:NFourier*i)*GZeta2DelDisc(:,1:i*NFourier).';

    
    %% Searching for a new location
    Zeta_Idx = 0;
    CostMat = abs(U*NMat).^2;
    sigCost = 1./(1+exp(-CostMat));
    min_val = sum((sigCost-Y).^2,'all');
    if min_val>min_val_Data(1,i) && length(r_inc_Data)>5
        r_end = r_inc_Data(end-1);
        r_inc_Data = r_inc_Data(1:end-1);
    end
    if r_end == r_inc_Data(5)
        BigN         = BigN(1:NFourier*i,1:NFourier*i);       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CenterCoord  = CenterCoord(:,1:i);                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DiskRad      = DiskRad(1,1:i);                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        min_val_Data = min_val_Data(:,1:i+1);                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        return;                                               %%%%%% RETURN CONDITION
    end                                                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    min_val_Data(1,i+1) = min_val;
    disp(['min_val is ',num2str(min_val)])
    
    
    if i==itermax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        break;    %%%%%% BREAK CONDITION
    end           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    UN =U*NMat;
    temp1 = U * NzZeta;
    temp2 = U * NzZeta1;
    temp3 = U * NzZeta2;
    for ii=1:NZetaDim
        UNMat_add = r_init^2*pi*(k^2*temp1(:,ii)*NxZeta(:,ii).' ...
            - 2*temp2(:,ii)*NxZeta1(:,ii).' - 2*temp3(:,ii)*NxZeta2(:,ii).');
        CostMattemp = abs(UN+UNMat_add).^2;
        sigCostTemp = 1./(1+exp(-CostMattemp));
        val_opt = sum((sigCostTemp-Y).^2,'all');
        
        if val_opt < min_val
            Zeta_Idx = ii;
            min_val = val_opt;
        end
    end
    if Zeta_Idx == 0
        disp(['Zeta_Idx is ',num2str(0)])
        BigN         = BigN(1:NFourier*i,1:NFourier*i);
        CenterCoord  = CenterCoord(:,1:i);
        DiskRad      = DiskRad(1,1:i);
        min_val_Data = min_val_Data(:,1:i+1);
        return;
    end
    disp(['min_val after location-finding is ',num2str(min_val)])
    disp(['r_end is ',num2str(r_end)])
    CenterCoord(:,i+1) = ZetaRow(:,Zeta_Idx);
 
    ZI=Zeta_Idx;
    
    NMat = NMat + r_init^2*pi*(k^2*NzZeta(:,ZI)*NxZeta(:,ZI).' ...
            - 2*NzZeta1(:,ZI)*NxZeta1(:,ZI).' ...
            - 2*NzZeta2(:,ZI)*NxZeta2(:,ZI).');

    %% R Vals  
    CC = ZetaRow(:,ZI);
    CC1 = CC(1); 
    CC2 = CC(2);
    surround_Coords = [CC1-r_end, CC1-r_end/2, CC1, CC1+r_end/2, CC1+r_end,...
        CC1-r_end, CC1-r_end/2, CC1, CC1+r_end/2, CC1+r_end,...
        CC1-r_end, CC1-r_end/2, CC1, CC1+r_end/2, CC1+r_end,...
        CC1-r_end, CC1-r_end/2, CC1, CC1+r_end/2, CC1+r_end,...
        CC1-r_end, CC1-r_end/2, CC1, CC1+r_end/2, CC1+r_end;...
                       CC2-r_end,CC2-r_end,CC2-r_end,CC2-r_end,CC2-r_end, ...
        CC2-r_end/2,CC2-r_end/2,CC2-r_end/2,CC2-r_end/2,CC2-r_end/2,...
        CC2, CC2, CC2, CC2, CC2,...
        CC2+r_end/2,CC2+r_end/2,CC2+r_end/2,CC2+r_end/2,CC2+r_end/2,...
        CC2+r_end,CC2+r_end,CC2+r_end,CC2+r_end,CC2+r_end];
    
    [GammaSurCirct, GammaSur1Circt, GammaSur2Circt, GammaSurCirc1t, GammaSurCirc2t, ...
        GammaSur1Circ1t, GammaSur1Circ2t, GammaSur2Circ1t, GammaSur2Circ2t, ~] ...
        = getGamma(surround_Coords, [ALL_DiscPts_x(1:NFourier*i);ALL_DiscPts_y(1:NFourier*i)], k);
    
    [RVals, ~, ~, nabla0nabla1RVals, nabla0nabla2RVals, ...
        nabla1nabla1RVals, nabla1nabla2RVals, nabla2nabla1RVals, nabla2nabla2RVals] ...
        = getEvalCircReum(surround_Coords, surround_Coords, BigN(1:NFourier*i,1:NFourier*i), ...
            CenterCoord(:,1:i), DiskRad(1:i), NFourier,k);
    
    GammaSurDelCirc  = GammaSurCirc1t.'.*repmat(nor(1,:),25,i)+GammaSurCirc2t.'.*repmat(nor(2,:),25,i);
    GammaSur1DelCirc = GammaSur1Circ1t.'.*repmat(nor(1,:),25,i)+GammaSur1Circ2t.'.*repmat(nor(2,:),25,i);
    GammaSur2DelCirc = GammaSur2Circ1t.'.*repmat(nor(1,:),25,i)+GammaSur2Circ2t.'.*repmat(nor(2,:),25,i);
    
    
    % This is actually very difficult logistical problem
    IpVals_R(:,:,:,:,1) = reshape(RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,2) = reshape(nabla0nabla1RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,3) = reshape(nabla0nabla2RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,4) = reshape(nabla1nabla1RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,5) = reshape(nabla1nabla2RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,6) = reshape(nabla2nabla1RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,7) = reshape(nabla2nabla2RVals.',5,5,5,5);
    RCC = RVals(13,13);
    
    %% Interpolated N0(Old Circle, New Circle)
    N0CircSur  = (GammaSurCirct.' - GammaSurDelCirc.*repmat(integConst,25,1)*(BigN(1:NFourier*i,1:NFourier*i)*2)).';
    N0CircSur1 = (GammaSur1Circt.' - GammaSur1DelCirc.*repmat(integConst,25,1)*(BigN(1:NFourier*i,1:NFourier*i)*2)).';
    N0CircSur2 = (GammaSur2Circt.' - GammaSur2DelCirc.*repmat(integConst,25,1)*(BigN(1:NFourier*i,1:NFourier*i)*2)).';  
    
    N0CircSurR  = zeros(NFourier*i,NFourier);
    N0CircSurR1 = zeros(NFourier*i,NFourier);
    N0CircSurR2 = zeros(NFourier*i,NFourier);
%     for n=1:NFourier
%         N0CircSurR(n,:)  = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircSur(n,:) ,5,5),cos(tPts_ip),sin(tPts_ip),'makima');
%         N0CircSurR1(n,:) = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircSur1(n,:),5,5),cos(tPts_ip),sin(tPts_ip),'makima');
%         N0CircSurR2(n,:) = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircSur2(n,:),5,5),cos(tPts_ip),sin(tPts_ip),'makima');
%     end
    % Reshaping
    reshN0CircSur=zeros(5,5,NFourier*i);
    reshN0CircSur1=zeros(5,5,NFourier*i);
    reshN0CircSur2=zeros(5,5,NFourier*i);
    for n=1:NFourier*i
        reshN0CircSur(:,:,n)=reshape(N0CircSur(n,:) ,5,5);
        reshN0CircSur1(:,:,n)=reshape(N0CircSur1(n,:) ,5,5);
        reshN0CircSur2(:,:,n)=reshape(N0CircSur2(n,:) ,5,5);
    end
    tempIV = [-1,-0.5,0,0.5,1];
    N0CircSurRip = interpn(tempIV,tempIV, linspace(0,1,NFourier*i), reshN0CircSur, repmat(cos(tPts_ip),1,NFourier*i),repmat(sin(tPts_ip),1,NFourier*i), kron(linspace(0,1,NFourier*i),ones(1,NFourier)),'makima');
    N0CircSurR1ip = interpn(tempIV,tempIV, linspace(0,1,NFourier*i), reshN0CircSur1, repmat(cos(tPts_ip),1,NFourier*i),repmat(sin(tPts_ip),1,NFourier*i), kron(linspace(0,1,NFourier*i),ones(1,NFourier)),'makima');
    N0CircSurR2ip = interpn(tempIV,tempIV, linspace(0,1,NFourier*i), reshN0CircSur2, repmat(cos(tPts_ip),1,NFourier*i),repmat(sin(tPts_ip),1,NFourier*i), kron(linspace(0,1,NFourier*i),ones(1,NFourier)),'makima');
    
    % BackShaping
    for n=1:NFourier*i
        N0CircSurR(n,:) = N0CircSurRip((1:NFourier) + (n-1)*NFourier);
        N0CircSurR1(n,:) = N0CircSurR1ip((1:NFourier) + (n-1)*NFourier);
        N0CircSurR2(n,:) = N0CircSurR2ip((1:NFourier) + (n-1)*NFourier);
    end
    
    N0CircDelSurR = N0CircSurR1.*repmat(nor(1,:),NFourier*i,1)+N0CircSurR2.*repmat(nor(2,:),NFourier*i,1);
end
BigN         = BigN(1:NFourier*i,1:NFourier*i);
CenterCoord  = CenterCoord(:,1:i);
DiskRad      = DiskRad(1,1:i);
min_val_Data = min_val_Data(:,1:i+1);
end

%Useful code
% % If some ipolation values are on the bdry, we set them zero. If some
% % interpolation values are close to another zero, we are not correting
% % them or avoiding them whatsoever.
% ZetaX_idx= floor((BallInfos(1,1)-1)/size(ZetaX,1))+1;
% ZetaY_idx= mod((BallInfos(1,1)-1),size(ZetaX,1))+1;
% CenterIdx=(ZetaX_idx-1)*size(ZetaX,1)+ZetaY_idx;
% 
% CC = CenterIdx;
% SX1 = size(ZetaX,1);
% SX2 = size(ZetaX,2);
% %Cornercases
% if ZetaX_idx == 1 && ZetaY_idx == 1
%     IpSquare_idx = [0, 2, SX1+2, 0, 1, SX1+1, 0, 0, 0];
% elseif  ZetaX_idx == SX1 && ZetaY_idx == SX2
%     IpSquare_idx = [SX1*(SX2-1)-1, SX1*SX2-1, 0,   SX1*(SX2-1), SX1*SX2, 0,   0, 0, 0];
% elseif ZetaX_idx == 1 && ZetaY_idx == SX2  %Down left corner
%     IpSquare_idx = [0, 0, 0,    0, SX1, 2*SX1,    0, 0, 0];
% elseif  ZetaX_idx == SX1 && ZetaY_idx == 1 %Up right corner
%     IpSquare_idx = [SX1*(SX2-2)+2, SX1*(SX2-1)+2, 0,   SX1*(SX2-2)+1, SX1*(SX2-1)+1, 0,   0, 0, 0];
%     %Edgecases
% elseif  ZetaY_idx == 1 %Upper edge
%     IpSquare_idx = [CC+1-SX1, CC+1, CC+SX1+1,   CC-SX1, CC, CC+SX1,   0, 0, 0];
% elseif ZetaY_idx == SX1 %Lower edge
%     IpSquare_idx = [0, 0, 0,   CC-SX1, CC, CC+SX1,   CC+1-SX1, CC+1, CC+SX1+1];
% elseif  ZetaX_idx == 1 %left edge
%     IpSquare_idx = [0, CC+1, CC+1+SX1,   0, CC, CC+SX1,   0, CC-1, CC-1+SX1];
% elseif  ZetaX_idx == SX2 %right edge
%     IpSquare_idx = [CC+1-SX1, CC+1, 0,   CC-SX1, CC, 0,   CC-1-SX1, CC-1, 0];
%     %NoBdryCase
% else
%     IpSquare_idx = [CC+1-SX1, CC+1, CC+1+SX1,   CC-SX1, CC, CC+SX1,   CC-1-SX1, CC-1, CC-1+SX1];
% end

