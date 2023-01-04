function BigN = getCircNeumann(r_center_Coords, r_init, r_end, r_amount, NFourier,k)
%PRE:  r_center_Coords is a 2 x NCenters real matrix
%      r_init, r_end are a real numbers
%      r_amount, NFourier is an pos integer
%      k is a pos real number
%Desc: Returns a (NFourier*NCenters)^2 matrix. Represents the Neumann Function
% NFourier must be even! Preferably has low prime divisors.

%% Notes
% We do the SMat for the sources. Additionally we compute the functions at
% the rec pts. We also need to evaluate the function at the 9 points
% surrounding the circle

%% ASSERTS
assert(all(r_init < r_end));
assert(length(r_init) == 1);
assert(length(r_end) == size(r_center_Coords,2));

%% Changeable configuration
r_end_unique = unique(r_end);
r_end_min = r_end_unique(1);
r_inc_Data = [(logspace(r_init,r_end_min,r_amount)-10^r_init)*(r_end_min-r_init)/(10^r_end_min-10^r_init)+r_init, r_end_unique(2:end)];

%% Vital Variables
eg=0.57721566490153286060651209008240243104215933593992;
twopi=2*pi;
fourpi=4*pi;
sigmaVal_nor=twopi/NFourier;

rcCx = r_center_Coords(1,:).';
rcCy = r_center_Coords(2,:).';
% surround_Coords_x = [rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.' 
%                      rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.' 
%                      rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.' 
%                      rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.'
%                      rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.'];
% surround_Coords_y = [rcCy-r_end.',   rcCy-r_end.',   rcCy-r_end.',   rcCy-r_end.',   rcCy-r_end.'
%                      rcCy-r_end.'/2, rcCy-r_end.'/2, rcCy-r_end.'/2, rcCy-r_end.'/2, rcCy-r_end.'/2
%                      rcCy,           rcCy,           rcCy,           rcCy,           rcCy,
%                      rcCy+r_end.'/2, rcCy+r_end.'/2, rcCy+r_end.'/2, rcCy+r_end.'/2, rcCy+r_end.'/2
%                      rcCy+r_end.',   rcCy+r_end.',   rcCy+r_end.',   rcCy+r_end.',   rcCy+r_end.'];
surround_Coords_x = [rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.',...
                     rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.',...
                     rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.',...
                     rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.',...
                     rcCx-r_end.', rcCx-r_end.'/2, rcCx, rcCx+r_end.'/2, rcCx+r_end.'];
surround_Coords_y = [rcCy-r_end.',   rcCy-r_end.',   rcCy-r_end.',   rcCy-r_end.',   rcCy-r_end.',...
                     rcCy-r_end.'/2, rcCy-r_end.'/2, rcCy-r_end.'/2, rcCy-r_end.'/2, rcCy-r_end.'/2,...
                     rcCy,           rcCy,           rcCy,           rcCy,           rcCy,...
                     rcCy+r_end.'/2, rcCy+r_end.'/2, rcCy+r_end.'/2, rcCy+r_end.'/2, rcCy+r_end.'/2,...
                     rcCy+r_end.',   rcCy+r_end.',   rcCy+r_end.',   rcCy+r_end.',   rcCy+r_end.'];

NCenters=size(r_center_Coords,2);

integConst = zeros(1,NFourier*NCenters); % BIG
for n = 1:NCenters
    integConst((1:NFourier)+(n-1)*NFourier) = r_end(n)*sigmaVal_nor;
end

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

% Extra variables
tPts_ip = linspace(-pi,pi,NFourier+1);
tPts_ip = tPts_ip(1:end-1);
dtPts = tPts_ip(2)-tPts_ip(1);
nor=[cos(tPts_ip);sin(tPts_ip)];

IpVals_R=zeros(5,5,5,5,7);
%N_MAT = zeros([NFourier,NFourier,NCenters,NCenters]);
BigN = zeros([NFourier*NCenters,NFourier*NCenters]); 

ALL_DiscPts_x = zeros(NFourier, NCenters); %BIG
ALL_DiscPts_y = zeros(NFourier, NCenters); %BIG
for n = 1:NCenters
    ALL_DiscPts_x(:,n) = r_end(n)*cos(tPts_ip).'+ r_center_Coords(1,n);
    ALL_DiscPts_y(:,n) = r_end(n)*sin(tPts_ip).'+ r_center_Coords(2,n);
end

CosM = cos((tPts_ip)-(tPts_ip).');
LMat = log(1-CosM)/fourpi;
for i = 1:NFourier
    LMat(i,i) = (log(tPts_ip(2)-tPts_ip(1))-2)/twopi;
end
% LMat(1,:)*ones(NFourier,1)*(tPts_ip(2)-tPts_ip(1)) should give -0.346574


%% First Step
circleDATA_Sin=zeros([NFourier, NFourier, NCenters]); %TOO BIG
circleDATA_Cos=zeros([NFourier, NFourier, NCenters]); %TOO BIG
circleDATA_Cos(:,:,1) = [(log(k*r_init)/4/pi + eg/4/pi - 1i/8 + 0/2)*ones(NFourier,1), zeros(NFourier,NFourier-1)];

for ii = 2:length(r_inc_Data)
    %Info: 'n' means 'new' and refers to the update
    % 'R' is only for circle values and refers to Pts at circle with
    % radius R but domain with circle r

    % Actual increase
    r=r_inc_Data(ii-1);
    R=r_inc_Data(ii);
    [circleDATA_Cos(:,:,1), circleDATA_Sin(:,:,1)] = getDiscInfl(circleDATA_Cos(:,:,1), circleDATA_Sin(:,:,1),...
        IpVals_R, r, R, k, r_end(1), DDG2_FCoeffs_PREC(ii,:), DDG12_FCoeffs_PREC(ii,:));
    RMatn =  circleDATA_Cos(:,:,1) * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]...%'n' stand for new
           + circleDATA_Sin(:,:,1) * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
    
    if abs(r_inc_Data(ii) - r_end(1) ) < 2*eps
        break;
    end
end
%N_MAT(:,:,1,1) = RMatn+LMat;
BigN(1:NFourier,1:NFourier) = RMatn+LMat;

%% RETURN CONDITION
if NCenters == 1
    return 
end

%% R Vals

surround_Coords = [surround_Coords_x(2,:); surround_Coords_y(2,:)];
DiscPts = [ALL_DiscPts_x(:,1),ALL_DiscPts_y(:,1)].';
[GammaZetaCirct, GammaZeta1Circt, GammaZeta2Circt, GammaZetaCirc1t, GammaZetaCirc2t, ...
    GammaZeta1Circ1t, GammaZeta1Circ2t, GammaZeta2Circ1t, GammaZeta2Circ2t, ~] ...
    = getGamma(surround_Coords, DiscPts, k);

[RVals, ~, ~, nabla0nabla1RVals, nabla0nabla2RVals, ...
    nabla1nabla1RVals, nabla1nabla2RVals, nabla2nabla1RVals, nabla2nabla2RVals] ...
    = getEvalCircReum(surround_Coords, surround_Coords, BigN(1:NFourier,1:NFourier), r_center_Coords(:,1:1), r_end(:,1:1), NFourier,k);

GammaZetaDelCirc  = GammaZetaCirc1t.'.*repmat(nor(1,:),25,1)+GammaZetaCirc2t.'.*repmat(nor(2,:),25,1);
GammaZeta1DelCirc = GammaZeta1Circ1t.'.*repmat(nor(1,:),25,1)+GammaZeta1Circ2t.'.*repmat(nor(2,:),25,1);
GammaZeta2DelCirc = GammaZeta2Circ1t.'.*repmat(nor(1,:),25,1)+GammaZeta2Circ2t.'.*repmat(nor(2,:),25,1);

    
% This is actually very difficult logistical problem
IpVals_R(:,:,:,:,1) = reshape(RVals,5,5,5,5);
IpVals_R(:,:,:,:,2) = reshape(nabla0nabla1RVals,5,5,5,5);
IpVals_R(:,:,:,:,3) = reshape(nabla0nabla2RVals,5,5,5,5);
IpVals_R(:,:,:,:,4) = reshape(nabla1nabla1RVals,5,5,5,5);
IpVals_R(:,:,:,:,5) = reshape(nabla1nabla2RVals,5,5,5,5);
IpVals_R(:,:,:,:,6) = reshape(nabla2nabla1RVals,5,5,5,5);
IpVals_R(:,:,:,:,7) = reshape(nabla2nabla2RVals,5,5,5,5);
RCC = RVals(13,13);

%% Interpolated N0(Old Circle, New Circle) 
N0CircZeta  = GammaZetaCirct -   (2*BigN(1:NFourier,1:NFourier))*(GammaZetaDelCirc.' *(sigmaVal_nor*r_end(1)));
N0CircZeta1 = GammaZeta1Circt -   (2*BigN(1:NFourier,1:NFourier))*(GammaZeta1DelCirc.' *(sigmaVal_nor*r_end(1)));
N0CircZeta2 = GammaZeta2Circt -   (2*BigN(1:NFourier,1:NFourier))*(GammaZeta2DelCirc.' *(sigmaVal_nor*r_end(1)));

N0CircZetaR  = zeros(NFourier,NFourier);
N0CircZetaR1 = zeros(NFourier,NFourier);
N0CircZetaR2 = zeros(NFourier,NFourier);
for n=1:NFourier
  N0CircZetaR(n,:)  = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircZeta(n,:) ,5,5),cos(tPts_ip),sin(tPts_ip),'makima');
  N0CircZetaR1(n,:) = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircZeta1(n,:),5,5),cos(tPts_ip),sin(tPts_ip),'makima');
  N0CircZetaR2(n,:) = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircZeta2(n,:),5,5),cos(tPts_ip),sin(tPts_ip),'makima');
end

N0CircDelZetaR = N0CircZetaR1.*repmat(nor(1,:),NFourier,1)+N0CircZetaR2.*repmat(nor(2,:),NFourier,1);

%% Main Loop
for i=2:(NCenters-1)
    disp(['i = ',num2str(i),' of ',num2str(NCenters)]);
    %% Initial New Radii
    circleDATA_Sin(:,:,i) = zeros([NFourier, NFourier]);
    circleDATA_Cos(:,:,i) = [(log(k*r_init)/4/pi + 0.57721566490153286060651209/4/pi - 1i/8 + RCC/2)...
            *ones(NFourier,1), zeros(NFourier,NFourier-1)];
    
    for ii = 2:length(r_inc_Data)
        %Info: 'n' means 'new' and refers to the update
        % 'R' is only for circle values and refers to Pts at circle with
        % radius R but domain with circle r
        
        % Actual increase
        r=r_inc_Data(ii-1);
        R=r_inc_Data(ii);
        [circleDATA_Cos(:,:,i), circleDATA_Sin(:,:,i)] = getDiscInfl(circleDATA_Cos(:,:,i), circleDATA_Sin(:,:,i),...
            IpVals_R, r, R, k, r_end(i), DDG2_FCoeffs_PREC(ii,:), DDG12_FCoeffs_PREC(ii,:));
        RMatn =  circleDATA_Cos(:,:,i) * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]... %'n' stand for new
               + circleDATA_Sin(:,:,i) * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
        
        if abs(R - r_end(i)) < 2*eps
            break;
        end
    end
    %N_MAT(:,:,i,i) = RMatn+LMat;
    BigN((1:NFourier)+(i-1)*NFourier,(1:NFourier)+(i-1)*NFourier) = RMatn+LMat;
    
    %% Update old circles
    for j=1:i-1
        temp = -0.5*(N0CircZetaR((1:NFourier) + NFourier*(j-1),:))*(N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:).' * (sigmaVal_nor*r_end(i)))...
               +0.5*(N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:)*(2*BigN((1:NFourier)+(i-1)*NFourier,(1:NFourier)+(i-1)*NFourier))*N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:).') *(sigmaVal_nor*r_end(i))^2;
        circleDATA_Cos(:,:,j) = circleDATA_Cos(:,:,j) + [getFourierConst(temp),zeros(NFourier,NFourier/2)];
        RMatn =  circleDATA_Cos(:,:,j) * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]...%'n' stand for new
               + circleDATA_Sin(:,:,j) * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
        %N_MAT(:,:,j,j) = RMatn+LMat;
        BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(j-1)*NFourier) = RMatn+LMat;
    end
    
    %% New CrossActions
    for j=1:(i-1)
        %N_MAT(:,:,j,i) = N_MAT(:,:,j,i) ...
        %    - (N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:))*(N_MAT(:,:,i,i).' *(sigmaVal_nor*r_end(i)));
        BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(i-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(i-1)*NFourier)...
            - (N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:))*(BigN((1:NFourier)+(i-1)*NFourier,(1:NFourier)+(i-1)*NFourier).' *(sigmaVal_nor*r_end(i)));
        %N_MAT(:,:,i,j) = N_MAT(:,:,j,i).';
        BigN((1:NFourier)+(i-1)*NFourier,(1:NFourier)+(j-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(i-1)*NFourier).';
    end
    
    %% Update Old CrossActions
    for j=1:(i-1)
        for jj=1:(j-1)
            %N_MAT(:,:,j,jj) = N_MAT(:,:,j,jj) ...
            %    - (N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:))*(N_MAT(:,:,jj,i).'*(sigmaVal_nor*r_end(i)));
            BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier)...
                - (N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:))*(BigN((1:NFourier)+(jj-1)*NFourier,(1:NFourier)+(i-1)*NFourier).'*(sigmaVal_nor*r_end(i)));
            %N_MAT(:,:,jj,j) = N_MAT(:,:,j,jj).';
            BigN((1:NFourier)+(jj-1)*NFourier,(1:NFourier)+(j-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier).';
        end
    end
    
    %% Get R Values  
%     BigN = zeros(NFourier*i, NFourier*i);
%     for ii =1:i
%         for iii =1:i
%             BigN((1:NFourier)+NFourier*(ii-1),(1:NFourier)+NFourier*(iii-1)) = N_MAT(:,:,ii,iii);
%         end
%     end
    surround_Coords = [surround_Coords_x(i+1,:); surround_Coords_y(i+1,:)];
    DiscPts = [reshape(ALL_DiscPts_x(:,1:i),NFourier*i,1),reshape(ALL_DiscPts_y(:,1:i),NFourier*i,1)].';
    [GammaZetaCirct, GammaZeta1Circt, GammaZeta2Circt, GammaZetaCirc1t, GammaZetaCirc2t, ...
        GammaZeta1Circ1t, GammaZeta1Circ2t, GammaZeta2Circ1t, GammaZeta2Circ2t, ~] ...
        = getGamma(surround_Coords, DiscPts, k);
    
    [RVals, ~, ~, nabla0nabla1RVals, nabla0nabla2RVals, ...
        nabla1nabla1RVals, nabla1nabla2RVals, nabla2nabla1RVals, nabla2nabla2RVals] ...
        = getEvalCircReum(surround_Coords, surround_Coords, BigN(1:(i)*NFourier,1:(i)*NFourier), r_center_Coords(:,1:i), r_end(1:i), NFourier,k);
            
    GammaZetaDelCirc = GammaZetaCirc1t.'.*repmat(nor(1,:),25,i)+GammaZetaCirc2t.'.*repmat(nor(2,:),25,i);
    GammaZeta1DelCirc = GammaZeta1Circ1t.'.*repmat(nor(1,:),25,i)+GammaZeta1Circ2t.'.*repmat(nor(2,:),25,i);
    GammaZeta2DelCirc = GammaZeta2Circ1t.'.*repmat(nor(1,:),25,i)+GammaZeta2Circ2t.'.*repmat(nor(2,:),25,i);        
    
    
    % This is actually a very difficult logistical problem
    IpVals_R(:,:,:,:,1) = reshape(RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,2) = reshape(nabla0nabla1RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,3) = reshape(nabla0nabla2RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,4) = reshape(nabla1nabla1RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,5) = reshape(nabla1nabla2RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,6) = reshape(nabla2nabla1RVals.',5,5,5,5);
    IpVals_R(:,:,:,:,7) = reshape(nabla2nabla2RVals.',5,5,5,5);
    RCC = RVals(13,13);
    
    % Interpolated N0(Old Circle, New Circle) 
    N0CircZeta  = (GammaZetaCirct.'  -  GammaZetaDelCirc.*repmat(integConst(1, 1:i*NFourier),25,1)*(BigN(1:NFourier*i,1:NFourier*i)*2)).';
    N0CircZeta1 = (GammaZeta1Circt.'  -  GammaZeta1DelCirc.*repmat(integConst(1, 1:i*NFourier),25,1)*(BigN(1:NFourier*i,1:NFourier*i)*2)).';
    N0CircZeta2 = (GammaZeta2Circt.'  -  GammaZeta2DelCirc.*repmat(integConst(1, 1:i*NFourier),25,1)*(BigN(1:NFourier*i,1:NFourier*i)*2)).';  
    
    N0CircZetaR  = zeros(NFourier*i, NFourier);
    N0CircZetaR1 = zeros(NFourier*i, NFourier);
    N0CircZetaR2 = zeros(NFourier*i, NFourier);
    % Reshaping
    reshN0CircZeta=zeros(5,5,NFourier*i);
    reshN0CircZeta1=zeros(5,5,NFourier*i);
    reshN0CircZeta2=zeros(5,5,NFourier*i);
    for n=1:NFourier*i
        reshN0CircZeta(:,:,n)=reshape(N0CircZeta(n,:) ,5,5);
        reshN0CircZeta1(:,:,n)=reshape(N0CircZeta1(n,:) ,5,5);
        reshN0CircZeta2(:,:,n)=reshape(N0CircZeta2(n,:) ,5,5);
    end
    tempIV = [-1,-0.5,0,0.5,1];
    N0CircZetaRip = interpn(tempIV,tempIV, linspace(0,1,NFourier*i), reshN0CircZeta, repmat(cos(tPts_ip),1,NFourier*i),repmat(sin(tPts_ip),1,NFourier*i), kron(linspace(0,1,NFourier*i),ones(1,NFourier)),'makima');
    N0CircZetaR1ip = interpn(tempIV,tempIV, linspace(0,1,NFourier*i), reshN0CircZeta1, repmat(cos(tPts_ip),1,NFourier*i),repmat(sin(tPts_ip),1,NFourier*i), kron(linspace(0,1,NFourier*i),ones(1,NFourier)),'makima');
    N0CircZetaR2ip = interpn(tempIV,tempIV, linspace(0,1,NFourier*i), reshN0CircZeta2, repmat(cos(tPts_ip),1,NFourier*i),repmat(sin(tPts_ip),1,NFourier*i), kron(linspace(0,1,NFourier*i),ones(1,NFourier)),'makima');
    
    % BackShaping
    for n=1:NFourier*i
        N0CircZetaR(n,:) = N0CircZetaRip((1:NFourier) + (n-1)*NFourier);
        N0CircZetaR1(n,:) = N0CircZetaR1ip((1:NFourier) + (n-1)*NFourier);
        N0CircZetaR2(n,:) = N0CircZetaR2ip((1:NFourier) + (n-1)*NFourier);
    end
    
%     %THIS IS THE ORIGINAL FORM; Withoutreshaping; very slow
%     for n=1:NFourier*i
%         N0CircZetaRp(n,:)  = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircZeta(n,:) ,5,5),cos(tPts_ip),sin(tPts_ip),'spline');
%         N0CircZetaR1p(n,:) = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircZeta1(n,:),5,5),cos(tPts_ip),sin(tPts_ip),'spline');
%         N0CircZetaR2p(n,:) = interpn([-1,-0.5,0,0.5,1],[-1,-0.5,0,0.5,1], reshape(N0CircZeta2(n,:),5,5),cos(tPts_ip),sin(tPts_ip),'spline');
%     end
    
    N0CircDelZetaR = N0CircZetaR1.*repmat(nor(1,:),NFourier*i,1)+N0CircZetaR2.*repmat(nor(2,:),NFourier*i,1);

end

%% End Step
circleDATA_Sin(:,:,NCenters) = zeros([NFourier, NFourier]);
circleDATA_Cos(:,:,NCenters) = [(log(k*r_init)/4/pi + 0.57721566490153286060651209/4/pi - 1i/8 + RCC/2)...
    *ones(NFourier,1), zeros(NFourier,NFourier-1)];

for ii = 2:length(r_inc_Data)
    %Info: 'n' means 'new' and refers to the update
    % 'R' is only for circle values and refers to Pts at circle with
    % radius R but domain with circle r
    
    % Actual increase
    r=r_inc_Data(ii-1);
    R=r_inc_Data(ii);
    [circleDATA_Cos(:,:,NCenters), circleDATA_Sin(:,:,NCenters)] = getDiscInfl(circleDATA_Cos(:,:,NCenters), circleDATA_Sin(:,:,NCenters),...
        IpVals_R, r, R, k, r_end(NCenters), DDG2_FCoeffs_PREC(ii,:), DDG12_FCoeffs_PREC(ii,:));
    RMatn =  circleDATA_Cos(:,:,NCenters) * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]...%'n' stand for new
        + circleDATA_Sin(:,:,NCenters) * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
    
    if abs(R - r_end(NCenters)) < 2*eps
        break;
    end
end
%N_MAT(:,:,NCenters,NCenters) = RMatn+LMat;
BigN((1:NFourier)+NFourier*(NCenters-1),(1:NFourier)+NFourier*(NCenters-1)) = RMatn+LMat;

% Update old circles
for j=1:NCenters-1
    temp = -0.5*(N0CircZetaR((1:NFourier) + NFourier*(j-1),:))*(N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:).' *(sigmaVal_nor*r_end(NCenters)))...
            +0.5*(N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:)*(2*BigN((1:NFourier)+NFourier*(NCenters-1),(1:NFourier)+NFourier*(NCenters-1)))*N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:).') *(sigmaVal_nor*r_end(NCenters))^2;
    circleDATA_Cos(:,:,j) = circleDATA_Cos(:,:,j) + [getFourierConst(temp),zeros(NFourier,NFourier/2)];
    RMatn =  circleDATA_Cos(:,:,j) * [ones(1,NFourier) ;cos((1:NFourier-1).'*tPts_ip)]...%'n' stand for new
        + circleDATA_Sin(:,:,j) * [zeros(1,NFourier);sin((1:NFourier-1).'*tPts_ip)];
    %N_MAT(:,:,j,j) = RMatn+LMat;
    BigN((1:NFourier)+NFourier*(j-1),(1:NFourier)+NFourier*(j-1)) = RMatn+LMat;
end

% New CrossActions
for j=1:(NCenters-1)
    %N_MAT(:,:,j,NCenters) = N_MAT(:,:,j,NCenters) - (N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:))*(N_MAT(:,:,NCenters,NCenters).' *(sigmaVal_nor*r_end(NCenters)));
    BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(NCenters-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(NCenters-1)*NFourier)...
            - (N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:))*(BigN((1:NFourier)+(NCenters-1)*NFourier,(1:NFourier)+(NCenters-1)*NFourier).' *(sigmaVal_nor*r_end(NCenters)));
    %N_MAT(:,:,NCenters,j) = N_MAT(:,:,j,NCenters).';
    BigN((1:NFourier)+(NCenters-1)*NFourier,(1:NFourier)+(j-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(NCenters-1)*NFourier).';
end

% Update Old CrossActions
for j=1:NCenters-1
    for jj=1:j-1
        %N_MAT(:,:,j,jj) = N_MAT(:,:,j,jj) - (N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:))*(N_MAT(:,:,jj,NCenters).'*(sigmaVal_nor*r_end(NCenters)));
        BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier)...
                - (N0CircDelZetaR((1:NFourier) + NFourier*(j-1),:))*(BigN((1:NFourier)+(jj-1)*NFourier,(1:NFourier)+(NCenters-1)*NFourier).'*(sigmaVal_nor*r_end(NCenters)));
        %N_MAT(:,:,jj,j) = N_MAT(:,:,j,jj).';
        BigN((1:NFourier)+(jj-1)*NFourier,(1:NFourier)+(j-1)*NFourier) = BigN((1:NFourier)+(j-1)*NFourier,(1:NFourier)+(jj-1)*NFourier).';
    end
end

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

