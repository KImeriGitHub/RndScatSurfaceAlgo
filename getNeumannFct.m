function [res, BdryDens, IsInfIdx, resR,resRBDRY] = getNeumannFct(Bdry, Source, Receiver, k)
%Pre:  Bdry is a C2boundary class, can be empty
%      Source is a (2 x N) array, Source is not part of the boundary,
%      Receiver is a (2 x M) array, Receiver is not part of the boundary,
%      k is a real positive number
%Post: res is (M x N) array
%      BdryDens is a (size(Bdry.points,2) x N) array
%      IsInfIdx is a column of positive integers
%
% res represents the Neumann Function at the N source points.
% BdryDens represents the Neumann Function at boundary points
%     One column is associated to one Source point
%     Thus it is N^k(z,bdry.points).
% IsInfIdx show almost inf entries in res, i.e. res(IsInfIdx)) is large
% Source is not considered to be inside the domain.

twopi=2*pi;
eg=0.57721566490153286060651209008240243104215933593992;

M=size(Receiver,2);
N=size(Source,2);

SmR=(Source(1,:).'-Receiver(1,:)).^2 +(Source(2,:).'-Receiver(2,:)).^2;
IsInfIdx=find(SmR < 1e-14);

if isempty(Bdry)
    res = -1i/4*besselh(0,k*sqrt(SmR.'));
    logSmRplus=log(sqrt(SmR(mod(IsInfIdx,M*N)+1)));
    logSmRminus=log(sqrt(SmR(mod(IsInfIdx-2,M*N)+1)));
    logSmRplus(isinf(logSmRplus))=logSmRminus(isinf(logSmRplus));
    logSmRminus(isinf(logSmRminus))=logSmRplus(isinf(logSmRminus));
    logSmR = (logSmRplus+logSmRminus)/2;
    res(IsInfIdx) = (logSmR-2)/twopi+((log(k/2)+eg)/twopi-1i/4); 
    % We assign to inf number values which allows to integrate over a
    % logarithm integral. Only thought about the case in which Source ==
    % Receiver
    
    BdryDens = [];
    return
end

pts=Bdry.points;
tvec=Bdry.tvec;
avec=Bdry.avec;
normal=Bdry.normal;
sepPtr=Bdry.seperationPtr;
sigma=Bdry.sigma;

Npts=size(pts,2);

besselTemp = besselh(1,k*sqrt((pts(1,:).'-pts(1,:)).^2+((pts(2,:).'-pts(2,:)).^2)));
NPMat=1i*k/4*besselTemp...
        .*((pts(1,:)-pts(1,:).').*repmat(normal(1,:),Npts,1)+(pts(2,:)-pts(2,:).').*repmat(normal(2,:),Npts,1))...
        ./sqrt((pts(1,:).'-pts(1,:)).^2+((pts(2,:).'-pts(2,:)).^2))...
        .*repmat(sigma,Npts,1);
for i=1:sepPtr(end)-1
    NPMat(i,i) = -1/(2*twopi)*normal(:,i).'*avec(:,i)/sum(tvec(:,i).^2)*sigma(i);
end

LHSMat=0.5*eye(size(NPMat)) + NPMat;

SmP=sqrt((Source(1,:)-pts(1,:).').^2+(Source(2,:)-pts(2,:).').^2);
RhS=-1i/4*besselh(0,k*SmP);
logi_SmP_zero = SmP<1e-10;

BdryDens=zeros(size(SmP));
resRBDRY=zeros(size(SmP));
RhS_Mat=zeros(size(SmP));
RhS_Mat_G=zeros(size(SmP));
plus_Index=2:Npts+1; %We need 1+(1:Npts) but while counting inside the same boundary
for s=2:length(sepPtr)
    plus_Index(sepPtr(s)-1)=sepPtr(s)-2;
end
for i=1:N
    if sum(logi_SmP_zero(:,i))<0.5
        RhS_Mat(:,i)=RhS(:,i);
    else
        RhS_Vec_G=(-1i/4*besselh(0,k*SmP(:,i)));
        findfirst_SmP_zero=find(logi_SmP_zero(:,i),1);
        RhS_Vec_G(findfirst_SmP_zero)=(log(SmP(plus_Index(findfirst_SmP_zero),i))-2)/twopi+((log(k/2)+eg)/twopi-1i/4); 
        %The reason for this kind of smoothness is quite difficult to write
        %here. See one of the paper. -> Star4
        RhS_Mat(:,i)=-NPMat*RhS_Vec_G;
        RhS_Mat_G(:,i)=RhS_Vec_G;
    end
end

LHSslashRhS=LHSMat\RhS_Mat; % We extracted this operation from the loop because of speed up.

for i=1:N
    if sum(logi_SmP_zero(:,i))<0.5
        BdryDens(:,i)=LHSslashRhS(:,i);
    else
        resRBDRY(:,i)=LHSslashRhS(:,i);
        BdryDens(:,i)=LHSslashRhS(:,i) + RhS_Mat_G(:,i); % !!!Actually we shouldnt put the smoothed out RhS_Vec_G but the singular one. But it looks nicer this way.
    end
end


% FROM HERE ON OUT WE ASSUME THAT SOURCE IS NOT ON THE BOUNDARY

res0   = -1i/4*besselh(0,k*sqrt(SmR));
resRec = 1i*k/4*besselh(1,k*sqrt((Receiver(1,:).'-pts(1,:)).^2+(Receiver(2,:).'-pts(2,:)).^2))...
    .* (repmat(pts(1,:).*normal(1,:),M,1)-(Receiver(1,:).'.*normal(1,:))...
       +repmat(pts(2,:).*normal(2,:),M,1)-(Receiver(2,:).'.*normal(2,:)))...
    ./ sqrt((Receiver(1,:).'-pts(1,:)).^2 + (Receiver(2,:).'-pts(2,:)).^2)...
    .* repmat(sigma,M,1);

resR = (- BdryDens.'*resRec.').';
res  = res0.' + resR;


end

%TEST VECTORISATION
% A=zeros(2,3);
% Dom=shape.Ellipse(0.1,0.1,1000)+[0;0.3];
% Source=[0.1,-0.1;0.1,0.1];
% Receiver=[0.1,-0.1, 0.05;0.5,0.5, 0.6];
% 
% AFull=getNeumannFct(Dom, Source, Receiver, k);
% k=1.5;
% for i=1:size(A,1)
%     for j=1:size(A,2)
%         A(i,j) = getNeumannFct(Dom, Source(:,i), Receiver(:,j), k);
%     end
% end