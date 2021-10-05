clear;
close all;
format compact

%% NOTE:
% Works very well. Very Stable (for big NFourier (>=2^6) and small k (<10)). 21 oct20.
% Works well with extra circle. However, error might not decrease anymore for too big NSteps (>50). 27Nov20.
% Does not yield O(dr^1) for k>5. (NF=2^8, r_end=1). Yields O(dr^2) for k=5
% with r_end=0.25, NF=2^8. Yields kinda O(dr^1) for k=5 with r_end=0.25, NF=2^7. 
% Yields log(error)=log(Delta r)^2 for k=20, r_end=0.25, NF=2^8. 
% Yields O(Delta r)^2 for k=20, r_end=0.1, NF=2^8. 08Apr21

% Check whether computing the big matrices is a time loss only due to the
% besselh functions
%%

NFourier=2^8; % Must be even! Preferably has low prime divisors.
k=25;
eg=0.57721566490153286060651209008240243104215933593992;
twopi=2*pi;
fourpi=4*pi;

NStepsDATA = 10:30;
rDeltaDATA = zeros(1, length(NStepsDATA));
errornorm = zeros(1, length(NStepsDATA));

r_init = 0.001;
r_end  = 0.5;

ZetaExt = [1;2.5];
rExt = 1;
DomExt = shape.Ellipse(rExt, rExt, NFourier)+ZetaExt;
theta_ip=linspace(0,twopi,NFourier+1);
theta_ip(end)=[];

Dom_Rend = shape.Ellipse(r_end, r_end, 1000)+[0;0];
Dom_end = merge(DomExt,Dom_Rend);
[~,~,~,~, RMat_end]= getNeumannFct(Dom_end, circshift(Dom_Rend.points,1000/2,2), [nan;nan], k);
RMat_end=RMat_end(NFourier+1:end,:);
RMat_end=circshift(RMat_end,1000/2);
CosM=cos((Dom_Rend.theta-pi)-(Dom_Rend.theta-pi).');
RMat_endp=RMat_end+(-1i/4*besselh(0,k*r_end*sqrt(2-2*CosM))-log(1-CosM)/fourpi);
for j=1:size(RMat_end,1)
    RMat_endp(j,j)=(RMat_end(j,j)+(2*eg-1i*pi+log(2)+2*log(k*r_end/2))/4/pi);
end
RMat_endp_ip=interp1(Dom_Rend.theta,RMat_endp,theta_ip);
RMat_endp_ip=interp1(Dom_Rend.theta,RMat_endp_ip.',theta_ip).';

for i = 1:length(NStepsDATA)
    disp(['i = ',num2str(i), ' of ', num2str(length(NStepsDATA))]);
    r_inc_Data = linspace(r_init,r_end,NStepsDATA(i));
    %r_inc_Data = (logspace(r_init,r_end,NSteps)-10^r_init)*(r_end-r_init)/(10^r_end-10^r_init)+r_init;
    
    RMat = getInflation(r_inc_Data, 1, [1;2.5], NFourier, k);
    
    rDeltaDATA(i) = r_inc_Data(2)-r_init;
    error=abs(RMat_endp_ip-RMat)./abs(RMat_endp_ip);
    errornorm(i) = sum(error,'all')/NFourier^2;
end


loglog(rDeltaDATA,(errornorm));hold on;
loglog(rDeltaDATA,rDeltaDATA.^2./rDeltaDATA(end).^2.*errornorm(end), '-.');hold on;
loglog(rDeltaDATA,rDeltaDATA.^1./rDeltaDATA(end).^1.*errornorm(end), '--');hold on;

title('Numerical Error');
xlabel('\Delta r');
ylabel('error');
legend('Actual Error','(\Delta r)^2', '(\Delta r)^1','Location','northwest')
% loglog(r_inc_DATA,r_inc_DATA.^2./r_inc_DATA(1).^2.*error_DATA(1));
% plot(r_inc_DATA,r_inc_DATA./r_inc_DATA(1).*error_DATA(1));
% title('error DATA');
%inspect=@(n,m)[term11(n,m),term12(n,m),term21(n,m),term22(n,m),term31(n,m),term32(n,m),term41(n,m),term42(n,m)];
%inspect2=@(n,m)[RMat(n,m),RInc(n,m),RMat_new(n,m)];
