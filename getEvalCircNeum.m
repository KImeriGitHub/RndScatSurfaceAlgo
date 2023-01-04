function [NVals, z1x0NVals, z2x0NVals, z0x1NVals, z0x2NVals, ...
    z1x1NVals, z2x1NVals, z1x2NVals, z2x2NVals, NValsEX] = getEvalCircNeum(z, x, BigN, r_center_Coords, r_end, NFourier,k)
%PRE:  z is a 2 x Nz matrix of real values
%      x is a 2 x Nx matrix of real values
%      BigN is a (NCenters*NFourier)^2 matrix
%      r_center_Coords is a 2 x NCenters real matrix
%      r_end is a NCenters array
%      NFourier is a pos integer. Ideally a power of 2
%      k is a pos real number
%POST: All matrices of size Nz x Nx
%Desc: Neumann function evaluation with source z and receiver x. 
%   Also including partial derivatives. NValsEx Excludes the Gamma function, but else is the same as NVals.

NCenters=length(r_end);
Nz = size(z,2);
Nx = size(x,2);

twopi=2*pi;
sigmaVal_nor=twopi/NFourier;

tPts_ip = linspace(-pi,pi,NFourier+1);
tPts_ip = tPts_ip(1:end-1);
nor=[cos(tPts_ip);sin(tPts_ip)];
DiscPts_x = zeros(1,NFourier*NCenters);
DiscPts_y = zeros(1,NFourier*NCenters);
integConst = zeros(1,NFourier*NCenters);
for n = 1:NCenters
    DiscPts_x((1:NFourier)+(n-1)*NFourier) = r_end(n)*cos(tPts_ip).'+ r_center_Coords(1,n);
    DiscPts_y((1:NFourier)+(n-1)*NFourier) = r_end(n)*sin(tPts_ip).'+ r_center_Coords(2,n);
    integConst((1:NFourier)+(n-1)*NFourier) = r_end(n)*sigmaVal_nor;
end

[GammaZetaCirc, GammaZeta1Circ, GammaZeta2Circ, GammaZetaCirc1, GammaZetaCirc2, ...
    GammaZeta1Circ1, GammaZeta1Circ2, GammaZeta2Circ1, GammaZeta2Circ2, ~] ...
    = getGamma([z, x], [DiscPts_x; DiscPts_y], k);

GammaZetaCirc   = GammaZetaCirc.';
GammaZeta1Circ  = GammaZeta1Circ.';
GammaZeta2Circ  = GammaZeta2Circ.';
GammaZetaCirc1  = GammaZetaCirc1.';
GammaZetaCirc2  = GammaZetaCirc2.';
GammaZeta1Circ1 = GammaZeta1Circ1.';
GammaZeta1Circ2 = GammaZeta1Circ2.';
GammaZeta2Circ1 = GammaZeta2Circ1.';
GammaZeta2Circ2 = GammaZeta2Circ2.';

GammaZetaDelCirc = GammaZetaCirc1.*repmat(nor(1,:),Nz+Nx,NCenters)+GammaZetaCirc2.*repmat(nor(2,:),Nz+Nx,NCenters);
GammaZeta1DelCirc = GammaZeta1Circ1.*repmat(nor(1,:),Nz+Nx,NCenters)+GammaZeta1Circ2.*repmat(nor(2,:),Nz+Nx,NCenters);
GammaZeta2DelCirc = GammaZeta2Circ1.*repmat(nor(1,:),Nz+Nx,NCenters)+GammaZeta2Circ2.*repmat(nor(2,:),Nz+Nx,NCenters);

Gamma_z_Circ = GammaZetaCirc(1:Nz,:);
%Gamma_x_Circ = GammaZetaCirc((1:Nx)+Nz,:);
%Gamma_z_Circ1 = GammaZetaCirc1(1:Nz,:);
%Gamma_x_Circ1 = GammaZetaCirc1((1:Nx)+Nz,:);
%Gamma_z_Circ2 = GammaZetaCirc2(1:Nz,:);
%Gamma_x_Circ2 = GammaZetaCirc2((1:Nx)+Nz,:);
Gamma_z1_Circ = GammaZeta1Circ(1:Nz,:);
%Gamma_x1_Circ = GammaZeta1Circ((1:Nx)+Nz,:);
Gamma_z2_Circ = GammaZeta2Circ(1:Nz,:);
% Gamma_x2_Circ = GammaZeta2Circ((1:Nx)+Nz,:);
% Gamma_z1_Circ1 = GammaZeta1Circ1(1:Nz,:);
% Gamma_x1_Circ1 = GammaZeta1Circ1((1:Nx)+Nz,:);
% Gamma_z2_Circ1 = GammaZeta2Circ1(1:Nz,:);
% Gamma_x2_Circ1 = GammaZeta2Circ1((1:Nx)+Nz,:);
% Gamma_z1_Circ2 = GammaZeta1Circ2(1:Nz,:);
% Gamma_x1_Circ2 = GammaZeta1Circ2((1:Nx)+Nz,:);
% Gamma_z2_Circ2 = GammaZeta2Circ2(1:Nz,:);
% Gamma_x2_Circ2 = GammaZeta2Circ2((1:Nx)+Nz,:);

Gamma_z_DelCirc = GammaZetaDelCirc(1:Nz,:).*repmat(integConst,Nz,1);
Gamma_x_DelCirc = GammaZetaDelCirc((1:Nx)+Nz,:).*repmat(integConst,Nx,1);
Gamma_z1_DelCirc = GammaZeta1DelCirc(1:Nz,:).*repmat(integConst,Nz,1);
Gamma_x1_DelCirc = GammaZeta1DelCirc((1:Nx)+Nz,:).*repmat(integConst,Nx,1);
Gamma_z2_DelCirc = GammaZeta2DelCirc(1:Nz,:).*repmat(integConst,Nz,1);
Gamma_x2_DelCirc = GammaZeta2DelCirc((1:Nx)+Nz,:).*repmat(integConst,Nx,1);

[Gamma_z_x, Gamma_z1_x, Gamma_z2_x, Gamma_z_x1, Gamma_z_x2, ...
    Gamma_z1_x1, Gamma_z1_x2, Gamma_z2_x1, Gamma_z2_x2, ~] ...
    = getGamma(z, x, k);

NVals = Gamma_z_x.' - Gamma_z_Circ*Gamma_x_DelCirc.' + Gamma_z_DelCirc * (2*BigN) * Gamma_x_DelCirc.';
z1x0NVals = Gamma_z1_x.' - Gamma_z1_Circ*Gamma_x_DelCirc.' + Gamma_z1_DelCirc * (2*BigN) * Gamma_x_DelCirc.';
z2x0NVals = Gamma_z2_x.' - Gamma_z2_Circ*Gamma_x_DelCirc.' + Gamma_z2_DelCirc * (2*BigN) * Gamma_x_DelCirc.';
z0x1NVals = Gamma_z_x1.' - Gamma_z_Circ*Gamma_x1_DelCirc.' + Gamma_z_DelCirc * (2*BigN) * Gamma_x1_DelCirc.';
z0x2NVals = Gamma_z_x2.' - Gamma_z_Circ*Gamma_x2_DelCirc.' + Gamma_z_DelCirc * (2*BigN) * Gamma_x2_DelCirc.';

z1x1NVals = Gamma_z1_x1.' - Gamma_z1_Circ*Gamma_x1_DelCirc.' + Gamma_z1_DelCirc * (2*BigN) * Gamma_x1_DelCirc.';
z2x1NVals = Gamma_z2_x1.' - Gamma_z2_Circ*Gamma_x1_DelCirc.' + Gamma_z2_DelCirc * (2*BigN) * Gamma_x1_DelCirc.';
z1x2NVals = Gamma_z1_x2.' - Gamma_z1_Circ*Gamma_x2_DelCirc.' + Gamma_z1_DelCirc * (2*BigN) * Gamma_x2_DelCirc.';
z2x2NVals = Gamma_z2_x2.' - Gamma_z2_Circ*Gamma_x2_DelCirc.' + Gamma_z2_DelCirc * (2*BigN) * Gamma_x2_DelCirc.';
NValsEX =Gamma_z_DelCirc * (2*BigN) * Gamma_x_DelCirc.';
end

