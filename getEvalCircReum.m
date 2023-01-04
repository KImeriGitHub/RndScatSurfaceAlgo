function [RVals, z1x0RVals, z2x0RVals, z0x1RVals, z0x2RVals, ...
    z1x1RVals, z2x1RVals, z1x2RVals, z2x2RVals] = getEvalCircReum(z, x, BigN, r_center_Coords, r_end, NFourier,k)
%DESC: Same as getEvalCircNeum but only evaluating the remainder function.

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

GammaZetaDelCirc  = GammaZetaCirc1.*repmat(nor(1,:),Nz+Nx,NCenters)+GammaZetaCirc2.*repmat(nor(2,:),Nz+Nx,NCenters);
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

% Gamma_z_x = getNeumannFct([], z, x, k).';
% 
% Gamma_z1_x = getNeumannFct([], z+[discDis;0], x, k).';
% Gamma_z2_x = getNeumannFct([], z+[0;discDis], x, k).';
% Gamma_z_x1 = getNeumannFct([], z, x+[discDis;0], k).';
% Gamma_z_x2 = getNeumannFct([], z, x+[0;discDis], k).';

% Gamma_z1_x1 = getNeumannFct([], z+[discDis;0], x+[discDis;0], k).';
% Gamma_z2_x1 = getNeumannFct([], z+[0;discDis], x+[discDis;0], k).';
% Gamma_z1_x2 = getNeumannFct([], z+[discDis;0], x+[0;discDis], k).';
% Gamma_z2_x2 = getNeumannFct([], z+[0;discDis], x+[0;discDis], k).';
% 
% Gamma_z1_x = (Gamma_z1_x-Gamma_z_x)/discDis;
% Gamma_z2_x = (Gamma_z2_x-Gamma_z_x)/discDis;
% Gamma_z_x1 = (Gamma_z_x1-Gamma_z_x)/discDis;
% Gamma_z_x2 = (Gamma_z_x2-Gamma_z_x)/discDis;

% Gamma_z1_x1 = (Gamma_z1_x1 + Gamma_z1_x + Gamma_z_x1 - Gamma_z_x)/discDis^2;
% Gamma_z2_x1 = (Gamma_z2_x1 + Gamma_z2_x + Gamma_z_x1 - Gamma_z_x)/discDis^2;
% Gamma_z1_x2 = (Gamma_z1_x2 + Gamma_z1_x + Gamma_z_x2 - Gamma_z_x)/discDis^2;
% Gamma_z2_x2 = (Gamma_z2_x2 + Gamma_z2_x + Gamma_z_x2 - Gamma_z_x)/discDis^2;

RVals =  - Gamma_z_Circ*Gamma_x_DelCirc.' + Gamma_z_DelCirc * (2*BigN) * Gamma_x_DelCirc.';
z1x0RVals =  - Gamma_z1_Circ*Gamma_x_DelCirc.' + Gamma_z1_DelCirc * (2*BigN) * Gamma_x_DelCirc.';
z2x0RVals =  - Gamma_z2_Circ*Gamma_x_DelCirc.' + Gamma_z2_DelCirc * (2*BigN) * Gamma_x_DelCirc.';
z0x1RVals =  - Gamma_z_Circ*Gamma_x1_DelCirc.' + Gamma_z_DelCirc * (2*BigN) * Gamma_x1_DelCirc.';
z0x2RVals =  - Gamma_z_Circ*Gamma_x2_DelCirc.' + Gamma_z_DelCirc * (2*BigN) * Gamma_x2_DelCirc.';

z1x1RVals =  - Gamma_z1_Circ*Gamma_x1_DelCirc.' + Gamma_z1_DelCirc * (2*BigN) * Gamma_x1_DelCirc.';
z2x1RVals =  - Gamma_z2_Circ*Gamma_x1_DelCirc.' + Gamma_z2_DelCirc * (2*BigN) * Gamma_x1_DelCirc.';
z1x2RVals =  - Gamma_z1_Circ*Gamma_x2_DelCirc.' + Gamma_z1_DelCirc * (2*BigN) * Gamma_x2_DelCirc.';
z2x2RVals =  - Gamma_z2_Circ*Gamma_x2_DelCirc.' + Gamma_z2_DelCirc * (2*BigN) * Gamma_x2_DelCirc.';

end

