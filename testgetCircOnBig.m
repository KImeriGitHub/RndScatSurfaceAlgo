clear

k=15;
NFourier = 2^7;
%parpool(24)

%                               |          *
%                               |       *  S  *
%                               |     *    *    *
%                               |       *     *
%                               |     *    *    *
%---------------------------------------*-----*---------------------------
%                               |     *    *    *
%                               |       *     *
%                               |          *
%                               |

% r_center_Coords_Data   = [0.5 0.5  0.5  0.2  0.8  0.2  0.8  0.2  0.8  0.2  0.8 -0.1 1.1 -0.1  1.1 -0.1 1.1 0.5  0.5;
%                           0   0.5 -0.5  0.3 -0.3 -0.3  0.3  0.8 -0.8 -0.8  0.8  0   0    0.6 -0.6 -0.6 0.6 1.2 -1.2];  
% z=[0.5;  0.9];
% x=[0.5; -0.9];
% 
% z=[0.5;  0.9];
% x=[0.5; -0.9];
% r_init = 0.001;
% 
% r_end_vals = linspace(0.01,0.03,10);
% r_end_Data = r_end_vals(randi(10,1,size(r_center_Coords_Data,2)));
% r_amount = 10;
% 
% xPlot=linspace(-1,2,200);
% nx=length(xPlot);
% yPlot=linspace(-2,2,nx);

%                               |         
%                               |       
%                               |     
%           * * * * * * * * * * * * * * * * * * * * * *
%                               |    
%-------------------------------S-----------------------------------------
%                               |     
%           * * * * * * * * * * * * * * * * * * * * * *                    
%                               |   
%                               |
%r_center_Coords_Data   = [-0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1  0    0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1  0    0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9;
%                           0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2 -0.2];                        
% r_center_Coords_Data = [-3:0.1:3,-3:0.1:3; 0.15*ones(1,length(-3:0.1:3)),-0.15*ones(1,length(-3:0.1:3))];
%                      
% z=[0;  0];
% x=[0.5; -0.9];
% 
% xPlot=linspace(-5,5,300);
% nx=length(xPlot);
% yPlot=linspace(-1.2,1.2,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                               |   *       
%         *                     |                  *
%                     *         |     
%                               |          *           
%              *                |    
%-------------------------------S----------------------------
%                           *   |      *
%                               |                                     
%      *                        |                    *
%                               |*

NDisks=20;
nd_loop = 1;
r_center_Coords_Data = zeros(2,NDisks);
while nd_loop <= NDisks
    newCoord = (rand(2,1)-0.5).*2; % Vector in [-1,1]^2
    if ~any(sum((r_center_Coords_Data-newCoord).^2)<0.12^2) % if newCoord are close to a existing disk or the origin
        r_center_Coords_Data(:,nd_loop) = newCoord;
        nd_loop=nd_loop+1;
    end
end

r_init = 0.001;
r_end_vals = linspace(0.01,0.05,15);
r_end_Data = r_end_vals(randi(length(r_end_vals),1,size(r_center_Coords_Data,2)));
r_amount = 10;

z=[-1.1;  0];
x=[1.1;  0];

xPlot=linspace(-1.5,1.5,200);
nx=length(xPlot);
yPlot=linspace(-1.5,1.5,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r_init = 0.001;
% r_end_Data = 0.03*ones(1,size(r_center_Coords_Data,2));
% r_amount = 20;
NCenter = size(r_center_Coords_Data,2);

tPts_Full = linspace(-pi,pi,NFourier+1);
tPts_Full = tPts_Full(1:end-1);

tPts_semi = linspace(-pi,pi,80+1);
tPts_semi = tPts_semi(1:end-1);

r_center_Coords = r_center_Coords_Data;
r_end = r_end_Data;

sepPtr=[1]; D=[];
for i=1:length(r_end)
    D=[D, r_end(i)*[cos(tPts_Full);sin(tPts_Full)]+r_center_Coords(:,i)];
    sepPtr=[sepPtr, size(D,2)+1];
end
Dom = shape.Spline(D, sepPtr, 3);

sepPtr=[1]; D=[];
for i=1:length(r_end)
    D=[D, r_end(i)*[cos(tPts_semi);sin(tPts_semi)]+r_center_Coords(:,i)];
    sepPtr=[sepPtr, size(D,2)+1];
end
Dom_semi = shape.Spline(D, sepPtr, 3);

tic;
NExact_Data = getNeumannFct(Dom_semi, z, x, k).';
NExactTime_Data = toc;

tic;
%try
    BigN = getCircNeumann(r_center_Coords, r_init, r_end, r_amount, NFourier, k);
%catch
%    warning('Using Inflation method did not work')
%end
[N_Data, ~, ~, ~, ~, ~, ~, ~, ~] = getEvalCircNeum(z, x, BigN, r_center_Coords, r_end, NFourier,k);
NApTime_Data = toc;

%colormap('hot'); imagesc(real(BigN)); colorbar;

%% Plots
[X,Y]=meshgrid(xPlot,yPlot);
Z=[reshape(X,1,nx^2);reshape(Y,1,nx^2)];

% RAM saving
NumSplits = max([floor((nx/100)^2*(NFourier/128)^2*(NCenter/16)^2),1]);
Splits = floor(linspace(1,nx^2+1,floor(NumSplits)+1));
NApZ = zeros(nx^2,1);
%NExZ = zeros(nx^2,1);
ZxSplits = cell(NumSplits, 1);
ZySplits = cell(NumSplits, 1);
NApZSplits = cell(NumSplits, 1);
for i=1:NumSplits
    ZxSplits(i)={Z(1,Splits(i):Splits(i+1)-1)};
    ZySplits(i)={Z(2,Splits(i):Splits(i+1)-1)};
end

for i=1:NumSplits
    disp(['i = ',num2str(i),' of ',num2str(NumSplits)]);
    Zi=[ZxSplits{i}; ZySplits{i}];
    %NExZ(Splits(i):Splits(i+1)-1,1) = getNeumannFct(Dom_semi,z,Zi,k);
    [a, ~, ~, ~, ~, ~, ~, ~, ~] = getEvalCircNeum(z, Zi, BigN, r_center_Coords, r_end, NFourier,k);
    NApZSplits(i) = {a};
end
for i=1:NumSplits
    NApZ(Splits(i):Splits(i+1)-1,1)= NApZSplits{i};
end
NExZ = getNeumannFct(Dom_semi,z,Z,k);

NExGrid=reshape(NExZ, nx, nx);
realNExGrid=real(NExGrid);
imagNExGrid=imag(NExGrid);
realNExGrid(abs(realNExGrid)>1/2)=0.5;
imagNExGrid(abs(imagNExGrid)>2)=2;

NApGrid=reshape(NApZ,nx,nx);
realNApGrid=real(NApGrid);
imagNApGrid=imag(NApGrid);
realNApGrid(abs(realNApGrid)>1/2)=0.5;
imagNApGrid(abs(imagNApGrid)>2)=2;

for i=1:nx
    for j=1:nx
        if any((xPlot(i)-r_center_Coords_Data(1,:)).^2 + (yPlot(j)-r_center_Coords_Data(2,:)).^2 <= r_end_Data.^2+2*eps)
            realNExGrid(j,i) = inf;
            imagNExGrid(j,i) = inf;
            realNApGrid(j,i) = inf;
            imagNApGrid(j,i) = inf;
        end
    end
end
figure(1)
subplot(2,2,1); 
surf(X,Y,realNExGrid);title('Exact N plot, real');view(2); shading interp; colorbar;
subplot(2,2,2); 
surf(X,Y,imagNExGrid);title('Exact N plot, imag');view(2); shading interp; colorbar;
subplot(2,2,3); 
surf(X,Y,realNApGrid);title('Approx N plot, real');view(2); shading interp; colorbar;
subplot(2,2,4); 
surf(X,Y,imagNApGrid);title('Approx N plot, imag');view(2); shading interp; colorbar;



%% End
disp('end')