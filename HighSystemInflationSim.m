clear

k=50;
NFourier = 2^6;
%parpool(24)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% xPlot=linspace(-1.5,2.5,300);
% nx=length(xPlot);
% yPlot=linspace(-2.5,2.5,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%                               |                             S
%                               |       
%                               |    
%                               |      
%                               |     
%-----------------------------------------------------------------------
%       *    *   *    *    *    *    *    *    *    *    *    *
%          *   *   *    *     * |  *    *    *    *    *    *
%       *    *   *    *    *    *    *    *    *    *    *    *
%          *   *   *    *     * |  *    *    *    *    *    *

r_center_Coords_Data   = [-4:0.2:4,        -3.9:0.2:3.9,    -4:0.2:4,        -3.9:0.2:3.9,    -4:0.2:4,        -3.9:0.2:3.9,    -4:0.2:4,        -3.9:0.2:3.9; 
                          -0.1*ones(1,41), -0.2*ones(1,40), -0.3*ones(1,41), -0.4*ones(1,40), -0.5*ones(1,41), -0.6*ones(1,40), -0.7*ones(1,41), -0.8*ones(1,40)];                        
z=[3;  1];
x=[-3; 1];
% 
xPlot=linspace(-4.5,4.5,300);
nx=length(xPlot);
yPlot=linspace(-1,1.5,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                               |                             
%                               |       
%                             * | *  
%                          *    *    * 
%                        *   *  |  *   *
%-----------------------*---*---S---*---*--------------------------------
%                        *   *  |  *   *
%                          *    *    *
%                             * | * 
%                               |
% 
% theta = linspace(-pi,pi,41);
% theta1 = theta(1:2:end-1);
% theta2 = theta(2:2:end-1);
% r_center_Coords_Data   = [0.5*cos(theta1), 0.6*cos(theta2), 0.7*cos(theta1), 0.8*cos(theta2), 0.9*cos(theta1), 1*cos(theta2); 
%                           0.5*sin(theta1), 0.6*sin(theta2), 0.7*sin(theta1), 0.8*sin(theta2), 0.9*sin(theta1), 1*sin(theta2)];                        
% z=[0;  0];
% x=[-1; 0];
% % 
% xPlot=linspace(-2.5,2.5,300);
% nx=length(xPlot);
% yPlot=linspace(-2.5,2.5,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_init = 0.001;
r_end_Data = 0.03*ones(1,size(r_center_Coords_Data,2));
r_amount = 20;
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

disp('Getting BigN');
tic;
%try
    BigN = getCircNeumann(r_center_Coords, r_init, r_end, r_amount, NFourier, k);
%catch
%    warning('Using Inflation method did not work')
%end
[N_Data, ~, ~, ~, ~, ~, ~, ~, ~] = getEvalCircNeum(z, x, BigN, r_center_Coords, r_end, NFourier,k);
NApTime_Data = toc;
disp(['RunTime to get Big N and 1 Eval = ',num2str(NApTime_Data)]);


%% Plots
[X,Y]=meshgrid(xPlot,yPlot);
Z=[reshape(X,1,nx^2);reshape(Y,1,nx^2)];

% RAM saving
NumSplits = max([floor((nx/100)^2*(NFourier/128)^2*(NCenter/50)^2),1]);
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

disp('Getting GridEvaluation');
for i=1:NumSplits
    disp(['i = ',num2str(i),' of ',num2str(NumSplits)]);
    Zi=[ZxSplits{i}; ZySplits{i}];
    [a, ~, ~, ~, ~, ~, ~, ~, ~] = getEvalCircNeum(z, Zi, BigN, r_center_Coords, r_end, NFourier,k);
    NApZSplits(i) = {a};
end
for i=1:NumSplits
    NApZ(Splits(i):Splits(i+1)-1,1)= NApZSplits{i};
end

NApGrid=reshape(NApZ,nx,nx);
realNApGrid=real(NApGrid);
imagNApGrid=imag(NApGrid);
realNApGrid(abs(realNApGrid)>1/2)=0.5;
imagNApGrid(abs(imagNApGrid)>2)=2;

for i=1:nx
    for j=1:nx
        if any((xPlot(i)-r_center_Coords_Data(1,:)).^2 + (yPlot(j)-r_center_Coords_Data(2,:)).^2 <= r_end_Data.^2+2*eps)
            realNApGrid(j,i) = inf;
            imagNApGrid(j,i) = inf;
        end
    end
end
figure(1)
subplot(1,2,1); 
surf(X,Y,realNApGrid);title('InflationMethod N plot, real');view(2); shading interp; colorbar;axis([xPlot(1), xPlot(end), yPlot(1), yPlot(end)]);
subplot(1,2,2); 
surf(X,Y,imagNApGrid);title('InflationMethod N plot, imag');view(2); shading interp; colorbar;axis([xPlot(1), xPlot(end), yPlot(1), yPlot(end)]);



%% End
disp('end')