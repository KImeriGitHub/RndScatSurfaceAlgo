clear

k=30;
NFourier = 2^6;

%delete(gcp('nocreate'))
%parpool(24)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *  S                            |   *                         R  *
% *  S      *                     |                  *          R  *
% *  S                  *         |                             R  *
% *  S                            |          *                  R  *
% *  S           *                |                             R  *
%-*--S----------------------------S-----------------------------R--*
% *  S                        *   |      *                      R  *
% *  S                            |                             R  *     
% *  S   *                        |                    *        R  * 
% *  S                            |*                            R  *

% NDisks=50;
% nd_loop = 1;
% r_center_Coords_Data = zeros(2,NDisks);
% while nd_loop <= NDisks
%     newCoord = (rand(2,1)-0.5).*2; % Vector in [-1,1]^2
%     if ~any(sum((r_center_Coords_Data(:,1:nd_loop-1)-newCoord).^2)<0.1^2) % if newCoord are close to a existing disk or the origin
%         r_center_Coords_Data(:,nd_loop) = newCoord;
%         nd_loop=nd_loop+1;
%     end
% end
% 
% 
% r_init = 0.001;
% r_end_vals = linspace(0.02,0.04,15);
% r_end_Data = r_end_vals(randi(10,1,size(r_center_Coords_Data,2)));
% r_amount = 10;
% 
% % r_center_Coords_Data = [r_center_Coords_Data, [-1.115*ones(1,20);linspace(-1,1,20)],[1.115*ones(1,20);linspace(-1,1,20)]];
% % r_end_Data = [r_end_Data, 0.025*ones(1,40)];
% 
% r_center_Coords_Data = [r_center_Coords_Data, [-1.115*ones(1,10);linspace(-1,1,10)],[1.115*ones(1,10);linspace(-1,1,10)]];
% r_end_Data = [r_end_Data, 0.025*ones(1,20)];
% 
% ztest=[-1.075;  0];
% xtest=[ 1.075;  0];
% 
% xPlot=linspace(-1.3,1.3,100);
% nx=length(xPlot);
% yPlot=linspace(-1.3,1.3,nx);
%
%%% Get Random matrix
% 
% NReceivers = 30;
% 
% nr_loop = 1;
% z = [-1.075 * ones(1,NReceivers);  -100 * ones(1,NReceivers)];
% x = [ 1.075 * ones(1,NReceivers);   100 * ones(1,NReceivers)];
% 
% while nr_loop<=NReceivers
%     zrnd = (rand(1,1)-0.5)*2;
%     xrnd = (rand(1,1)-0.5)*2;
%     
%     if ~any((z(2,:)-zrnd).^2<0.005^2) && ~any((x(2,:)-xrnd).^2<0.005^2)
%         z(2,nr_loop) = zrnd;
%         x(2,nr_loop) = xrnd;
%         nr_loop = nr_loop + 1;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    S                          |   *                         R
%    S     *                    |                  *         R
%   S                  *        |                             R
% S                             |          *                R
%  S           *                |                             R
%------------------------------------------------------------R--
%   S                        *  |      *                      R
%  S                            |                           R       
%    S   *                      |                    *        R
%  S                            |*                           R

NDisks=200;
Disk2DiskThresh=0.03;

r_init = 0.001;
r_end_max = 0.03;
r_end_min = 0.02;
r_amount = 15;

nd_loop = 1;
r_center_Coords_Data = zeros(2,NDisks);
while nd_loop <= NDisks
    newCoord = (rand(2,1)-0.5).*2; % Vector in [-1.0,1.0]^2
    if ~any(sum((r_center_Coords_Data(:,1:nd_loop-1)-newCoord).^2)<(r_end_max*2+Disk2DiskThresh)^2) % if newCoord are close to an existing disk or the origin
        r_center_Coords_Data(:,nd_loop) = newCoord;
        nd_loop=nd_loop+1;
    end
end

r_end_vals = linspace(r_end_min,r_end_max,r_amount);
r_end_Data = r_end_vals(randi(r_amount,1,size(r_center_Coords_Data,2)));

ztest=[-1.1;  0];
xtest=[ 1.1;  0];

xPlot=linspace(-1.15,1.15,500);
nx=length(xPlot);
yPlot=linspace(-1.15,1.15,nx);

% Get Random matrix
NReceivers = 200;
Pt2PtThresh = 0.02;

nr_loop = 1;
z = zeros(2,NReceivers);
while nr_loop<=NReceivers
    zrnd2 = (rand(1,1)-0.5)*2;
    zrnd1 = (rand(1,1)*0.8-1);
    if ~any(sum((r_center_Coords_Data-[zrnd1;zrnd2]).^2)<(r_end_max+Disk2DiskThresh/2)^2) ...
      && ~any(sum((z(:,1:nr_loop-1)-[zrnd1;zrnd2]).^2)<Pt2PtThresh^2)
        z(:,nr_loop) = [zrnd1;zrnd2];
        nr_loop = nr_loop + 1;
    end
end

nr_loop = 1;
x = zeros(2,NReceivers);
while nr_loop<=NReceivers
    xrnd2 = (rand(1,1)-0.5)*2;
    xrnd1 = (rand(1,1)*0.8+0.2);
    if  ~any(sum((r_center_Coords_Data-[xrnd1;xrnd2]).^2)<(r_end_max+Disk2DiskThresh/2)^2) ...
      && ~any(sum((x(:,1:nr_loop-1)-[xrnd1;xrnd2]).^2)<Pt2PtThresh^2)
        x(:,nr_loop) = [xrnd1;xrnd2];
        nr_loop = nr_loop + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
[N_Data, ~, ~, ~, ~, ~, ~, ~, ~] = getEvalCircNeum(ztest, xtest, BigN, r_center_Coords, r_end, NFourier,k);
NApTime_Data = toc;
disp(['Elapsed time = ', num2str(NApTime_Data)]);

%colormap('hot'); imagesc(real(BigN)); colorbar;

%% Plots
%delete(gcp('nocreate'))
%parpool(24)

[X,Y]=meshgrid(xPlot,yPlot);
Z=[reshape(X,1,nx^2);reshape(Y,1,nx^2)];

% RAM saving
NumSplits = max([floor((nx/100)^2*(NFourier/64)^2*(NCenter/30)^2),1]);
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
parfor i=1:NumSplits
    disp(['i = ',num2str(i),' of ',num2str(NumSplits)]);
    Zi=[ZxSplits{i}; ZySplits{i}];
    [a, ~, ~, ~, ~, ~, ~, ~, ~] = getEvalCircNeum(ztest, Zi, BigN, r_center_Coords, r_end, NFourier,k);
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
surf(X,Y,realNApGrid);title('InflationMethod N plot, real');view(2); shading interp; colorbar;
subplot(1,2,2); 
surf(X,Y,imagNApGrid);title('InflationMethod N plot, imag');view(2); shading interp; colorbar;


[NRM, ~, ~, ~, ~, ~, ~, ~, ~,~] = getEvalCircNeum(z, x, BigN, r_center_Coords, r_end, NFourier,k);
UnbiasedNRM = (NRM-sum(NRM,'all')/(size(NRM,1)*size(NRM,2)));
RndMat = UnbiasedNRM/sqrt(var(UnbiasedNRM,0,'all'));

figure(2)
v = svd(RndMat)/sqrt(size(RndMat,1));
histogram(v,floor(length(v)/5),'Normalization','probability')

figure(3)
colormap('hot'); imagesc(real(RndMat)); colorbar;

figure(4)
colormap('hot'); imagesc(imag(BigN)); colorbar;

figure(5)
v = svd(real(RndMat))/sqrt(size(real(RndMat),1));
histogram(v,7, 'Normalization','probability')

figure(8)
v = svd(imag(RndMat))/sqrt(size(real(RndMat),1));
histogram(v,7, 'Normalization','probability')

figure(6)
histogram(real(RndMat(:)), 'Normalization','probability')

figure(7)
histogram(imag(RndMat(:)), 'Normalization','probability')

% figure(8)
% v = svd(BigN(1:1000,1:1000))/sqrt(size(BigN(1:1000,1:1000),1));
% histogram(v(1:20),8,'Normalization','probability')

%% End
disp('end')


%Notes
%What we want:
% t = 1 ;
% r = 1 ;
% n = 50 ;
% m  =   n ;
% v = [ ] ;
% dx = .05 ;
% a   =   0 ;   
% b   =   2 ;
% X = (randn(n)+1i*randn(n))/sqrt(2);
%
% figure(1)
% v = svd(X)/sqrt(size(X,1));
% histogram(v, 'Normalization','probability')
% 
% figure(2)
% colormap('hot'); imagesc(real(X)); colorbar;