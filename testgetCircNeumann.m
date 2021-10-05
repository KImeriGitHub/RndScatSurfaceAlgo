clear

k=30;
NFourier = 2^6;
                                                   
r_center_Coords_Data   = [-2:0.8:2,-1.6:0.8:1.6, -2:0.8:2,-1.6:0.8:1.6; 
                           -0.1*ones(1,6),-0.2*ones(1,5),-0.3*ones(1,6),-0.4*ones(1,5)];                        
z=[2;  1];
x=[-2; 1];
% 
xPlot=linspace(-2.5,2.5,100);
nx=length(xPlot); %NOT RAM EFFICIENT
yPlot=linspace(-1,1.5,nx);     


r_init = 0.001;
r_end_Data = 0.03*ones(1,size(r_center_Coords_Data,2));
r_amount = 20;

NData= size(r_center_Coords_Data,2);

tPts_ip = linspace(-pi,pi,NFourier+1);
tPts_ip = tPts_ip(1:end-1);

%FOR PLOTTING DOM
sepPtr=[1]; D=[];
for i=1:size(r_center_Coords_Data, 2)
    D=[D, r_end_Data(i)*[cos(tPts_ip);sin(tPts_ip)]+r_center_Coords_Data(:,i)];
    sepPtr=[sepPtr, size(D,2)+1];
end
Dom = shape.Spline(D, sepPtr, 3);
figure(1)
plot(Dom);axis equal;
pause(1);

NExact_Data     = zeros(1,NData);
N_Data          = zeros(1,NData);
NExactTime_Data = zeros(1,NData);
NApTime_Data    = zeros(1,NData);
error_Data      = zeros(1,NData);
timeQuot_Data   = zeros(1,NData);

for n=1:NData
    disp(['n = ',num2str(n), ' out of ',num2str(NData)])
    r_center_Coords = r_center_Coords_Data(:,1:n);
    r_end = r_end_Data(:,1:n);
    
    sepPtr=[1]; D=[];
    for i=1:length(r_end)
        D=[D, r_end(i)*[cos(tPts_ip);sin(tPts_ip)]+r_center_Coords(:,i)];
        sepPtr=[sepPtr, size(D,2)+1];
    end
    Dom = shape.Spline(D, sepPtr, 3);
    
    tic;
    NExact_Data(1,n) = getNeumannFct(Dom, z, x, k).';
    NExactTime_Data(1,n) = toc;  
    
    tic;
    try
        BigN = getCircNeumann(r_center_Coords, r_init, r_end, r_amount, NFourier, k);
    catch
        warning('Using Inflation method did not work')
        NExact_Data     = NExact_Data(1,1:n-1);
        N_Data          = N_Data(1,1:n-1);
        NExactTime_Data = NExactTime_Data(1,1:n-1);
        NApTime_Data    = NApTime_Data(1,1:n-1);
        error_Data      = error_Data(1,1:n-1);
        timeQuot_Data   = timeQuot_Data(1,1:n-1);
        NData = n-1;
        r_center_Coords_Data = r_center_Coords_Data(:,1:n-1);
        r_end_Data = r_end_Data(1,1:n-1);
        BigN = getCircNeumann(r_center_Coords_Data, r_init, r_end_Data, r_amount, NFourier, k);
        
        sepPtr=[1]; D=[];
        for i=1:n-1
            D=[D, r_end(1,i)*[cos(tPts_ip);sin(tPts_ip)]+r_center_Coords(:,i)];
            sepPtr=[sepPtr, size(D,2)+1];
        end
        Dom = shape.Spline(D, sepPtr, 3);
        
        break;
    end
    [N_Data(1,n), ~, ~, ~, ~, ~, ~, ~, ~] = getEvalCircNeum(z, x, BigN, r_center_Coords, r_end, NFourier,k);
    NApTime_Data(1,n) = toc;
    
    error_Data(1,n) = abs(N_Data(1,n)-NExact_Data(1,n));
    timeQuot_Data(1,n) = NExactTime_Data(1,n)./NApTime_Data(1,n);
    
    disp(['n time = ',num2str(NApTime_Data(1,n)+NExactTime_Data(1,n))])
end
%colormap('hot'); imagesc(real(BigN)); colorbar;

%% Analysis
figure(2)
subplot(3,1,1)
plot(3:NData, error_Data(3:NData))
title('Error')
xlabel('Number of circles')
ylabel('diff')
subplot(3,1,2)
plot(3:NData, timeQuot_Data(1,3:NData))
title('Time Quotient (Exact/Approx)')
xlabel('Number of circles')
ylabel('Quotient')
subplot(3,1,3)
plot(3:NData, NApTime_Data(1,3:NData))
title('Runtime Inflation Method')
xlabel('Number of circles')
ylabel('Seconds')

%% Plots
[X,Y]=meshgrid(xPlot,yPlot);
Z=[reshape(X,1,nx^2);reshape(Y,1,nx^2)];

NExZ=getNeumannFct(Dom,z,Z,k);
[NApZ, ~, ~, ~, ~, ~, ~, ~, ~] = getEvalCircNeum(z, Z, BigN, r_center_Coords, r_end, NFourier,k);

NExGrid=reshape(NExZ,nx,nx);
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
figure(3)
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
% 
% Zr=reshape(N,100,100);
% 
% figure(1);
% A=diff(Zr,2,1);
% A(:,[1, end])=[];
% B=diff(Zr,2,2);
% B([1, end],:)=[];
% C=(0.1)^2*Zr;
% C([1, end],:)=[];
% C(:,[1, end])=[];
% D=(A+B)/min(abs(diff(x)))^2+C;
% 
% figure(1)
% subplot(2,1,1)
% surf(real(A/min(abs(diff(x)))^2+B/min(abs(diff(y)))^2+C));
% subplot(2,1,2)
% surf(imag((A+B)/min(abs(diff(x)))^2+C));
% 
% figure(2)
% subplot(2,1,1)
% surf(X,Y,real(Zr));
% subplot(2,1,2)
% surf(X,Y,imag(Zr));
