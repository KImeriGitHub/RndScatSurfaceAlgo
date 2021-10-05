clear
format compact

%% Initializing

%% MNIIST WORKING SET
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest,YTest,anglesTest] = digitTest4DArrayData;

disp(size(XTrain))
disp(size(YTrain))
disp(size(XTest))

Y = zeros(5000,10);
Ytest = zeros(5000,10);
labels = 0:9;
for n=1:5000
    Y(n,labels(YTrain(n))+1) = 1;
    Ytest(n,labels(YTest(n))+1) = 1;
end

U = reshape(XTrain,28*28,5000).';
Util = reshape(XTest,28*28,5000).';

%% AUDIO DIGITS RECODING WORKING SET
% pathToRecordingsFolder = fullfile('C:\Users\KI-GPC\Desktop\CurrentProject\Matlab Code','free-spoken-digit-dataset-master','recordings');
% location = pathToRecordingsFolder;
% ads = audioDatastore(location);
% tic;
% ads.Labels = helpergenLabels(ads);
% %summary(ads.Labels); %Summary
% 
% ads = shuffle(ads);
% [adsTrain,adsTest] = splitEachLabel(ads,0.8);
% countEachLabel(adsTrain)
% countEachLabel(adsTest)
% 
% XTrain = [];
% XTest = [];
% scatds_Train = transform(adsTrain,@(x)helperReadSPData(x));
% scatds_Test = transform(adsTest,@(x)helperReadSPData(x));
% while hasdata(scatds_Train)
%     smat = read(scatds_Train);
%     XTrain = cat(2,XTrain,smat);
% end
% while hasdata(scatds_Test)
%     smat = read(scatds_Test);
%     XTest = cat(2,XTest,smat);
% end
% 
% adsTrain.Labels = helpergenLabels(adsTrain);
% adsTest.Labels = helpergenLabels(adsTest);
% YTrain = double(adsTrain.Labels)-1;
% YTest = double(adsTest.Labels)-1;
% Y = zeros(size(YTrain,1),10);
% labels = 0:9;
% for n=1:size(YTrain,1)
%     Y(n,labels(YTrain(n)+1)+1) = 1;
% end
% 
% U = XTrain.';
% Util = XTest.';
% 
% U = U(:,2000:2:6000);
% Util = Util(:,2000:2:6000);
% 
% toc
% % Reducing Labels to 3
% labels = 0:2;
% redTrainIdx = YTrain<3;
% redTestIdx = YTest<3;
% 
% Ured = U(redTrainIdx,:);
% Utilred = Util(redTestIdx,:);
% 
% U = Ured;
% Y = Y(redTrainIdx,1:3);
% Util = Utilred;
% YTest = YTest(YTest<3);

%% Getting Opt Circles
k=30;

r_init = 0.0005;
r_end = 0.03;
r_amount = 30;
r_inc_Data = (logspace(r_init,r_end,r_amount)-10^r_init)*(r_end-r_init)/(10^r_end-10^r_init)+r_init;
%r_inc_Data = (exp(linspace(0,log(2),r_amount))-1).*(r_end-r_init)+r_init;
%r_inc_Data = (2*r_end+r_init)/2+0.5*(2*r_end-r_init).*cos(pi*(2*(1:(2*r_amount+1))-1)/2/(2*r_amount+1));
%r_inc_Data = [r_init, fliplr(r_inc_Data(r_amount+2:end-1)), r_end];
itermax = 150;
NZeta = 300;
NFourier = 2^6;

%z = [kron(ones(1,28), linspace(-1.15,-1.05,28)); kron(linspace(0.5,-0.5,28), ones(1,28))];
z = [-1.1*ones(1,size(U,2)); linspace(-0.5,0.5,size(U,2))];
x = [+1.1*ones(1,size(Y,2)); linspace(-0.5,0.5,size(Y,2))];
y = optimalDiskLoc(U, Y, z, x, 10, 100);

tic
[NMat, CenterCoord, DiskRad, min_val_Data, BigN] = getOptCircles(U,Y, k, z, x, r_inc_Data, itermax, NZeta, NFourier);
elaTime = toc;

Gz = getGamma(z,[0;0],k);
Gx = getGamma(x,[0;0],k);

temp = Gx*(pinv(Y)*U)*Gz.';


BigNExact = getCircNeumann(CenterCoord, r_init, DiskRad, r_amount, NFourier,k);

disp(['Former Norm:',num2str(min_val_Data(1,1))]);
disp(['New Norm:',num2str(min_val_Data(1,end))]);
disp(['Elapsed time: ',num2str(elaTime)]);
%Make NMatExact
sepPtr=1; D=[];
tPts_ip = linspace(-pi,pi,NFourier+1);
tPts_ip = tPts_ip(1:end-1);
for i=1:size(CenterCoord,2)
    D=[D, DiskRad(i)*[cos(tPts_ip);sin(tPts_ip)]+CenterCoord(:,i)];
    sepPtr=[sepPtr, size(D,2)+1]; %#ok<*AGROW>
end
Dom = shape.Spline(D, sepPtr, 3);
%NMatExact = getNeumannFct(Dom, x, z, k);

colormap('hot'); imagesc(real(BigNExact-BigN)); colorbar;
figure(2); plot(Dom); axis([-1,1,-1,1]);
%figure(2); colormap('hot'); imagesc(real(NMatExact-NMat)); colorbar;

%% Applying Test Set

accuracy = zeros(1,length(labels)); %preload

Ytil = abs(Util*NMat).^2;
sigCostTemp = 1./(1+exp(-Ytil));
[~,estIdxLabel] = max((Ytil).');
estLabel = (estIdxLabel-1).';
for q=1:length(labels)
%     qIdx = (Ytest(:,q)==1);
%     accuracy(q) = sum(estLabel(qIdx) == q-1)/length(find(qIdx==1));
    qIdx = (YTest==q-1);
    accuracy(q) = sum(estLabel(qIdx) == YTest(qIdx))/length(YTest(qIdx));
    disp(['Accuracy of label ',num2str(labels(q)), ' is ', num2str(accuracy(q))]);
end
disp(['Mean Accuracy : ', num2str(mean(accuracy))]);
disp('end');