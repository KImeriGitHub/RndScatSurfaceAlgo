clear
format compact

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

%% Reducing
labels = 0:2;
colYidx = Y(:,1)==1 | Y(:,2)==1 | Y(:,3)==1;
U = U(colYidx,:);
Y = Y(colYidx,1:3);
colYtestidx = Ytest(:,1)==1 | Ytest(:,2)==1 | Ytest(:,3)==1;
Util = Util(colYtestidx,:);
Ytest = Ytest(colYtestidx,1:3);

%randomize test set
randidx = randperm(size(Ytest,1));
Util = Util(randidx,:);
Ytest = Ytest(randidx,:);

% reducing U
% redidx = [];
% for i = 0:2:26
%     redidx = [redidx,(1:2:27)+i*28]; %#ok<AGROW>
% end
% U = U(:,redidx);
% Util = Util(:,redidx);
%% Test

z = [-1.1*ones(1,size(U,2)); linspace(-1.0,1.0,size(U,2))];
%z = [kron(ones(1,28), linspace(-1.35,-1.05,28)); kron(linspace(0.5,-0.5,28), ones(1,28))];
x = [+1.1*ones(1,size(Y,2)); linspace(-0.7,0.7,size(Y,2))];

kRange = linspace(34.1,34.3,24*5);
NRange = [100:10:200, 700,1000];
absDy = @(y)sqrt((y(1,:).'-y(1,:)).^2+(y(2,:).'-y(2,:)).^2);
absDx = @(y)sqrt((y(1,:).'-x(1,:)).^2+(y(2,:).'-x(2,:)).^2);
absDz = @(y)sqrt((y(1,:).'-z(1,:)).^2+(y(2,:).'-z(2,:)).^2);

parpool(24);
fData=zeros(length(NRange),length(kRange));
for i=1:length(NRange)
    disp(['i = ', num2str(i)]);
    N = NRange(i);
    parfor j = 1:length(kRange)
        disp(['j = ', num2str(j)]);
        
        k = kRange(j);
        accuracy = zeros(1,length(labels)); %preload
        
        [y, ~] = optimalDiskLoc(U, Y, z, x, k, N);
        
        M = -1i/4*besselh(0,k*absDy(y))/N;
        M(1:1+size(M,1):end) = 1;
        Nk = -(-1i/4*besselh(0,k*absDx(y))).'*(M\(-1i/4*besselh(0,k*absDz(y))))/N;
        
        Ytil = abs(Util*Nk.');
        [~,estIdxLabel] = max((Ytil).');
        estLabel = (estIdxLabel-1).';
        for q=labels
            qIdx = (Ytest(:,q+1) == 1);
            accuracy(q+1) = sum(estLabel(qIdx) == q)/sum(qIdx);
        end
        fData(i,j) = mean(accuracy);
    end
end

%% Test Outcome 
k=34.1;
N=700;
disp(['k : ', num2str(k)]);
%y=rand(2,N);
[y, nor] = optimalDiskLoc(U, Y, z, x, k, N);

Gamma = @(y)-1i/4*besselh(0,k*absDy(y));
Gammax = @(y)-1i/4*besselh(0,k*absDx(y));
Gammaz = @(y)-1i/4*besselh(0,k*absDz(y));
M = Gamma(y)/N;
M(1:1+size(M,1):end) = 1;

Nk = -Gammax(y).'*(M\Gammaz(y))/N;
norUY = sum(abs(U*Nk.'-Y), 'all'); 
disp(['nor = ', num2str(norUY)]);
plot(y(1,:),y(2,:),'.g','MarkerSize',5); hold on;
plot(z(1,:),z(2,:),'.r','MarkerSize',5); 
plot(x(1,:),x(2,:),'.b','MarkerSize',5); 
axis([-1.2, 1.2, -1.2, 1.2]);hold off;

% Applying Test Set

accuracy = zeros(1,length(labels)); %preload

Ytil = abs(Util*Nk.');
sigCostTemp = 1./(1+exp(-abs(Ytil)));
[~,estIdxLabel] = max((Ytil).');
estLabel = (estIdxLabel-1).';
for q=labels
    qIdx = (Ytest(:,q+1) == 1);
    accuracy(q+1) = sum(estLabel(qIdx) == q)/sum(qIdx);
    disp(['Accuracy of label ',num2str(q), ' is ', num2str(accuracy(q+1))]);
end
disp(['Mean Accuracy : ', num2str(mean(accuracy))]);
disp('end');
