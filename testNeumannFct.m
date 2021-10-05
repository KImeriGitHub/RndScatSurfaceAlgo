clear

sepPtr=[1]; D=[];

D=[D, [0.1*cos(linspace(0,15*pi/8,32))+0.8;0.1*sin(linspace(0,15*pi/8,32))+0.2]];
sepPtr=[sepPtr, size(D,2)+1];

D=[D, [0.1*cos(linspace(0,15*pi/8,32))+0.4;0.1*sin(linspace(0,15*pi/8,32))+0.4]];
sepPtr=[sepPtr, size(D,2)+1];

Dom = shape.Spline(D, sepPtr, 3);

x=linspace(0.4,0.8,100);
y=linspace(0.1,0.5,100);
assert(length(x)==length(y) && abs(x(end)-x(1))==abs(y(end)-y(1)));
[X,Y]=meshgrid(x,y);

Z=[reshape(X,1,100^2);reshape(Y,1,100^2)];
N=getNeumannFct(Dom,[0.6; 0.2],Z,0.1);

Zr=reshape(N,100,100);

figure(1);
A=diff(Zr,2,1);
A(:,[1, end])=[];
B=diff(Zr,2,2);
B([1, end],:)=[];
C=(0.1)^2*Zr;
C([1, end],:)=[];
C(:,[1, end])=[];
D=(A+B)/min(abs(diff(x)))^2+C;

figure(1)
subplot(2,1,1)
surf(real(A/min(abs(diff(x)))^2+B/min(abs(diff(y)))^2+C));
subplot(2,1,2)
surf(imag((A+B)/min(abs(diff(x)))^2+C));

figure(2)
subplot(2,1,1)
surf(X,Y,real(Zr));
subplot(2,1,2)
surf(X,Y,imag(Zr));
