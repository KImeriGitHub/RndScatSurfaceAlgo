k=10;
Source=[-0.1;0];
Receiver=[0.1;0];
Zeta=[0.5;0.5]; 
discDis = 0.000001;

D0=shape.Ellipse(0.1,0.1,1000)+[0;-0.3];
r=0.0001;
Dr=shape.Ellipse(r,r,1000)+Zeta;
Dom=merge(D0,Dr);


x=linspace(-r,r,100)+Zeta(1);
nx=length(x);
y=linspace(r*0.75,r*1.25,nx)+Zeta(2);
[X,Y]=meshgrid(x,y);

Gamma = -1i/4*besselh(0,k*norm(Zeta+[0;r]-Source));
Z=[reshape(X,1,nx^2);reshape(Y,1,nx^2)];
%NS=getNeumannFct(Dom,Source,Z,k)-(-1i/4*besselh(0,k*sqrt((Z(1,:)-Source(1,:)).^2+(Z(2,:)-Source(2,:)).^2))).';
NS=getNeumannFct(Dom,Source,Z,k);
NSp=getNeumannFct(Dom,Source,Z+[0; discDis],k);
dNS = (NSp - NS)/discDis;

%NR=getNeumannFct(Dom,Receiver,Z,k)-(-1i/4*besselh(0,k*sqrt((Z(1,:)-Receiver(1,:)).^2+(Z(2,:)-Receiver(2,:)).^2))).';
NR=getNeumannFct(Dom,Receiver,Z,k);

ZS=reshape(dNS,nx,nx);
%ZR=reshape(NR,nx,nx);
realZS=real(ZS);
imagZS=imag(ZS);
realZS(abs(realZS)>1/2)=0.5;
imagZS(abs(imagZS)>2)=2;

figure(1)
subplot(2,2,1); 
surf(X,Y,realZS);title('Source plot, real');
subplot(2,2,2); 
surf(X,Y,imagZS);title('Source plot, imag');

% subplot(2,2,3); 
% surf(X,Y,real(ZR));title('Receiver plot, real');
% subplot(2,2,4); 
% surf(X,Y,imag(ZR));title('Receiver plot, imag');