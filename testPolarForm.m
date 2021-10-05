clear

sepPtr=[1]; D=[];

D=[D, [0.1*cos(linspace(0,15*pi/8,32))+0.8;0.1*sin(linspace(0,15*pi/8,32))+0.2]];
sepPtr=[sepPtr, size(D,2)+1];

D=[D, [0.1*cos(linspace(0,15*pi/8,32))+0.4;0.1*sin(linspace(0,15*pi/8,32))+0.4]];
sepPtr=[sepPtr, size(D,2)+1];

Dom = shape.Spline(D, sepPtr, 3);

k=2;
NFourier=2^6;
Source = [-1;-1];
discDis = 0.0001;
r=0.005;
tPts_ip = linspace(-pi,pi,NFourier+1);
tPts_ip = tPts_ip(1:end-1);
nor=[cos(tPts_ip);sin(tPts_ip)];
Pts=r*nor;

NF = getNeumannFct(Dom,Source,Pts,k).';
dNF = diff([NF, NF(1)])/(tPts_ip(2)-tPts_ip(1));
ddNF = diff([dNF, dNF(1)])/(tPts_ip(2)-tPts_ip(1));

nabla1NF = getNeumannFct(Dom, Source,Pts+[discDis;0],k).';
nabla1NF = (nabla1NF-NF)/discDis;
nabla2NF = getNeumannFct(Dom, Source,Pts+[0;discDis],k).';
nabla2NF = (nabla2NF-NF)/discDis;

delrNF = (getNeumannFct(Dom,Source,Pts + discDis*nor,k).'-NF)/discDis;
del2rNF = ((getNeumannFct(Dom,Source,Pts + 2*discDis*nor,k).'-getNeumannFct(Dom,Source,Pts + discDis*nor,k).')...
          -(getNeumannFct(Dom,Source,Pts + discDis*nor,k).'-getNeumannFct(Dom,Source,Pts,k).'))/discDis^2;
     

A = (k^2*NF + del2rNF + delrNF/r +(ddNF)/r^2).' ;
B = (ddNF/r^2).';
C = (k^2*NF + del2rNF + delrNF/r).';