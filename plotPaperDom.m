tPts = linspace(-pi,pi,1000+1);
tPts = tPts(1:end-1);

N = 20;
r_end = rand(1,N)/4+0.1;
ksep=randi([6,10],N);
nd_loop = 1;
r_center_Coords = zeros(2,N);
while nd_loop <= N
    newCoord = (rand(2,1)-0.5).*2;
    if nd_loop==1 || ~any(sum((r_center_Coords(:,1:nd_loop-1)-newCoord).^2,1)<0.5*(r_end(1:nd_loop-1)+r_end(nd_loop)).^2)
        r_center_Coords(:,nd_loop) = newCoord;
        nd_loop=nd_loop+1;
    end
end

sepPtr=[1]; D=[];
for i=1:length(r_end)
    D=[D, r_end(i)*((cos(ksep(i)*tPts+rand)*r_end(i)+1.2)/3.5.*[cos(tPts);sin(tPts)])+r_center_Coords(:,i)];
    sepPtr=[sepPtr, size(D,2)+1];
end
Dom = shape.Spline(D, sepPtr, 100);

plot(Dom,'-k'); hold on;

z = [-1.2*ones(1,30); linspace(-1.0,1.0,30)];
%z = [kron(ones(1,28), linspace(-1.35,-1.05,28)); kron(linspace(0.5,-0.5,28), ones(1,28))];
x = [+1.2*ones(1,15); linspace(-0.7,0.7,15)];

plot(z(1,:),z(2,:),'.r','MarkerSize',11); 
plot(x(1,:),x(2,:),'.b','MarkerSize',11); 
axis([-1.3, 1.3, -1.3, 1.3]);hold off;