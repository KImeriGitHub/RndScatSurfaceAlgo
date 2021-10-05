f=@(t) [1+3i*cos(2*t)+1*cos(3*t);...
        -0.5+-0.2i*sin(2*t)+8*cos(5*t)];

N=2^5; % Must be even! Preferably has low prime divisors.
assert(mod(N,2)==0);

t= linspace(-pi,pi,N+1);
t=t(1:end-1);

F=f(t);
[CosConst, SinConst] = getFourierConst(F);

hold on;
plot(t,imag(CosConst*cos((0:N/2-1).'*t)+SinConst*sin((0:N/2-1).'*t)),'b')
plot(t,imag(F),'r')
