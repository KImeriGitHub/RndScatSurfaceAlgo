function [CosConst, SinConst] = getFourierConst(ft)
%PRE: ft is a matrix of size(M x N). It is the evaluation of a function .
%   f: [-pi,pi)^M-> C^M- at the points linspace(-pi,pi,N+1) (without end)
%   N must be even. N should have low prime divisors (opt: N=2^m)
%   f should be 2pi periodic.
%POST: CosConst and SinConst are matrices of size M x N/2
%Desc: 
%   If f(t)=sum_{n=0}^(K) (a_n cos(nt)+ b_n sin(nt))
%   Then this function returns a_n and b_n as long as K<=N/2
%   Function can be vectorized along 1 dim.
%NOTE:
% - ft can be complex. Then this program divides real part and imag part
%   and combines it again

M=size(ft,1);
N=size(ft,2);

rFfft=fft(real(ft),[],2);
iFfft=fft(imag(ft),[],2);
cjrFfft=conj(rFfft);
cjiFfft=conj(iFfft);

signmat=repmat((-1).^(0:N-1),M,1);
ar=(rFfft+cjrFfft).*signmat/(N*2);
br=(rFfft-cjrFfft).*signmat/(N*2/1i);

ai=(iFfft+cjiFfft).*signmat/(N*2);
bi=(iFfft-cjiFfft).*signmat/(N*2/1i);

%Test
%exact=f(t);
%formu=a*cos(t.*(0:N-1).')+b*sin(t.*(0:N-1).');
%
%All the thing above work but a and b are wrong, and I do not know why. The
%correct thing is the following below:

CosConst=[ar(:,1),ar(:,2:N/2)+fliplr(ar(:,N/2+2:end))]+1i*[ai(:,1),ai(:,2:N/2)+fliplr(ai(:,N/2+2:end))]; %Cutting a one short!
SinConst=[br(:,1),br(:,2:N/2)-fliplr(br(:,N/2+2:end))]+1i*[bi(:,1),bi(:,2:N/2)-fliplr(bi(:,N/2+2:end))]; %Cutting b one short!



end

