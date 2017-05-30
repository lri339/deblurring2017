function [g,N] = fastBilateral(f, sigmar, W, eps, flag)

% See [1] K. N. Chaudhury, S. Dabhade, ''Fast and provably accurate bilateral filtering'', IEEE Transactions on Image Processing, vol 26, no. 5, pp. 2519-2528, 2016. (arXiv: http://arxiv.org/abs/1603.08109)

% Gaussian Bilateral filter:
% [g,Ng] = GPA(f, sigmar, sigmas, eps, 'Gauss')
% f             : input image 
% sigmar        : width of range Gaussian
% sigmas        : width of spatial Gaussian
% eps           : kernel approximation accuracy
% g             : output image
% Ng            : approximation order
%
% Box bilateral filter:
% [b,Nb] = GPA(f, sigmar, B, eps, 'box')
% f             : input image 
% sigmar        : width of range Gaussian
% B             : width of box kernel
% eps           : kernel approximation accuracy
% g             : output image
% Nb            : approximation order
%
if strcmp(flag,'Gauss')
    L=round(3*W);
	Hs=fspecial('gaussian',2*L+1,W);
elseif  strcmp(flag,'box')
    L=W;
	Hs=fspecial('average',2*L+1);
else
	error('not enough arguments');
end 
T = 128;
N = estN(sigmar,T,eps);
f=padarray(f,[L,L]);
% changed this to add single()
H=(single(f)-T)/sigmar;     
F=exp(-0.5*H.^2);   
G=ones(size(H));
P=zeros(size(H));  
Q=zeros(size(H));   
Fbar=imfilter(F,Hs);   
for n = 1 : N
	Q=Q+G.*Fbar;
    F=H.*F/sqrt(n);
    Fbar=imfilter(F,Hs);
    P=P+G.*Fbar*sqrt(n);
    G=H.*G/sqrt(n);
end
g= T+sigmar*(P(L+1:end-L,L+1:end-L) ...
    ./Q(L+1:end-L,L+1:end-L));      
g(g<0)=0;
g(g>255)=255;
end

function Nest=estN(sigmar,T,eps)
% sigmar        : width of range Gaussian
% T             : dynamic range of image is [0,2T]
% eps           : kernel approximation accuracy 
% Nest          : approximation order
if  sigmar > 70    
    N=10;
elseif sigmar < 5
    N=800;
else 
    lam=(T/sigmar)^2;
    p = log(exp(1)*lam);
    q = -lam - log(eps);
    t = q*exp(-1)/lam;
    W = t - t^2 + 1.5*t^3 - (8/3)*t^4; 
    N = min(max(q/W,10),300);
    if sigmar < 30
        for iter = 1:5  
            N = N - (N*log(N)-p*N-q)/(log(N)+1-p);
        end
    end 
end
Nest = ceil(N);
end