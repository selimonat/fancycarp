% TFALS Truncated Fourier Asymmetric Least Squares method of baseline estimation
% TFALS estimates the baseline of a signal vector xab.
% Returns the estmated baseline xb and basline corrected vector xa %
% Input:
% xab=Raw signal vector [variable,sample];Column vector
% nf=Total number of frequency components needs to accomodate baseline
% p=Asymmetry parameter (0.001>=p<=0.1) %
% Choice of nf depends on the shape/frequency of baseline; linear baseline usually
% requires nf=2, more variable baseline requires higher value of nf. %
% Value of p depends on the noise level of raw signal vector; Reletivly lower p
% value require for signal with high noise or higher value for vice versa %
% Output:
% xb=(Estimated baseline)
% xa=(Baseline corrected data)
%
function [xa,xb]=TFALS(xab,nf,p)
if nargin==2
    p=0.001;
end
%% Normalized Fourier basis sets
N=size(xab,1);
t=0:N-1;
nb=nf*2+1;
z=zeros(nb,N);
z(1,:)=ones(1,N);
fs=0.25;
for i=1:nf
    z(2*i,:)=cos((2*pi*fs*t)/N);
    z(2*i+1,:)=sin((2*pi*fs*t)/N);
    if fs < 0.5
        fs = 0.25+fs;
    elseif fs == 0.5;
        fs = 0.5+fs;
    else
        fs = 1+fs;
    end
    % Calculate number of data points in xab
    % calculate basis functions % calculate basis functions    
end
for m=1:nb;
    z(m,:)=z(m,:)/norm(z(m,:));
end
%% Orthonormal Basis sets 
[~ , ~, P]=svd((z),'econ');
%% Asymmetric least squares 
w=ones(N,1);
target=1; 
while target
    W=spdiags(w,0,N,N); 
    bw=(P'*W); 
    q=(bw*P)\(bw*xab); 
    xb=P*q;
    w0=w;
    w(xab>(xb))=p; w(xab<=(xb))=(1-p);
    target = sum(abs(w - w0)) > 0;
end
xa=xab-xb;
% Normalize basis functions
% Orthogonalize basis functions
% Generate sparse matrix of asymmetric weights
% perform least squares fit % Estimate baseline
% Calculate baseline corrected signal vector