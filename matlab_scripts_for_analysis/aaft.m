function [zV] = aaft(xV, nsur)
% [zV] = aaft(xV, nsur)
% AAFT (Amplitude Adjusted Fourier Transform )
% For more details see: D. Kugiumtzis, “Surrogate data test for nonlinearity including monotonic transformations,” Phys. Rev. E, vol. 62, no. 1, 2000.
% generates surrogate data for a given time series 'xV' using the AAFT  
% 
% INPUT:
% - xV        : The original time series (column vector)
% - nsur      : The number of surrogate data to be generated 
% OUTPUT:
% - wM        : The CAAFT surrogate data (matrix of 'nsur' columns)
% 

n = length(xV);

% The following gives the rank order, ixV
[oxV,T1] = sort(xV);
[T,ixV]=sort(T1);

% ===== AAFT algorithm 
zV=zeros(n,nsur);
for count=1:nsur;
   rV = randn(n,1); % Rank ordering white noise with respect to xV 
   [orV,T]=sort(rV);   
   yV = orV(ixV); % Y
% >>>>> Phase randomisation (Fourier Transform): yV -> yftV 
   if rem(n,2) == 0
      n2 = n/2;
   else
      n2 = (n-1)/2;
   end
   tmpV = fft(yV,2*n2);
   magnV = abs(tmpV);
   fiV = angle(tmpV);
   rfiV = rand(n2-1,1) * 2 * pi;
   nfiV = [0; rfiV; fiV(n2+1); -flipud(rfiV)];
% New Fourier transformed data with only the phase changed
   tmpV = [magnV(1:n2+1)' flipud(magnV(2:n2))']';
   tmpV = tmpV .* exp(nfiV .* i); 
% Transform back to time domain
   yftV=real(ifft(tmpV,n)); % 3-step AAFT;
% <<<<<
   [T,T2] = sort(yftV); % Rank ordering xV with respect to yftV 
   [T,iyftV] = sort(T2);
   zV(:,count) = oxV(iyftV);  % zV is the AAFT surrogate of xV
end;

return;
