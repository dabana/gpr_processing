function [mumat,sigmamat] = localnorm(stp,windowSize)
% Given an 2D array returns two 2D arrays containing the parameters (mu and
% sigma)
% of the normal distribution fitting the data

% Input:
%   trace - a 2D array
%   windowSize – a 1x2 array [r c]  defining the num rows and num cols 
%                of image windows.
%   nmesh - the number of bins
% Output:
%   sigmamat – an NxM array.
%   mumat – an NxM array.

% For every (i,j) the 1D array histArray(i,j,:) is the histogram
% of an image window of size windowSize whose top left corner 
% pixel is (i,j). Histograms of windows that exceed the boundary 
% of the of img are not included in histArray (thus N = number of 
% rows of img less the number of rows of window +1. 
% Similarly  M=size(img,2)-windowSize(2)+1 ).
%


N = size(stp,1) - windowSize(1) + 1;
M = size(stp,2) - windowSize(2) + 1;
sigmamat = zeros(N,M);
mumat=zeros(N,M);

for i = 1:N
 for j = 1:M
   [mu,sigma]=normfit(reshape(stp(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1));
   sigmamat(i,j)=sigma;
   mumat(i,j)=mu;
 end
end

end