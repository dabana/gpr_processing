function [kdemat,xmesh] = localkde(trace,windowSize,nmesh,MIN,MAX)
% Given an 2D array returns a 3D array of probability densities using kernel 
% density estimation for every window in the array

% Input:
%   trace - a 2D array
%   windowSize – a 1x2 array [r c]  defining the num rows and num cols 
%                of image windows.
%   nmesh - the number of mesh points used for the kernel density
%   estimation
%   MIN - minimum value for kde
%   MAX - maximum value for kde
% Output:
%   kdemat – an NxM cellarray.
%
% For every (i,j) the 1D array histArray(i,j,:) is the histogram
% of an image window of size windowSize whose top left corner 
% pixel is (i,j). Histograms of windows that exceed the boundary 
% of the of img are not included in histArray (thus N = number of 
% rows of img less the number of rows of window +1. 
% Similarly  M=size(img,2)-windowSize(2)+1 ).
%


N = size(trace,1) - windowSize(1) + 1;
M = size(trace,2) - windowSize(2) + 1;
kdemat = cell(N ,M);

for i = 1:N
 for j = 1:M
   [bandwidth,kdemat{i,j},xmesh,cdf] = kde(reshape(trace(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1),nmesh,MIN,MAX);
 end
end

end