function [pdfmat,xmesh] = localpdfs(trace,windowSize,nmesh)
% Given an 2D array returns a 3D array of probability densities using 
% histogram normalisation

% Input:
%   trace - a 2D array
%   windowSize – a 1x2 array [r c]  defining the num rows and num cols 
%                of image windows.
%   nmesh - the number of bins
% Output:
%   pdfmat – an NxM cellarray.
%
% For every (i,j) the 1D array histArray(i,j,:) is the histogram
% of an image window of size windowSize whose top left corner 
% pixel is (i,j). Histograms of windows that exceed the boundary 
% of the of img are not included in histArray (thus N = number of 
% rows of img less the number of rows of window +1. 
% Similarly  M=size(img,2)-windowSize(2)+1 ).


N = size(trace,1) - windowSize(1) + 1;
M = size(trace,2) - windowSize(2) + 1;
pdfmat = cell(N ,M);

for i = 1:N
 for j = 1:M
   [histo,xmesh]=hist(reshape(trace(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1),nmesh);
   pdfmat{i,j}=histo/(windowSize(1)*windowSize(2))/diff(xmesh(1:2));
 end
end

end