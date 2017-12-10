function [ histArray  ] = localHistograms ( img,windowSize )
% Given an image returns a 3D array of histograms – 
% one histogram per window in image.
% Input:
%   img - a grayscale image in the range [0..255]
%   windowSize – a 1x2 array [r c]  defining the num rows and num cols 
%                of image windows.      
% Output:
%   histArray – an NxMx256 array.
%
% For every (i,j) the 1D array histArray(i,j,:) is the histogram
% of an image window of size windowSize whose top left corner 
% pixel is (i,j). Histograms of windows that exceed the boundary 
% of the of img are not included in histArray (thus N = number of 
% rows of img less the number of rows of window +1. 
% Similarly  M=size(img,2)-windowSize(2)+1 ).
%
% Method:   Scans the img pixel by pixel. For each scanned pixel,
% determines the histogram of the image window starting at the 
% pixel and extending  windowSize(1) rows and windowSize(2).

N = size(img,1) - windowSize(1) + 1;
M = size(img,2) - windowSize(2) + 1;
histArray = zeros(N ,M,256);

for i = 1:N
 for j = 1:M
   histArray(i,j,:) = histImage(img(i:i+windowSize(1)-1,j:j+windowSize(2)-1));
 end
end

end