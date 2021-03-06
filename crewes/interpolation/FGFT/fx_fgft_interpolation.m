function Y=fx_fgft_interpolation(X,pctg)
% This function interpolates regularly sampled seismic traces using
% a fast generalized Fourier transform. The Spitz method of using the low
% frequency portion for the interpolation of high frequencies is used.
%
%% Input
%   X: Data in t-x domain
%   pctg: thresholding percentage (0.8 means 80 percent forced to be zero)
%
%% Output
%   Y: Interpolated data
%
% Author: Mostafa Naghizadeh; Copyright (C) 2010
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its author (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

[nt,nx]=size(X);

% Find the nearest power of two to padd with zeros
NT=2^(nextpow2(nt));
NX=2^(nextpow2(nx));

TFX1=fft(X,2*NT,1);
FX1=zeros(2*NT,NX);
FX1(:,1:nx)=TFX1;

TFX2=fft(X,NT,1);
FX2=zeros(NT,2*NX);
FX2(:,1:2:2*nx)=TFX2;


% Starting loop over frequencies to be reconstructed
FY=zeros(NT,2*NX);
for ia=1:(NT/2)+1
    %Selecting frequency slices form a frequency and its half
    x1=FX1(ia,:); %used for estimating Mask function
    x2=FX2(ia,:); %used for interpolation using mask

    % Forward fgft of low half frequency
    [D1]=forward_fgft(x1);

    % Making mask function by threshholding D1
    cfs = sort(abs(D1));
    nb = round(pctg*length(cfs));
    cutoff = cfs(nb);
    % Set small coefficients to zero
    for w=1:length(D1)
        if (abs(D1(w))>cutoff)
            MASK1(w)= 1;
        else
            MASK1(w)=0;
        end
    end

    MASK2=matrix_scaling(MASK1,1,2*NX);
    
    % Fitting mask to available data of to be interpolated section
    INTD=ls_mask_fitting_fgft(x2(1:2:2*nx),[1:2:2*nx],MASK2,2*NX,3);
    
    % Inverting D2 using inverse fgft and saving it in final array
    FY(ia,:)=INTD;
end

% Making the symmetric part of the frequency domain
FY(NT/2+2:NT,:)=conj(flipud(FY(2:NT/2,:)));

% Inverse Fourier transform from f-x to t-x domain
INTD=real(ifft(FY,NT,1));
Y=INTD(1:nt,1:2*nx);

% Normalizing trace values
for is=1:2*nx
    Y(:,is)= (Y(:,is) - mean(Y(:,is)))./std(Y(:,is));
end
