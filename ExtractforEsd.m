function [Signal] = ExtractforEsd(Es_d,SignNoiseMeasure,out_n_th,d_clip,smthMat)
%{
Function used to perform (3 part 2): Smoothly zero samples containing
significant noise and a recommanded (but facultative) step intended to reduce artifacts

OUTPUT

Signal - Extracted signal matrix

INPUT

Es_d - E(s|d) the esperance value of the signal s knowing the
samples d and a noise estimate n (use function Es_d.m)

d_clip - the samples d but cliped by WindowSize at their edge. WindowSize
is not an input of the function, so beware.

SignNoiseMeasure - Matrix containing the signal to noise measure used to
create ZeroMask. it can be either SignNoiseMeasure or SignalNoiseMeasure

out_n_th - outter noise threashold used to create ZeroMask using either
SignNoiseMeasure or SignNoiseMeasure2.
%}
%-------------------------------------------------------------------------%

%(3 part1) Zero all samples containing significant percentage of noise

%creating the zeroing mask
ZeroMask=zeros(size(Es_d,1),size(Es_d,2));
for i = 1:size(Es_d,1)
    for j = 1:size(Es_d,2)
        if SignNoiseMeasure(i,j)<out_n_th;
        ZeroMask(i,j)=0;
        else
        ZeroMask(i,j)=1;
        end
    end
end
%smoothing it
h=fspecial('gaussian',[5 5]);
ZeroMaskAnalytic=imfilter(ZeroMask,h);
%applying it
Emask=Es_d.*ZeroMaskAnalytic;

%(Recommanded step) Smooth an array of E(s|d)/d values both spatially and temporally before
%multiplying by d vector

    Eond=Emask./d_clip;
    [Inan,Jnan] = find(isnan(Eond));%kill NaNs
    [Iinf,Jinf] = find(isinf(Eond));%kill INFs
    [Ibig,Jbig] = find(abs(Eond)>100);%kill big values
    Iout=[Inan;Iinf;Ibig];
    Jout=[Jnan;Jinf;Jbig];
    for i=1:size(Iout,1)
        for j=1:size(Jout,1)
            Eond(Iout(i),Jout(j))=0;
        end
    end
    Signal=smooth2a(Eond,smthMat(1),smthMat(2)).*d_clip;


end

