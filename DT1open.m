fclose all;
clear all
fn='LINE03_bg+dw+SEC2.dt1';	%file name
i=1;
fid=fopen(fn);
while feof(fid)~=1
    header1 = fread(fid,32,'float');
    headers(:,i)=header1;
    if i==1
        time=(1:1:header1(3))'.*0.8-52.427;
    end
    trace(:,i)=fread(fid,header1(3),'int16');
    i=i+1;
end
fclose(fid);
clear headers header1 i fid fn 
