function [isi ta]=freq(t,v)
%Calcuclate frequency variation with time

v=v-mean(v);
crs=[];
isi=[];
j=1;
for i = 1:length(t)-1
    if v(i)<0 && v(i+1)>0
        crs(j)=t(i);
        j=j+1;
    end
end

i=0;
for i = 1:length(crs)-1
    isi(i)=crs(i+1)-crs(i);
end

ta=crs(1:length(crs)-1);



