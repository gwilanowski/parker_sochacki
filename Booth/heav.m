function h=heav(t)
% Heaviside function (Kurian et al. 2010)    
% Define the heaviside function, since the command `heaviside(t)' 
% in matlab gives NaN when t = 0.

        h=zeros(size(t));
        h(t>0)=1;
end  