function cercle(x,y,r)
%
%   cercle(x,y,r)
%
%   Cette fonction trace un cercle de rayon r et de centre (x,y)
%
% Copyright :  G. GARCIA

clear xx yy zz
xx = linspace((x-r),(x+r),200);
for i=1:length(xx)
      tt(i) = abs(r^2-(xx(i)-x)^2);
      yy(i) = sqrt(tt(i));
      zz(i) = -sqrt(tt(i));
end
plot(xx,zz,'r-.',xx,yy,'r-.')

%xx = tt';
%yy = yy';
%zz = zz';
