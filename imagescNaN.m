function [h,hcb] = imagescNaN(x,y,a,cm,nanclr,colorlim)
% IMAGESC with NaNs assigning a specific color to NaNs

%# find minimum and maximum
if nargin < 6
    amin=min(a(:));
    amax=max(a(:));
else
    amin = colorlim(1);
    amax = colorlim(2);
end
%# size of colormap
n = size(cm,1);
%# color step
dmap=(amax-amin)/n;

%# standard imagesc
him = imagesc(x,y,a);
%# add nan color to colormap
colormap([nanclr; cm]);
%# changing color limits
if amax > amin
    caxis([amin-dmap amax]);
end
%# place a colorbar
hcb = colorbar;
%# change Y limit for colorbar to avoid showing NaN color
if amax > amin
    ylim(hcb,[amin amax])
end

if nargout > 0
    h = him;
end
