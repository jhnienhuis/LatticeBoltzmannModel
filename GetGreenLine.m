function [x,t,z,fitresult] = GetGreenLine(img_list)
%GREENLINEEXTRACTOR ( folder ) retrieves and corrects the green line from
%our SLR flume photos. 

% Input: [img list]
% A numerical list of IMG_xxxx.JPG photos.
%
% Output: [x,y,z,fitresult]
% x,y,z is along flume vector, time vector and height matrix in meters,
% including NaN's for all points where the green line hasn't been
% succesfully extracted
% fitresult is a sfit matlab structure, where the x,y and z have been
% fitted to a surface cubic interpolator (for the missing values).
% 
% Method: 
% very simple imread, distortion correction, angle correction and
% perspective correction (only in x) calibrated with a picture of a grid,
% placed in the flume.
%
%

% if the matlabpool isnt open, open it.
%if ~matlabpool('size'),
%    matlabpool open
%end

% configure my parfor progress monitor

filepath = [dropbox filesep 'work' filesep '_Tools' filesep 'wavewatch' filesep 'ParforProgMon'];

addpath(filepath)
pctRunOnAll javaaddpath C:\Users\jaap\Dropbox\work\_Tools\wavewatch\ParforProgMon\


% read first image, retrieve information/timestamp
img = imread(['IMG_' num2str(img_list(1),'%04.0f') '.JPG']);
exif =imfinfo(['IMG_' num2str(img_list(1),'%04.0f') '.JPG']);

% get image size
[M N ~] = size(img);
n = numel(img_list);

%let user select a portion to extract
fig = figure(1); image(img);
frame = uint16(ginput(4));
close(fig);
drawnow
%pre allocate arrays
[x y] = meshgrid(1:N,zeros(n,1));
t = zeros(n,1);

%extract data from first photo
[~,y(1,:)] = max(img(:,:,2)-img(:,:,1));
t(1) = datenum(exif.DateTime,'yyyy:mm:dd HH:MM:SS');

%start parfor loop
ppm = ParforProgMon( 'Running', n );

%loop through arrays
if n>2,
    for i=2:n
        
        %load photo
        img = imread(['IMG_' num2str(img_list(i),'%04.0f') '.JPG']);
        exif =imfinfo(['IMG_' num2str(img_list(i),'%04.0f') '.JPG']);
        t(i) = datenum(exif.DateTime,'yyyy:mm:dd HH:MM:SS');
        
        %display progress
        ppm.increment();
        progressbar(i/n);
        %extract green line
        [~,y(i,:)] = max(img(:,:,2)-img(:,:,1));
    end
end

%keep only within limits
frame = [min(frame(:,1)),min(frame(:,2));...
    max(frame(:,1)),max(frame(:,2))];

y = y.*(x > frame(1) & x < frame(2));
x = x.*(x > frame(1) & x < frame(2));

x = x.*(y > frame(3) & y < frame(4));
y = y.*(y > frame(3) & y < frame(4));

x(x==0)=NaN;
y(y==0)=NaN;

%distortion parameters
k = 0.2;
angle = -2;
center = [M/2 N/2];

%converts the mesh into a colum vector of coordiantes relative to the center
xt = x - center(1);
yt = y - center(2);

%converts the x-y coordinates to polar coordinates
[theta,r] = cart2pol(xt,yt);

%calculate the maximum vector (image center to image corner) to be used for normalization
R = sqrt(center(1)^2 + center(2)^2);

%normalize the polar coordinate r to range between 0 and 1
r = r/R;

%apply the r-based transformation
s = r.*(1+k.*(r.^2));

%un-normalize s
s2 = s * R;

%rotate, while we're at it
theta = theta-deg2rad(angle);

%return to cartesian
[xt,yt] = pol2cart(theta,s2);

%convert to original coordinates
xt = xt + center(1);
yt = yt + center(2);


%perspective correction simple but fast function to change perspective and to convert to meters
xt = 0.000145.*xt;
xt = xt-0.0002*xt;
yt = -0.00013742*yt;

% center data, fit data and save result to struct
x = xt(1,:);
z = yt;
x = x-mean(x(~isnan(x)));
z = z-mean(z(~isnan(z)));

%save figure
t = inpaint_nans(t);
x = inpaint_nans(x);
z = inpaint_nans(z);
t = t-t(1);
x = x-(x(end)+x(1))/2;
imagesc(x,t,z,[-0.03 0.03]);
colorbar, axis xy, datetick('y'), axis tight;



%save as matlab file to my dropbox
save([dropbox filesep 'work' filesep 'WaveRipple' filesep 'data_lb' filesep 'IMG_' num2str(img_list(1),'%04.0f') '.mat'],'x','t','z');




end