function PlotFrame(out)

n = 150;
n = n+1;

u0 = 2*pi*0.5*out.p.d0/out.p.Tosc;

%ux = PreprocessVelocity(out.velocity.ux(:,:,n),out);
%uy = PreprocessVelocity(out.velocity.uy(:,:,n),out);
ux = out.velocity.ux(:,:,n);
uy = out.velocity.uy(:,:,n);
q = sqrt(ux.*ux + uy.*uy);
%q = PreprocessQuantity(q,out);    


% determine extreme values over entire run for color plot limits

qmax = max(q(:))+0.1;
qmin = 0;


           
%g.colorlim = [out.p.qmin out.p.qmax];


out.p.doTracers = 0;
out.p.plotVectors = 1;
out.p.Vectorspacing = 10;

q(out.p.obst==1) = NaN; % set solid nodes to NaN, which lets us draw them a different color
NaNcolor = [.75 .75 .75]; % gray
customcmap = GetCmap;
colorlim = [qmin qmax];

h =imagescNaN(out.p.x0*(0:out.p.lx-1),out.p.x0*(0:out.p.ly-1),q',customcmap,NaNcolor,colorlim); % a wrapper for imagesc that makes the first row of color scale different for NaNs

set(gca,'ydir','normal')
axis image
hold on
xint = 4;
yint = 4;
    h = quiver(out.p.x0*(0:xint:out.p.lx-1),out.p.x0*(0:yint:out.p.ly-1),ux(1:xint:out.p.lx,1:yint:out.p.ly).',uy(1:xint:out.p.lx,1:yint:out.p.ly).',0.5);
    set(h,'color','k')


end

function q = PreprocessQuantity(q,out)

% force free-slip boundaries to have same value as neighboring row/column
if strcmp(out.p.bdy.left,'freeslip')
    q(1,:) = q(2,:);
end

if strcmp(out.p.bdy.right,'freeslip')
    q(out.p.lx,:) = q(out.p.lx-1,:);
end

if strcmp(out.p.bdy.lower,'freeslip')
    q(:,1) = q(:,2);
end

if strcmp(out.p.bdy.upper,'freeslip')
    q(:,out.p.ly) = q(:,out.p.ly-1);
end
end

function u = PreprocessVelocity(u,out)

% force obstacles to have zero velocity
u(repmat(out.p.obst,[1 1 size(u,3)])==1) = 0;

% force free-slip boundaries to have same velocity as neighboring row/column
if strcmp(out.p.bdy.left,'freeslip')
    u(1,:) = u(2,:);
end

if strcmp(out.p.bdy.right,'freeslip')
    u(out.p.lx,:) = u(out.p.lx-1,:);
end

if strcmp(out.p.bdy.lower,'freeslip')
    u(:,1) = u(:,2);
end

if strcmp(out.p.bdy.upper,'freeslip')
    u(:,out.p.ly) = u(:,out.p.ly-1);
end

% pad according to boundary conditions 
u = PadBoundaries(u,out);
end

function M = PadBoundaries(M,out)

switch out.p.bdy.left % (actually the first row of the matrix)
    
    case 'noslip'
        M = cat(1,M(1,:),M); % pad with extension
        
    case 'freeslip'
        M = cat(1,M(1,:),M); % pad with extension
        
    otherwise % boundary will be periodic, the default
        M = cat(1,M(out.p.lx,:),M); 
        
end

switch out.p.bdy.right % (actually the last row of the matrix)
    
    case 'noslip'
        M = cat(1,M,M(out.p.lx,:)); % pad with extension
        
    case 'freeslip'
        M = cat(1,M,M(out.p.lx,:)); % pad with extension
        
    otherwise % boundary will be periodic, the default
        M = cat(1,M,M(1,:)); 
        
end

switch out.p.bdy.lower % (actually the first column of the matrix)
    
    case 'noslip'
        M = cat(2,M(:,1),M); % pad with extension
        
    case 'freeslip'
        M = cat(2,M(:,1),M); % pad with extension
        
    otherwise % boundary will be periodic, the default
        M = cat(2,M(:,out.p.ly),M); 
        
end

switch out.p.bdy.upper % (actually the last column of the matrix)
    
    case 'noslip'
        M = cat(2,M,M(:,out.p.ly)); % pad with extension
        
    case 'freeslip'
        M = cat(2,M,M(:,out.p.ly)); % pad with extension
        
    otherwise % boundary will be periodic, the default
        M = cat(2,M,M(:,1)); 
        
end
end

function customcmap = GetCmap


% verylightjet256

customcmap = [0.67,0.92,1;0.68,0.92,0.99;0.68,0.92,0.98;...
    0.68,0.92,0.97;0.68,0.92,0.96;0.69,0.92,0.96;...
    0.69,0.92,0.95;0.69,0.92,0.94;0.69,0.926,0.937;...
    0.701,0.927,0.929;0.704,0.928,0.921;0.706,0.928,0.913;...
    0.709,0.929,0.906;0.711,0.93,0.898;0.714,0.93,0.89;...
    0.716,0.931,0.882;0.719,0.931,0.874;0.721,0.932,0.866;...
    0.724,0.933,0.858;0.727,0.933,0.85;0.729,0.934,0.843;...
    0.732,0.935,0.835;0.734,0.935,0.827;0.737,0.936,0.819;...
    0.739,0.936,0.811;0.742,0.937,0.803;0.744,0.938,0.795;...
    0.747,0.938,0.787;0.749,0.939,0.78;0.752,0.939,0.772;...
    0.754,0.94,0.764;0.757,0.941,0.756;0.759,0.941,0.748;...
    0.762,0.942,0.74;0.765,0.943,0.732;0.767,0.943,0.724;...
    0.77,0.944,0.717;0.772,0.944,0.709;0.775,0.945,0.701;...
    0.777,0.946,0.693;0.78,0.946,0.685;0.782,0.947,0.677;...
    0.785,0.948,0.669;0.787,0.948,0.661;0.79,0.949,0.654;...
    0.792,0.949,0.646;0.795,0.95,0.638;0.797,0.951,0.63;...
    0.80,0.951,0.622;0.803,0.952,0.614;0.805,0.952,0.606;...
    0.808,0.953,0.598;0.81,0.954,0.591;0.813,0.954,0.583;...
    0.815,0.955,0.575;0.818,0.956,0.567;0.82,0.956,0.559;...
    0.823,0.957,0.551;0.825,0.957,0.543;0.828,0.958,0.535;...
    0.83,0.959,0.528;0.833,0.959,0.52;0.835,0.96,0.512;...
    0.838,0.96,0.504;0.84,0.961,0.496;0.843,0.962,0.488;...
    0.846,0.962,0.48;0.848,0.963,0.472;0.851,0.964,0.465;...
    0.853,0.964,0.457;0.856,0.965,0.449;0.858,0.965,0.441;...
    0.861,0.966,0.433;0.863,0.967,0.425;0.866,0.967,0.417;...
    0.868,0.968,0.409;0.871,0.969,0.402;0.873,0.969,0.394;...
    0.876,0.97,0.386;0.878,0.97,0.378;0.881,0.971,0.37;...
    0.884,0.972,0.362;0.886,0.972,0.354;0.889,0.973,0.346;...
    0.891,0.973,0.339;0.894,0.974,0.331;0.896,0.975,0.323;...
    0.899,0.975,0.315;0.901,0.976,0.307;0.904,0.977,0.299;...
    0.906,0.977,0.291;0.909,0.978,0.283;0.911,0.978,0.276;...
    0.914,0.979,0.268;0.916,0.98,0.26;0.919,0.98,0.252;...
    0.922,0.981,0.244;0.924,0.981,0.236;0.927,0.982,0.228;...
    0.929,0.983,0.22;0.932,0.983,0.213;0.934,0.984,0.205;...
    0.937,0.985,0.197;0.939,0.985,0.189;0.942,0.986,0.181;...
    0.944,0.986,0.173;0.947,0.987,0.165;0.949,0.988,0.157;...
    0.952,0.988,0.15;0.954,0.989,0.142;0.957,0.99,0.134;...
    0.959,0.99,0.126;0.962,0.991,0.118;0.965,0.991,0.11;...
    0.967,0.992,0.102;0.97,0.993,0.094;0.972,0.993,0.087;...
    0.975,0.994,0.079;0.977,0.994,0.071;0.98,0.995,0.063;...
    0.982,0.996,0.055;0.985,0.996,0.047;0.987,0.997,0.039;...
    0.99,0.998,0.031;0.992,0.998,0.024;0.995,0.999,0.016;...
    0.997,0.999,0.080;1,1,0;1,0.992,0;1,0.984,0;1,0.977,0;...
    1,0.969,0;1,0.961,0;1,0.953,0;1,0.945,0;1,0.938,0;...
    1,0.93,0;1,0.922,0;1,0.914,0;1,0.906,0;1,0.898,0;...
    1,0.891,0;1,0.883,0;1,0.875,0;1,0.867,0;1,0.859,0;...
    1,0.852,0;1,0.844,0;1,0.836,0;1,0.828,0;1,0.82,0;...
    1,0.813,0;1,0.805,0;1,0.797,0;1,0.789,0;1,0.781,0;...
    1,0.773,0;1,0.766,0;1,0.758,0;1,0.75,0;1,0.742,0;...
    1,0.734,0;1,0.727,0;1,0.719,0;1,0.711,0;1,0.703,0;...
    1,0.695,0;1,0.688,0;1,0.68,0;1,0.672,0;1,0.664,0;...
    1,0.656,0;1,0.648,0;1,0.641,0;1,0.633,0;1,0.625,0;...
    1,0.617,0;1,0.609,0;1,0.602,0;1,0.594,0;1,0.586,0;...
    1,0.578,0;1,0.57,0;1,0.563,0;1,0.555,0;1,0.547,0;...
    1,0.539,0;1,0.531,0;1,0.523,0;1,0.516,0;1,0.508,0;...
    1,0.50,0;1,0.492,0;1,0.484,0;1,0.477,0;1,0.469,0;...
    1,0.461,0;1,0.453,0;1,0.445,0;1,0.438,0;1,0.43,0;...
    1,0.422,0;1,0.414,0;1,0.406,0;1,0.398,0;1,0.391,0;1,0.383,0;...
    1,0.375,0;1,0.367,0;1,0.359,0;1,0.352,0;1,0.344,0;1,0.336,0;...
    1,0.328,0;1,0.32,0;1,0.313,0;1,0.305,0;1,0.297,0;1,0.289,0;...
    1,0.281,0;1,0.273,0;1,0.266,0;1,0.258,0;1,0.25,0;1,0.242,0;...
    1,0.234,0;1,0.227,0;1,0.219,0;1,0.211,0;1,0.203,0;1,0.195,0;...
    1,0.188,0;1,0.18,0;1,0.172,0;1,0.164,0;1,0.156,0;1,0.148,0;...
    1,0.141,0;1,0.133,0;1,0.125,0;1,0.117,0;1,0.109,0;1,0.102,0;...
    1,0.094,0;1,0.086,0;1,0.078,0;1,0.070,0;1,0.063,0;1,0.055,0;...
    1,0.047,0;1,0.039,0;1,0.031,0;1,0.023,0;1,0,0;1,0,0;1,0,0];
end