function PlotTauB(out)
%simple function to plot the bedshearstress on the bed surface, and to save
%the result in a png


lx = out.p.lx;
d0 = out.p.d0;
x0 = out.p.x0;    
name = out.name;



%for i=1:numel(out.bedshearstress),
    
    
%bedshearstress = out.bedshearstress{i};
N = min(5,out.p.Nosc);
Tcycle = out.p.Tosc/(out.p.saveint*out.p.t0);

%do cycle-averaged bedshearstress
taub = mean(reshape(out.bedshearstress(:,:,end-(N*Tcycle)+1:end),lx,Tcycle,N),3)';
pcolor(taub),shading flat
set(gca,'CLim',[-1 1]),colormapeditor

[~,taub_grad] = gradient(bedshearstress(:,1,end-(N*Tcycle)+1:end),out.p.x0);
taub_grad = mean(reshape(taub_grad,lx(i),Tcycle,N),3)';

fig = figure;
s1 = subplot(1,1,1,'Parent',fig,'Layer','top','CLim',[-0.2 0.2]);

imagesc((1:lx(i))*x0,out.time(1:out.p.Tosc/(out.p.saveint*out.p.t0)),taub,'Parent',s1,'CDataMapping','scaled');
%imagesc((1:out.p.lx)*out.p.x0,out.time(1:out.p.Tosc/(out.p.saveint*out.p.t0)),taub_grad,'Parent',s2,'CDataMapping','scaled');
shading flat
box on;

set(s1,'Layer','top','CLim',[-0.2 0.2]);


axis(s1,'xy')


xlabel(s1,'Length along bed (m)','HorizontalAlignment','left');
ylabel(s1,'Time (s)','HorizontalAlignment','right');
title(s1,[name out.variable num2str(out.(out.variable)(i)) '  \tau_b [Nm^{-2}]']);


colorbar('peer',s1,'location','East');
%colorbar('peer',s2,'location','East');
colormap ([0,0,0.5625;0,0,0.625;0,0,0.6875;0,0,0.750;0,0,0.8125;...
    0,0,0.875;0,0,0.937;0,0,1;0,0.04166,1;0,0.0833,1;...
    0,0.125,1;0,0.16666,1;0,0.208,1;0,0.250,1;0,0.29166,1;...
    0,0.333,1;0,0.375,1;0,0.4166,1;0,0.4583,1;0,0.500,1;...
    0,0.5416,1;0,0.583,1;0,0.625,1;0,0.666,1;0,0.70833,1;...
    0,0.750,1;0,0.7916,1;0,0.833,1;0,0.875,1;0,0.916682,1;0,0.958333,1;...
    0,1,1;1,1,0;1,0.9565,0;1,0.913,0;1,0.8695651,0;1,0.8260,0;1,0.782,0;...
    1,0.7391,0;1,0.6956,0;1,0.65,0;1,0.6086,0;1,0.5652,0;1,0.52173,0;...
    1,0.4782,0;1,0.4347,0;1,0.3913,0;1,0.34782,0;1,0.3043,0;1,0.260,0;...
    1,0.2173,0;1,0.1739,0;1,0.13043,0;1,0.0869,0;1,0.04347,0;1,0,0;0.9375,0,0;...
    0.875,0,0;0.8125,0,0;0.750,0,0;0.6875,0,0;0.625,0,0;0.5625,0,0;0.500,0,0;]);
saveas(fig,[dropbox filesep 'work' filesep 'WaveRipple' filesep 'data_lb' filesep name '_' out.variable num2str(out.(out.variable)(i)) '_taub_av.png'],'png')

close all;

%end

%% Statistics
%ripple_lx2
clearvars -except out
%lx = out.p.lx;
hold off
lx = out.(out.variable);
%Nr = out.(out.variable);
Nr = 1;
d0 = out.p.d0;
x0 = out.p.x0;
%x0 = out.(out.variable);
name = 'ripple';
%names = {'ripple_eq2_x00.00015.mat','ripple_eq2_x00.0002.mat','ripple_eq2_x00.00025.mat','ripple_eq2_x00.0003.mat','ripple_eq2_x00.00035.mat'};
%hold on
%y = zeros(length(out.(out.variable)),1);


for i=1:numel(out.bedshearstress),
    

bedshearstress = out.bedshearstress{i};

%%


n = size(bedshearstress,1);
if strcmpi(out.variable,'Nr'),
    cresti = round(n/Nr(i)/2:(n/Nr(i)):n);
elseif strcmpi(out.variable,'lx') && rem(lx(i),2),
    cresti = round(n/2);
else,
    cresti = n/2;
end
%%
%bedshearstress = out.bedshearstress;
N = min(5,out.p.Nosc);
Tcycle = round(out.p.Tosc/(out.p.saveint*out.p.t0));
%do cycle-averaged bedshearstress
taub = mean(reshape(bedshearstress(:,:,end-(N*Tcycle)+1:end),lx(i),Tcycle,N),3)';


[~,taub_grad] = gradient(bedshearstress(:,1,end-(N*Tcycle)+1:end),out.p.x0);
taub_grad = mean(reshape(taub_grad,lx(i),Tcycle,N),3)';
%%
l = taub(1:(Tcycle/2),:);
r = taub([1:(Tcycle/2)]+(Tcycle/2),:);
%y(i) = sum(sum(l<0))+sum(sum(r>0));
l1 = abs(sum(l(l(:)<0)));
r1 = abs(sum(r(r(:)>0)));
y(i) = l1+r1;

%y(i) = sum(sum(abs(taub_grad)));
%y(i) = mean(mean(abs(taub)));

%y(i) = abs(mean(mean(taub(1:Tcycle/2,:))))+ abs(mean(mean(taub(Tcycle/2:end,:))));
%y(i) = sum(bedshearstress(:));
%y(i) = mean(mean(taub_grad(:,cresti)));
% use interp1!
%y(i) = mean(interp1(mean(taub_grad),cresti,'pchip'));

end
%y=y/(N*Tcycle);
%y=y/(out.p.x0*out.p.lx);

%yc = sqrt(-y/out.p.rho0)*out.p.Tosc/(pi*out.p.d0);
yc = y./(lx*Tcycle);
scatter((x0.*lx)/(Nr.*out.p.d0),yc), grid on
hold on
plot([1.53 1.53],get(gca,'YLim'),'r')
%%
% %blub = [min(taub1,[],2), max(taub1,[],2)];
% bar(abs(blub))
% %plot(taub(:,[50 150 255 345 450])))

%plot(mean(taub))
%plot(mean(taub)), grid on
%blub = mean(taub);
%scatter(out.p.d0,blub)


%plot(out.p.FBdy*2/out.p.rho0/((2*pi/out.p.Tosc)^2),abs(blub1)-abs(blub2))

%set(gca,'YLim',[-0.1 0.1])

%end


% 

fig = figure;
s1 = subplot(1,1,1,'Parent',fig,'Layer','top','CLim',[-0.2 0.2]);

%s2 = subplot(2,1,2,'Parent',fig,'Layer','top','CLim',[-50 50]);
%do reshape from 4th cycle until last cycle. cycles are 40 saved timesteps

imagesc((1:out.p.lx)*x0,out.time(1:out.p.Tosc/(out.p.saveint*out.p.t0)),taub,'Parent',s1,'CDataMapping','scaled');
%imagesc((1:out.p.lx)*out.p.x0,out.time(1:out.p.Tosc/(out.p.saveint*out.p.t0)),taub_grad,'Parent',s2,'CDataMapping','scaled');
shading flat
box on;

set(s1,'Layer','top','CLim',[-0.2 0.2]);
%set(s2,'Layer','top','CLim',[-50 50]);

axis(s1,'xy')
%axis(s2,'xy')

xlabel(s1,'Length along bed (m)','HorizontalAlignment','left');
ylabel(s1,'Time (s)','HorizontalAlignment','right');
title(s1,[name out.variable num2str(out.(out.variable)(i)) '  \tau_b [Nm^{-2}]']);

%xlabel(s2,'Length along bed (m)','HorizontalAlignment','left');
%ylabel(s2,'Time (s)','HorizontalAlignment','right');
%title(s2,'Bed Shear Stress Gradient [Nm^{-2}m^{-1}]');

colorbar('peer',s1,'location','East');
%colorbar('peer',s2,'location','East');
colormap ([0,0,0.5625;0,0,0.625;0,0,0.6875;0,0,0.750;0,0,0.8125;...
    0,0,0.875;0,0,0.937;0,0,1;0,0.04166,1;0,0.0833,1;...
    0,0.125,1;0,0.16666,1;0,0.208,1;0,0.250,1;0,0.29166,1;...
    0,0.333,1;0,0.375,1;0,0.4166,1;0,0.4583,1;0,0.500,1;...
    0,0.5416,1;0,0.583,1;0,0.625,1;0,0.666,1;0,0.70833,1;...
    0,0.750,1;0,0.7916,1;0,0.833,1;0,0.875,1;0,0.916682,1;0,0.958333,1;...
    0,1,1;1,1,0;1,0.9565,0;1,0.913,0;1,0.8695651,0;1,0.8260,0;1,0.782,0;...
    1,0.7391,0;1,0.6956,0;1,0.65,0;1,0.6086,0;1,0.5652,0;1,0.52173,0;...
    1,0.4782,0;1,0.4347,0;1,0.3913,0;1,0.34782,0;1,0.3043,0;1,0.260,0;...
    1,0.2173,0;1,0.1739,0;1,0.13043,0;1,0.0869,0;1,0.04347,0;1,0,0;0.9375,0,0;...
    0.875,0,0;0.8125,0,0;0.750,0,0;0.6875,0,0;0.625,0,0;0.5625,0,0;0.500,0,0;]);
saveas(fig,[name out.variable num2str(out.(out.variable)(i)) '_taub_av.png'],'png')

%saveas(fig,[out.p.savename '_taub_av.fig'],'fig')
close all;
%end
%% Figure 3: Average BedShearStress
%names = {'ripple6_d00.06.mat','ripple6_d00.08.mat','ripple6_d00.1.mat','ripple6_d00.12.mat','ripple6_d00.14.mat','ripple6_d00.16.mat','ripple6_d00.18.mat','ripple6_d00.2.mat'};
%hold on
%average
%plot(squeeze(mean(out.bedshearstress,3)),'k')
%average flow +
%plot(squeeze(mean(out.bedshearstress(:,:,diff(out.p.Fmag(out.iteration(2:end)))>0),3)),'r')
%average flow -
%plot(squeeze(mean(out.bedshearstress(:,:,diff(out.p.Fmag(out.iteration(2:end)))<0),3)),'b')


%% algorithm to find Lsep
%ripple_lx2
clearvars -except out
time = 1:numel(out.time);
hold on;
N = min(10,out.p.Nosc);
Tcycle = out.p.Tosc/(out.p.saveint*out.p.t0);
d0 = out.p.d0;
x0 = out.p.x0;
lx = out.(out.variable);

for j=1:numel(out.bedshearstress),
    mid = lx(j)/2;
    bedshearstress = out.bedshearstress{j};
   % bedshearstress = mean(reshape(bedshearstress(:,:,end-(N*Tcycle)+1:end),out.p.lx,Tcycle,N),3);
   shear  = squeeze(bedshearstress([(mid+1):end 1:mid],:));
   crest_shear = smooth(shear(end,:),1);
crest_shear(crest_shear == 0) = 0.00001;
idx = time(crest_shear(1:end-1).*crest_shear(2:end)<0);
%idx = Tcycle*[0.5:0.5:out.p.Nosc]-1;
%idx = [15 30];


for i=1:length(idx)

    if idx(i)<2,
        continue,
    end
    diff_at_i = shear(1:end-1,idx(i)-1).*shear(2:end,idx(i)-1);
    if diff(crest_shear([idx(i)-1 idx(i)]))>0,
        blub = x0*(lx(j)-find(diff_at_i<0,1,'first'));
        %disp(['a ' num2str(blub)])
    else,
        blub = x0*find(diff_at_i<0,1,'last');
        %disp(['b ' num2str(blub)])
    end


    if isempty(blub),
        Lsep(j,i) = lx(j)*x0;
    else,
        Lsep(j,i) = blub;
    end
    %LsepL = find(shear(1:end-1,cup(1)).*shear(2:end,cup(1))<0,1,'last');
end

Lsep(Lsep==0) = NaN;
%scatter(ones(size(Lsep(j,:))).*(lx(j)*x0)./d0,Lsep(j,:)./d0,'r','.')


end
%plot(d0./lx./x0/2,nanmean(Lsep,2)./(lx.*x0/2)','r-')

%plot(d0./(lx*x0),nanmean(Lsep,2)./(lx.*x0./2)','r-')
LsepMedian = nanmedian(Lsep,2)';
plot(lx*x0/d0,(LsepMedian)./(d0),'+-k')


%plot([0.65 0.65],get(gca,'YLim'),':b')
%errorbar(out.d0,nanmean(Lsep,2),nanstd(Lsep,1,2))

%% algorithm to find Lsep for Bump
clearvars -except out
time = 1:numel(out.time);
hold on;
N = min(10,out.p.Nosc);
Tcycle = out.p.Tosc/(out.p.saveint*out.p.t0);
d0 = out.p.d0;
x0 = out.p.x0;

peaks = [1, 244, 313, 420, 554 out.p.lx]; %for bump329_bed100
bedshearstress = out.bedshearstress;
for j=2:5

       % bedshearstress = mean(reshape(bedshearstress(:,:,end-(N*Tcycle)+1:end),out.p.lx,Tcycle,N),3);
       shear  = squeeze(bedshearstress([(peaks(j)+1):peaks(j+1) peaks(j-1):peaks(j)],:));
       crest_shear = smooth(shear(end,:),1);
       crest_shear(crest_shear == 0) = 0.00001;
       idx = time(crest_shear(1:end-1).*crest_shear(2:end)<0);
%idx = Tcycle*[0.5:0.5:out.p.Nosc]-1;
%idx = [15 30];


for i=1:length(idx)

    if idx(i)<2,
        continue,
    end
    diff_at_i = shear(1:end-1,idx(i)-1).*shear(2:end,idx(i)-1);
    if diff(crest_shear([idx(i)-1 idx(i)]))>0,
        blub = x0*(out.p.lx-find(diff_at_i<0,1,'first'));
        %disp(['a ' num2str(blub)])
    else,
        blub = x0*find(diff_at_i<0,1,'last');
        %disp(['b ' num2str(blub)])
    end


    if isempty(blub),
        Lsep(j,i) = NaN;
    else,
        Lsep(j,i) = blub;
    end
    %LsepL = find(shear(1:end-1,cup(1)).*shear(2:end,cup(1))<0,1,'last');
end

Lsep(Lsep==0) = NaN;
scatter(ones(size(Lsep(j,:))).*peaks(j),Lsep(j,:),'r','.')

end
%plot(out.p.lx.*x0./d0/2,nanmedian(Lsep,2)./(out.p.lx.*x0/2)','r-')
%plot([0.65 0.65],get(gca,'YLim'),':b')
%errorbar(out.d0,nanmean(Lsep,2),nanstd(Lsep,1,2))

%% check_asym

%% Statistics
clearvars -except out
lx = out.p.lx;
d0 = out.p.d0;
x0 = out.p.x0;



u0 = pi*d0/out.p.Tosc
    

bedshearstress = out.bedshearstress;

N = min(30,out.p.Nosc);

%N=9;
%[~,taub_grad] = gradient(bedshearstress(:,1,2:361),x0);
%taub = mean(reshape(bedshearstress(:,1,2:361),lx,Tcycle,N),3)';
Tcycle = out.p.Tosc/(out.p.saveint*out.p.t0);
%do cycle-averaged bedshearstress
taub = mean(reshape(bedshearstress(:,:,end-(N*Tcycle)+1:end),lx,Tcycle,N),3)';


[~,taub_grad] = gradient(bedshearstress(:,1,end-(N*Tcycle)+1:end),x0);
taub_grad = mean(reshape(taub_grad,lx,Tcycle,N),3)';
plot(mean(taub_grad))

%% algorithm to find Lsep for grid check

savemat = NaN(8,50);
nr = 50:10:120;

for j=1:numel(nr),
    

blub=0;
idx = 0;



out = load(['grid_check_lx' num2str(nr(j)) '.mat']);



time = 1:numel(out.time);
hold on;
N = min(10,out.p.Nosc);

Tcycle = out.p.Tosc/(out.p.saveint*out.p.t0);
d0 = out.p.d0;
x0 = out.p.x0;

    mid = out.p.lx/2;
    bedshearstress = out.bedshearstress(:,1,end-(N*Tcycle)+1:end);
   %bedshearstress = mean(reshape(bedshearstress(:,:,end-(N*Tcycle)+1:end),out.p.lx,Tcycle,N),3);
   shear  = squeeze(bedshearstress([(mid+1):end 1:mid],:));
   crest_shear = smooth(shear(end,:),1);
crest_shear(crest_shear == 0) = 0.00001;
idx = time(crest_shear(1:end-1).*crest_shear(2:end)<0);
%idx = Tcycle*[0.5:0.5:out.p.Nosc]-1;
%idx = [15 30];

%pcolor(bedshearstress), shading flat

%
for i=1:length(idx)

    if idx(i)<2,
        continue,
    end
    diff_at_i = shear(1:end-1,idx(i)-1).*shear(2:end,idx(i)-1);
    if diff(crest_shear([idx(i)-1 idx(i)]))>0,
        blub = x0*(out.p.lx-find(diff_at_i<0,1,'first'));
        %disp(['a ' num2str(blub)])
    else,
        blub = x0*find(diff_at_i<0,1,'last');
        %disp(['b ' num2str(blub)])
    end


    if isempty(blub),
        Lsep(i) = NaN; %out.p.lx*x0;
    else,
        Lsep(i) = blub;
    end
    %LsepL = find(shear(1:end-1,cup(1)).*shear(2:end,cup(1))<0,1,'last');
end

Lsep(Lsep==0) = NaN;
%plot(Lsep(:),'r-')

scatter(ones(size(Lsep))*j,Lsep)
savemat(j,1:numel(Lsep)) = Lsep;
savematmean = mean(Lsep);
end