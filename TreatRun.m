function out = TreatRun(name,variable,range)
%save all data in cells struct



%name = 'ripple_eq2';

%variable = 'x0';

%range = 0.00015:0.000025:0.00045;


out.(variable) = range;

for i=1:length(range)
    
   t = load([name '_' variable num2str(range(i)) '.mat']);
   
   out.bedshearstress(i) = {t.bedshearstress(:,:,:)};
   disp([t.p.savename ' loaded'])
   
       
   

end

out.p = t.p;
out.iteration = t.iteration;
out.time = t.time;
out.variable = variable;
out.name = name;
if ~exist(out.p.savedir,'dir'),
    mkdir(out.p.savedir)
end


out.p.Nosc = 9;

out.p.savedir = [dropbox filesep 'work' filesep 'WaveRipple' filesep 'data_lb'];

save([out.p.savedir filesep name '.mat'],'-struct','out','p','bedshearstress','iteration','time',variable,'variable')