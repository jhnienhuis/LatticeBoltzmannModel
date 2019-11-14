function DoWrite(p,g)

varnames = 'p';
structname = 'g';
fieldnames = '''iteration'',''time''';

for i = 1:length(p.savequant)
    fieldnames = [fieldnames ',''' p.savequant{i} '''']; % Couldn't use this because velocity has two components with different names
end

if ~exist(p.savedir,'dir')
    mkdir(p.savedir)
end

filename  = [p.savedir filesep p.savename '.mat'];
eval(['save(filename,''-struct'',structname,' fieldnames ',''-v7.3'')']) % save specified fields of struct
save(filename,'-append',varnames) % append additional non-struct variables

disp(['Writing to disk: ', p.savename ,num2str(g.n(end))])