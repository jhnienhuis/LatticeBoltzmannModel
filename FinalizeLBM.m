function [p,g] = FinalizeLBM(p,g)

if isfield(p,'saveint') && p.saveint % if we're saving solutions    
    DoWrite(p,g); % write to disk
end

if isfield(p,'plotint') && p.plotint % if we're plotting solution
    DoPlot(p,g); % draw plot
end
