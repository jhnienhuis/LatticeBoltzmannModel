function [p,g] = RunLBM(p,g)

if isfield(p,'plotint') && p.plotint
    figure('color','w')
    DoPlot(p,g);
end

for n = 1:p.N

    g.n = n;

    % Flow
    g = StepLBM(p,g);
    
    % If it's time to do so:
    
    % draw plot
    if isfield(p,'plotint') && p.plotint && ~mod(n,p.plotint)
        DoPlot(p,g);
    end
    
    
    if isfield(p,'saveint') && p.saveint 
        
        % save solution
        if ~mod(n,p.saveint)
            g = DoSave(p,g);
        end
        
        % write to disk        
        if isfield(p,'writeint') && p.writeint && ~mod(n,p.writeint)
            DoWrite(p,g);
        end
        
    end

end
