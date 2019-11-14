
function LBM_ParameterRun
% LBM_ParameterRun


%which parameters
%q.eta = [0.10 0.15 0.2];
q.u0 = [0.6];
q.wlength = [100];

%put parameter in struct
q_names = fieldnames(q);
q_cell = struct2cell(q);
sizes = cellfun('prodofsize', q_cell);

%if ~matlabpool('size'), matlabpool(min(12,sizes)), end

%parallel for-loop
for i = 1:prod(sizes)
    [a] = ind2subv(sizes,i);
    [p] = struct_fun(q_names,q_cell,a);
    
    LBM(p)
    
end

%if matlabpool('size'), matlabpool('close'), end

end

function [p] = struct_fun(q_names,q_cell,a)

for j=1:length(q_names), p.(q_names{j}) = q_cell{j}(a(j));  end

end

function v = ind2subv(siz,ndx)
[out{1:length(siz)}] = ind2sub(siz,ndx);
v = cell2mat(out);
end