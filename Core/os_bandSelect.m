% Selects modes in desired bandwidth
function [p,Psi] = os_bandSelect(P,PSI,bnd,lbnd)

if (nargin < 4);lbnd = [];end

for m = 1:size(P,2)
    [I,~] = find(abs(P{1,m})/(2*pi)<bnd);
    p{1,m} = P{1,m}(I);
    Psi{1,m} = PSI{1,m}(:,I);
end

if ~isempty(lbnd)
for m = 1:size(p,2)
    [I,~] = find(abs(p{1,m})/(2*pi)>lbnd);
    p{1,m} = p{1,m}(I);
    Psi{1,m} = Psi{1,m}(:,I);
end
end

end