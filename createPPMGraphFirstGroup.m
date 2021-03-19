function [W] = createPPMGraphFirstGroup(alpha,p,w)
% Creates a PPM Graph: The vector alpha contains the sizes of the groups.
% Then, for each pair (v_i,v_j) of vertices, with v_i belongig to the first
% group and v_j belonging to the second group, there is a weight w created
% with probability p.
if nargin<3
    w = 1;
end
n = sum(alpha);
W = zeros(n,n);
index = 1;
for k=1:length(alpha)
    s = alpha(k);
    W(index:index+s-1,index:index+s-1) = ones(s,s);
    index = index+s;
end
% Delete diagonal entries
for k=1:length(n)
    W(k,k) = 0;
end

for j=1:alpha(1)+alpha(2)
    for k=j+1:alpha(1)+alpha(2)
        if W(j,k)==0
            x = rand;
            if x<p
                W(j,k) = w;
                W(k,j) = w;
            end
        end
    end
end
% 
end

