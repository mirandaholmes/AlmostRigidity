% Get 3 particles to fix from adjacency matrix

% Created August 8, 2012

function pfix = get_pfix(A)

ind = find(triu(A)==1); % indices of bonds
siz = size(A);
pfix = [];

for ji=1:length(ind-2)
    
    [p1,p2] = ind2sub(siz,ind(ji)); % two particles connected together
    i3 = find(A(p2,:));  % third particle connected to second
    for p3 = i3
        if(A(p3,p1)==1)
            pfix = [p1 p2 p3];
            return;
        end
    end
end

