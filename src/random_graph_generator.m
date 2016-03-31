function [L] = random_graph_generator(n,p)
%This function generates the graph Laplacian for a random Erdos-Renyi Graph
%with n nodes, and parameter p (the probability that any given edge exists)
L = zeros(n,n);
for i = 1:n-1
    for j = i:n
        if(rand(1) < p)
            L(i,i) = L(i,i) + 1;
            L(j,j) = L(j,j) + 1;
            L(i,j) = L(i,j) - 1;
            L(j,i) = L(j,i) - 1;
        end
    end
end

end

