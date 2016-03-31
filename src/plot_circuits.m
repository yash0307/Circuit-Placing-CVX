function [] = plot_circuits(circuits_pos, circuits_size,L,xbound,ybound,varargin)
%%some code for plotting circuits
% circuits_pos is the position of the center of the components
%circuits_size is a 2 by n vector with the heigh and width of the
%components
%L is the laplacian of the connectivity graph
%xbound is the x plotting bound
%yboudn is the y plotting bound
%varargin contains the objective, violations, Kappa, trial_number, with_text
numvarargs = length(varargin);
if(numvarargs > 5)
    error('plot_circuits::TooManyArguments','required at most 4 additioanl arguments');
end
optargs = {-1,-1,-1,-1,-1};
optargs(1:numvarargs) = varargin;
[objective violation Kappa trial_number with_text] = optargs{:};
N = length(circuits_pos);

cla
hold on
for i = 1:N
    x = [circuits_pos(1,i)+circuits_size(1,i)/2; circuits_pos(1,i)+circuits_size(1,i)/2; circuits_pos(1,i)-circuits_size(1,i)/2; circuits_pos(1,i)-circuits_size(1,i)/2; circuits_pos(1,i)+circuits_size(1,i)/2];
    y = [circuits_pos(2,i)-circuits_size(2,i)/2; circuits_pos(2,i)+circuits_size(2,i)/2; circuits_pos(2,i)+circuits_size(2,i)/2; circuits_pos(2,i)-circuits_size(2,i)/2; circuits_pos(2,i)-circuits_size(2,i)/2];
    plot(x,y,'k');
end
for i = 1:N-1
    for j = i:N
        if(L(i,j) < 0)
            x = [circuits_pos(1,i) circuits_pos(1,j)];
            y = [circuits_pos(2,i) circuits_pos(2,j)];
            plot(x,y,':r');
        end
    end
end
axis equal
axis([-xbound xbound -ybound ybound]);
title1 = '';
if(objective ~= -1)
    title1 = ['objective= ' num2str(objective) ' '];
end
if(violation ~= -1)
    title1 = [title1  'violation=' num2str(violation) ' '];
end
if(Kappa ~= -1)
    title1 = [title1 'tau=' num2str(Kappa)];
end
title2 = '';
if(trial_number ~= -1)
    title2 = ['test:' num2str(trial_number) ' '];
end
title({title1,title2});
axis equal
axis([-xbound xbound -ybound ybound]);
if(with_text ~= -1)
    for i = 1:N
        text(circuits_pos(1,i),circuits_pos(2,i),num2str(i));
    end
end
drawnow

end

