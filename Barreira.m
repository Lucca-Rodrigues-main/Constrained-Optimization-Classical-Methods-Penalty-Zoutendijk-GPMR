clc, clear all, close all
warning('off','MATLAB:nearlySingularMatrix');

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Funcao objetivo
syms x1 x2 real
X = [x1 x2];
fo = @(x1,x2) x1.^2 + 4.*x2.^2 - 8.*x1 - 16.*x2;
% Restricoes
rest = [x1 + x2 - 5; x1 - 3; -x1; -x2];

% Parametro de reducao
gama = 0.1;
% Parametro de parada
epsilon = 1e-3;
% Parametro de barreira
r(1) = 10;
x(1,:) = [2 2];
f(1) = fo(x(1,1), x(1,2));
% Funcao objetivo aumentada
B = @(r0, X) fo(X(1), X(2)) - r0 * ...
    sum(1 ./ double(subs(rest, [x1 x2], [X(1) X(2)])));

k = 1;
fprintf(['|   k   |      r(k)      |   x1(k)   |'...
    '   x2(k)   |   fo(X(k))   |   b(X(k))   |'...
    '   B(X(k))   |   r(k)*b(X(k))   |']);
while 1
    fprintf(['\n| %5d | %14d | %9.4f | %9.4f | %12.4f |'...
        ' %11.4e | %11.4e | %16.4e |'],...
        k,r(k),x(k,:),f(k),...
        -sum(1./double(subs(rest,[x1 x2],x(k,:)))),...
        B(r(k),x(k,:)),...
        r(k)*-sum(1./double(subs(rest,[x1 x2],x(k,:)))));
    
    % Busca restrita a g(x) <= 1e-3 + TolCon < 0
    [x(k+1,:), foval] = fmincon(@(X) B(r(k), X), x(k,:), [],[], ...
        [],[],[],[],@mycon, optimoptions('fmincon','Display','off',...
        'MaxIterations',100,'TolCon',1e-6));
    % Funcao objetivo aumentada
    B = @(r0, X) fo(X(1), X(2)) - r0 * ...
        sum(1 ./ double(subs(rest, [x1 x2], [X(1) X(2)])));
    % Novo parametro de penalidade
    r(k+1) = gama * r(k);
    % Novo valor da funcao objetivo
    f(k+1) = fo(x(k+1,1), x(k+1,2));
        
    if ~(abs(r(k+1)*-sum(1./double(subs(rest,...
            [x1 x2],x(k+1,:))))) > epsilon)
        break
    end
    
    k = k + 1;
end
fprintf(['\n| %5d | %14d | %9.4f | %9.4f | %12.4f |'...
    ' %11.4e | %11.4e | %16.4e |'],...
    k+1,r(k+1),x(k+1,:),f(k+1),...
    -sum(1./double(subs(rest,[x1 x2],x(k+1,:)))),...
    B(r(k+1),x(k+1,:)),...
    r(k+1)*-sum(1./double(subs(rest,[x1 x2],x(k+1,:)))));

% Saida
xmin = x(k+1,:);
fmin = f(k+1);

fprintf('\n\nx0 = [%.6f %.6f]', x(1,1), x(2,1));
fprintf('\nx* = [%.6f %.6f]', xmin(1), xmin(2));
fprintf('\nf(x*) = %.6f\n', fmin);

function [c,ceq] = mycon(X)
    c = [X(1) + X(2) - 5 + 1e-4; ...
        X(1) - 3 + 1e-4; ...
        -X(1) + 1e-4; ...
        -X(2) + 1e-4];
    ceq = [];
end