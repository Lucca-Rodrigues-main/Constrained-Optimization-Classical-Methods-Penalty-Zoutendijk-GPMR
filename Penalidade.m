clc, clear all, close all

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Funcao objetivo
syms x1 x2 real
X = [x1 x2];
fo = @(x1,x2) x1.^2 + 4.*x2.^2 - 8.*x1 - 16.*x2;
% Restricoes
rest = [x1 + x2 - 5; x1 - 3; -x1; -x2];

% Parametro de crescimento
eta = 10;
% Parametro de parada
epsilon = 1;
% Parametro de penalidade
c(1) = 1;
x(1,:) = [4 1];
f(1) = fo(x(1,1), x(1,2));
% Checa restricoes violadas
violated = double(subs(rest,[x1 x2], x(1,:))) > 0;
% Funcao objetivo aumentada
teta = @(c0, X) fo(X(1), X(2)) + c0 * ...
    sum((max(0, double(subs(rest(violated), [x1 x2], [X(1) X(2)])))).^2);

k = 1;
fprintf(['|   k   |   c(k)   |   x1(k)   |'...
    '   x2(k)   |   fo(X(k))   |   RV   |   P(X(k))   |'...
    '   teta(X(k))   |   c(k)*P(X(k))   |']);
while 1
    fprintf(['\n| %5d | %8d | %9.4f | %9.4f | %12.4f | %6d |'...
        ' %11.4f | %14.4f | %16.4f |'],...
        k,c(k),x(k,:),f(k),sum(violated),...
        sum((max(0,double(subs(rest(violated),[x1 x2],x(k,:))))).^2),...
        teta(c(k),x(k,:)),...
        c(k)*sum((max(0,double(subs(rest(violated),[x1 x2],x(k,:))))).^2));
    
    % Busca com fmincon
    [x(k+1,:), foval] = fmincon(@(X) teta(c(k), X), x(k,:), [],[], ...
        [],[],[],[],@mycon, optimoptions('fmincon','Display','off'));
    
    % Novas restricoes violadas
    violated = double(subs(rest,[x1 x2], x(k+1,:))) > 0;
    % Funcao objetivo aumentada
    teta = @(c0, X) fo(X(1), X(2)) + c0 * ...
        sum((max(0, double(subs(rest(violated), ...
        [x1 x2], [X(1) X(2)])))).^2);
    % Novo parametro de penalidade
    c(k+1) = eta * c(k);
    % Novo valor da funcao objetivo
    f(k+1) = fo(x(k+1,1), x(k+1,2));
        
    if ~(norm(x(k+1,:) - x(k,:)) > epsilon && abs(f(k+1) - f(k)) > epsilon)
        break
    end
    
    k = k + 1;
end
fprintf(['\n| %5d | %8d | %9.4f | %9.4f | %12.4f | %6d |'...
    ' %11.4f | %14.4f | %16.4f |'],...
    k+1,c(k+1),x(k+1,:),f(k+1),sum(violated),...
    sum((max(0,double(subs(rest(violated),[x1 x2],x(k+1,:))))).^2),...
    teta(c(k+1),x(k+1,:)),...
    c(k+1)*sum((max(0,double(subs(rest(violated),[x1 x2],x(k+1,:))))).^2));

% Saida
xmin = x(k+1,:);
fmin = f(k+1);

fprintf('\n\nx0 = [%.6f %.6f]', x(1,1), x(2,1));
fprintf('\nx* = [%.6f %.6f]', xmin(1), xmin(2));
fprintf('\nf(x*) = %.6f\n', fmin);

function [c,ceq] = mycon(X)
    c = [X(1) + X(2) - 5; X(1) - 3; -X(1); -X(2)];
    ceq = [];
end