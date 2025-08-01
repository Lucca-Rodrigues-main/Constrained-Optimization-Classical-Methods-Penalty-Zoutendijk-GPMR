clc, clear all, close all
warning('off','MATLAB:singularMatrix');

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Funcao objetivo
syms x1 x2 real
fo = @(x1,x2) 2.*x1.^2 + 2.*x2.^2 - 2.*x1.*x2 - 4.*x1 - 6.*x2;
% Restricoes
rest = @(x1,x2) [x1 + x2 - 2; x1 + 5.*x2 - 5; -x1; -x2];
A = [1 1; 1 5; -1 0; 0 -1];
b = [2; 5; 0; 0];
bd = zeros(length(A),1);
% Gradiente
% gradient(fo, [x1 x2])
go = @(x1,x2) [4*x1 - 2*x2 - 4; 4*x2 - 2*x1 - 6];

% Ponto inicial factivel
x(:,1) = [0; 0];
f(1) = fo(x(1,1), x(2,1));

T = table([], [], [], [], [], [], [], [], [], [],...
    'VariableNames', {'k','xk','fxk','gxk','I','dk','lambda',...
    'amax','ak','xk1'});
step = 1;
k = 1;
while 1
    % Calcula o gradiente no ponto
    g(:,k) = go(x(1,k), x(2,k));
    
    % Checa restricoes ativas
    if step == 1
        ativa = abs(rest(x(1,k), x(2,k))) < 1e-2;
    end
    
    % Matriz de projecao
    P = eye(length(A(ativa,:))) - A(ativa,:).' * ...
        inv(A(ativa,:) * A(ativa,:).') * A(ativa,:);
    
    % Calcula a direcao
    d(:,k) = -P * g(:,k);
    
    if 0 <= norm(d(:,k)) && norm(d(:,k)) <= 1e-2
        % Calcula os multiplicadores de Lagrange
        lambda = -inv(A(ativa,:) * A(ativa,:).') * A(ativa,:) * g(:,k);
        
        x(:,k+1) = x(:,k);
        f(k+1) = fo(x(1,k), x(2,k));
        newRow = {k, x(:,k).', f(k), g(:,k).', num2str(find(ativa).'),...
            d(:,k).', strtrim(sprintf('%.4f ',lambda.')), nan, nan, x(:,k+1).'};
        
        if all(lambda >= 0)
            break;
        else
            % Remove restricao associada ao menor valor
            temp = find(ativa);
            ativa(temp(lambda == min(lambda))) = false;
        end
        
        step = 2;
    else
        % Calculando alpha max
        alphaMax = fmincon(@(a) -a, 0, [],[],[],[],0,[],...
            @(a) mycon(a,x(:,end),d(:,end),rest),...
            optimoptions('fmincon','Display','off'));
        
        % Line search
        [alpha(k), foval] = fmincon(@(a) fo(x(1,k)+(a*d(1,k)), ...
            x(2,k)+(a*d(2,k))), 1, [-1; 1], [0; alphaMax], ...
            [],[],[],[],[], optimoptions('fmincon','Display','off'));
        
        % Calcula o novo x
        x(:,k+1) = x(:,k) + alpha(k) * d(:,k);
        f(k+1) = fo(x(1,k+1), x(2,k+1));
        
        step = 1;
        newRow = {k, x(:,k).', f(k), g(:,k).', num2str(find(ativa).'),...
            d(:,k).', '', alphaMax, alpha(k), x(:,k+1).'};
    end
    
    T = [T; newRow];
    
    k = k + 1;
end
newRow = {k, x(:,k).', f(k), g(:,k).', num2str(find(ativa).'),...
    d(:,k).', strtrim(sprintf('%.4f ',lambda.')), nan, nan, [nan nan]};
T = [T; newRow];
disp(T);

% Saida
xmin = x(:,k);
fmin = f(k);

fprintf('\n\nx0 = [%.6f %.6f]', x(1,1), x(2,1));
fprintf('\nx* = [%.6f %.6f]', xmin(1), xmin(2));
fprintf('\nf(x*) = %.6f\n', fmin);

function [c,ceq] = mycon(a,x,d,rest)
    c = rest(x(1)+a*d(1), x(2)+a*d(2));
    
    ceq = [];
end