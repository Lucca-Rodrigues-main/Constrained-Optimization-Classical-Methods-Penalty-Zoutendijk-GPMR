clc, clear all, close all
warning('off','MATLAB:singularMatrix');

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Funcao objetivo
syms x1 x2 real
fo = @(x1,x2) x1.^2 + 2.*x2.^2;
% Restricoes
rest = @(x1,x2) [4.*x1 + x2 - 6; x1 + x2 - 3; -x1 - x2 + 3; -x1; -x2];
A = [4 1; 1 1; -1 -1; -1 0; 0 -1];
b = [6; 3; -3; 0; 0];
bd = zeros(length(A),1);
% Gradiente
% gradient(fo, [x1 x2])
go = @(x1,x2) [2*x1; 4*x2];

% Ponto inicial factivel
x(:,1) = [0; 3];
f(1) = fo(x(1,1), x(2,1));
% Direcoes
d(:,1) = x(:,1);

T = table([], [], [], [], [], [], [], [], [], [],...
    'VariableNames', {'k','xk','fxk','gxk','I','dk','gxkdk',...
    'amax','ak','xk1'});
k = 1;
while 1
    % Calcula o gradiente no ponto
    g(:,k) = go(x(1,k), x(2,k));
    
    % Checa restricoes ativas
    ativa = abs(rest(x(1,k), x(2,k))) < 1e-2;
    
    % Problema de localizacao da direcao
    [d(:,k), foval] = fmincon(@(X) (g(1,k)*X(1) + g(2,k)*X(2)), d(:,end), ...
        [A(ativa,:); -1 0; 1 0; 0 -1; 0 1], [bd(ativa); 1; 1; 1; 1], ...
        [],[],[],[],[], optimoptions('fmincon','Display','off'));
    if abs(g(:,k).'*d(:,k)) < 1e-3
        break
    end
    
    % Calculando alpha max
    bC = b(~ativa) - A(~ativa,:) * x(:,k);
    dC = A(~ativa,:) * d(:,k);
    if all(dC <= 0)
        alphaMax = inf;
    else
        dC = dC(dC > 0);
        bC = bC(dC > 0);
        alphaMax = min(bC./dC);
    end
    
    % Line search
    [alpha(k), foval] = fmincon(@(a) fo(x(1,k)+(a*d(1,k)), ...
        x(2,k)+(a*d(2,k))), 1, [-1; 1], [0; alphaMax], ...
        [],[],[],[],[], optimoptions('fmincon','Display','off'));
    
    % Calcula o novo x
    x(:,k+1) = x(:,k) + alpha(k) * d(:,k);
    f(k+1) = fo(x(1,k+1), x(2,k+1));
    
    newRow = {k, x(:,k).', f(k), g(:,k).', num2str(find(ativa).'),...
        d(:,k).', g(:,k).'*d(:,k), alphaMax, alpha(k), x(:,k+1).'};
    T = [T; newRow];
    
    k = k + 1;
end
newRow = {k, x(:,k).', f(k), g(:,k).', num2str(find(ativa).'),...
    d(:,k).', g(:,k).'*d(:,k), nan, nan, [nan nan]};
T = [T; newRow];
disp(T);

% Saida
xmin = x(:,k);
fmin = f(k);

fprintf('\n\nx0 = [%.6f %.6f]', x(1,1), x(2,1));
fprintf('\nx* = [%.6f %.6f]', xmin(1), xmin(2));
fprintf('\nf(x*) = %.6f\n', fmin);