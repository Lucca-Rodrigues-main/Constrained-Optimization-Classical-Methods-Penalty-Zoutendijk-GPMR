clc, clear all, close all
warning('off','MATLAB:singularMatrix');

% Lucca Rodrigues Pinto
% https://github.com/Lucca-Rodrigues-main

% Funcao objetivo
syms x1 x2 real
fo = @(x1,x2) 5.*x1.^2 - 10.*x1 - 10.*x2.*log10(x2);
% Restricoes
rest = @(x1,x2) [x1.^2 + 2.*x2.^2 - 4; -x1; -x2];
A = [x1 2.*x2; -1 0; 0 -1];
b = [4; 0; 0];
bd = zeros(length(A),1);
% Gradiente
% gradient(fo, [x1 x2])
go = @(x1,x2) [10*x1 - 10; -(10*(log(x2) + 1))/log(10)];
grest = @(x1,x2) [2*x1 -1 0; 4*x2 0 -1];

% Ponto inicial factivel
x(:,1) = [0; sqrt(2)];
f(1) = fo(x(1,1), x(2,1));
% Direcoes
d(:,1) = x(:,1);

T = table([], [], [], [], [], [], [], [], [], [],...
    'VariableNames', {'k','xk','fxk','gxk','I','dk','zk',...
    'amax','ak','xk1'});
k = 1;
while 1
    % Calcula o gradiente no ponto
    g(:,k) = go(x(1,k), x(2,k));
    
    % Checa restricoes ativas
    ativa = abs(rest(x(1,k), x(2,k))) < 1e-2;

    % Problema de localizacao da direcao
    dz = fmincon(@(Dz) Dz(3), [d(:,end).' 0], [],[],...
        [],[],[],[], @(Dz) mycon(Dz,x(:,end),ativa,go,grest),...
        optimoptions('fmincon','Display','off'));
    d(:,k) = dz(1:2);
    z(k) = dz(3);
    
    if 0 <= abs(z(k)) && abs(z(k)) <= 1e-3
        break
    end
    
    % Calculando alpha max
    alphaMax = fmincon(@(a) -a, 0, [],[],[],[],0,[],...
        @(a) mycon2(a,x(:,end),d(:,end),rest),...
        optimoptions('fmincon','Display','off'));
    
    % Line search
    [alpha(k), foval] = fmincon(@(a) fo(x(1,k)+(a*d(1,k)), ...
        x(2,k)+(a*d(2,k))), 1, [-1; 1], [0; alphaMax], ...
        [],[],[],[],[], optimoptions('fmincon','Display','off'));
    
    % Calcula o novo x
    x(:,k+1) = x(:,k) + alpha(k) * d(:,k);
    f(k+1) = fo(x(1,k+1), x(2,k+1));
    
    newRow = {k, x(:,k).', f(k), g(:,k).', num2str(find(ativa).'),...
        d(:,k).', z(k), alphaMax, alpha(k), x(:,k+1).'};
    T = [T; newRow];
    
    k = k + 1;
end
newRow = {k, x(:,k).', f(k), g(:,k).', num2str(find(ativa).'),...
    d(:,k).', z(k), nan, nan, [nan nan]};
T = [T; newRow];
disp(T);

% Saida
xmin = x(:,k);
fmin = f(k);

fprintf('\n\nx0 = [%.6f %.6f]', x(1,1), x(2,1));
fprintf('\nx* = [%.6f %.6f]', xmin(1), xmin(2));
fprintf('\nf(x*) = %.6f\n', fmin);

function [c,ceq] = mycon(dz,x,ativa,go,grest)
    temp = grest(x(1), x(2));
    
    c = [go(x(1), x(2)).'*[dz(1); dz(2)] - dz(3);... % fun obj
        temp(:,ativa).'*[dz(1); dz(2)] - dz(3);... % rest
        -dz(1)-1; dz(1)-1; -dz(2)-1; dz(2)-1]; % etc
    
    ceq = [];
end
function [c,ceq] = mycon2(a,x,d,rest)
    c = rest(x(1)+a*d(1), x(2)+a*d(2));
    
    ceq = [];
end