# Constrained-Optimization-Classical-Methods-Penalty-Zoutendijk-GPMR
This repository contains a collection of MATLAB scripts that implement some of the classical optimization methods for constrained optimization models: Penalty and Barrier methods, linear and non-linear Zoutendijk and also the Gradient Projection Method of Rosen (GMPR).

---

## Algoritmos de Otimização Restrita
Esta seção do repositório aborda algoritmos projetados para resolver problemas de otimização restrita, onde o objetivo é minimizar uma função sujeito a um conjunto de restrições de igualdade e/ou desigualdade. Os métodos implementados aqui se dividem em duas categorias principais: **Métodos de Penalidade e Barreira** e **Métodos Primais**.

---

### Métodos de Penalidade e Barreira
A ideia central destes métodos é transformar um problema de otimização restrita em uma sequência de problemas de otimização irrestrita. Isso é feito através da criação de uma função auxiliar que combina a função objetivo original com termos que medem as violações das restrições.

#### 1. Método de Penalidade (Exterior)
O método de penalidade adiciona um termo à função objetivo que impõe uma "penalidade" por qualquer violação das restrições. A abordagem gera uma sequência de pontos que são tipicamente **infactíveis** (exteriores à região viável) e que convergem para a solução ótima do problema original à medida que o parâmetro de penalidade aumenta.

A função objetivo aumentada (ou auxiliar) é formulada como:

$$\theta(c, x) = f(x) + c \cdot P(x)$$

onde $f(x)$ é a função objetivo original, $c$ é o parâmetro de penalidade ($c > 0$), e $P(x)$ é a função de penalidade que satisfaz $P(x) = 0$ se $x$ for viável e $P(x) > 0$ se $x$ for inviável. Uma forma comum para restrições de desigualdade ($g_i(x) \le 0$) é a penalidade quadrática:

$$P(x) = \sum_{i=1}^{m} (\max[0, g_i(x)])^2$$

O algoritmo, conhecido como **Técnica de Minimização Irrestrita Sequencial (SUMT)**, resolve uma sequência de problemas irrestritos, aumentando o valor de $c$ a cada iteração ($c_k \to \infty$), o que força a solução a se aproximar da região viável.

##### Implementação: `Penalidade.m`
Este script implementa o método de penalidade exterior utilizando a abordagem SUMT e resolve o seguinte exemplo de problema de otimização restrita:

$$
\begin{aligned}
\text{Minimizar } & f(x) = x_1^2 + 4x_2^2 - 8x_1 - 16x_2 \\
\text{Sujeito a:} \\
& x_1 + x_2 \le 5 \\
& x_1 \le 3 \\
& x_1 \ge 0 \\
& x_2 \ge 0
\end{aligned}
$$

#### 2. Método de Barreira (Interior)
Diferente do método de penalidade, o método de barreira adiciona um termo à função objetivo que cria uma "barreira" para impedir que os pontos gerados saiam da região viável. O algoritmo, portanto, trabalha com uma sequência de pontos estritamente **factíveis** (interiores à região viável) que convergem para a solução ótima.

A função objetivo aumentada é formulada como:

$$\mathcal{B}(r, x) = f(x) + r \cdot \beta(x)$$

onde $r$ é o parâmetro de barreira ($r > 0$) e $\beta(x)$ é a função de barreira, que é finita no interior da região viável e tende ao infinito à medida que se aproxima da fronteira. Uma função de barreira comum é a função de barreira inversa:

$$\beta(x) = -\sum_{i=1}^{m} \frac{1}{g_i(x)}$$

O processo iterativo resolve uma sequência de problemas irrestritos, reduzindo o valor de $r$ a cada passo ($r_k \to 0$), permitindo que a solução se aproxime da fronteira da região viável, onde a solução ótima geralmente se encontra.

##### Implementação: `Barreira.m`
Este script implementa o método de barreira interior e resolve o mesmo exemplo de problema de otimização restrita do método de penalidade:

$$
\begin{aligned}
\text{Minimizar } & f(x) = x_1^2 + 4x_2^2 - 8x_1 - 16x_2 \\
\text{Sujeito a:} \\
& x_1 + x_2 \le 5 \\
& x_1 \le 3 \\
& x_1 \ge 0 \\
& x_2 \ge 0
\end{aligned}
$$

##### Limitações com desigualdade estrita
Um único problema a se destacar é que, dado que o problema de otimização resolvido por `fmincon` é:

$$
\begin{aligned}
    \text{Minimizar } & f(x)\\
    \text{Sujeito a:} \\
    & c(x) \leq 0, \\
    & ceq(x) = 0, \\
    & A \cdot x \leq b, \\
    & Aeq \cdot x = beq, \\
    & lb \leq x \leq ub
\end{aligned}
$$

Não há como utilizar restrições de desigualdade estrita, mas é possível utilizar algo como:

$$
c(x) \leq -\epsilon < 0
$$

Porém, na realidade o que teremos com `fmincon` é:

$$
c(x) \leq -\epsilon + \text{TolFun}
$$

Então é possível atribuir um valor mínimo para `TolFun` por meio de `optimoptions` que seja menor que $\epsilon$, porém, não é garantido que $TolFun$ será menor ou igual ao valor que escolhermos, o que pode ocasionar $g(x) \geq 0$. Entretando, isso apenas ocasiona um número maior de iterações, o valor ótimo será atingido da mesma maneira.

---

### Métodos Primais
Métodos primais são algoritmos que trabalham diretamente no problema original, realizando uma busca pela solução ótima através da região viável. Uma característica fundamental é que cada ponto gerado pelo processo é viável, e o valor da função objetivo decresce a cada iteração. Uma grande vantagem é que, caso o algoritmo seja interrompido prematuramente, o ponto atual é uma solução viável para o problema.

#### 1. Método das Direções Viáveis de Zoutendijk
Este método pertence à classe dos **Métodos de Direções Viáveis (FDM)**. Cada iteração consiste em dois passos principais:
1.  **Determinação de uma Direção:** Encontra-se uma direção de busca $d_k$ que seja simultaneamente *viável* (não sai da região viável para um passo pequeno) e *de descida* (reduz o valor da função objetivo, ou seja, $\nabla f(x_k)^T d_k < 0$). Essa direção é tipicamente encontrada resolvendo-se um subproblema de programação linear que busca alinhar $d_k$ o mais próximo possível do gradiente negativo, respeitando as restrições ativas.
2.  **Busca Unidimensional (Line Search):** Determina-se o tamanho do passo $\alpha_k$ que minimiza a função objetivo na direção $d_k$, garantindo que o novo ponto $x_{k+1} = x_k + \alpha_k d_k$ permaneça na região viável.

##### Implementação: `Zoutendijk_Linear.m`
Este script aplica o método de Zoutendijk e resolve o seguinte exemplo de problema de otimização com restrições lineares:

$$
\begin{aligned}
\text{Minimizar } & f(x) = x_1^2 + 2x_2^2 \\
\text{Sujeito a:} \\
& 4x_1 + x_2 \le 6 \\
& x_1 + x_2 = 3 \\
& x_1 \ge 0 \\
& x_2 \ge 0
\end{aligned}
$$

##### Implementação: `Zoutendijk_Nao_Linear.m`
Este script adapta o método de Zoutendijk para lidar com problemas de otimização que possuem restrições não lineares e resolve o seguinte exemplo:

$$
\begin{aligned}
\text{Minimizar } & f(x) = 5x_1^2 - 10x_1 - 10x_2 \log_{10}(x_2) \\
\text{Sujeito a:} \\
& x_1^2 + 2x_2^2 \le 4 \\
& x_1 \ge 0 \\
& x_2 \ge 0
\end{aligned}
$$

#### 2. Método de Projeção de Gradientes de Rosen (GPMR)
Este método pertence à classe dos **Métodos de Conjunto Ativo (Active Set Methods)** A ideia principal é projetar o gradiente negativo sobre a superfície definida pelas restrições ativas no ponto atual.
* A direção de busca $d_k$ é calculada como $d_k = -P \nabla f(x_k)$, onde $P$ é a matriz que projeta um vetor sobre o subespaço tangente formado pelas restrições do "conjunto de trabalho" (working set).
* Se a direção projetada $d_k$ for não nula, ela é uma direção de descida, e uma busca unidimensional é realizada para encontrar o próximo ponto.
* Se a direção projetada for nula, o ponto atual é um mínimo sobre a superfície das restrições ativas. Os multiplicadores de Lagrange são então calculados para verificar a otimalidade KKT. Se todos os multiplicadores forem não-negativos, a solução ótima foi encontrada. Caso contrário, a restrição associada ao multiplicador mais negativo é removida do conjunto de trabalho, permitindo que a próxima iteração explore uma direção que se afaste dessa fronteira para continuar diminuindo a função objetivo.

##### Implementação: `GPMR.m`
Este script implementa o Método de Projeção de Gradientes de Rosen e resolve o seguinte problema de otimização:

$$
\begin{aligned}
\text{Minimizar } & f(x) = 2x_1^2 + 2x_2^2 - 2x_1x_2 - 4x_1 - 6x_2 \\
\text{Sujeito a:} \\
& x_1 + x_2 \le 2 \\
& x_1 + 5x_2 \le 5 \\
& x_1 \ge 0 \\
& x_2 \ge 0
\end{aligned}
$$
