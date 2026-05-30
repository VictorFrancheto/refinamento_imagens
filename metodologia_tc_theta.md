# Metodologia: estimativa de $t_c(N)$ e do expoente $\theta$ por colapso de dados

Este documento descreve, em detalhe matemático, **como** o tempo pseudo-crítico $t_c(N)$ foi estimado a partir das curvas $p(t,N)$ e **como** o expoente $\theta$ do colapso de dados foi obtido. Todas as fórmulas refletem exatamente o que está implementado no notebook [plots_random_attack_refined.ipynb](plots_random_attack_refined.ipynb).

---

## 1. Hipótese de finite-size scaling

Em torno de uma transição contínua de tamanho finito, postula-se que o parâmetro de ordem $p(t,N)$ satisfaz a forma de escala

$$
p(t,N) \;=\; \mathcal{F}\!\Big(\,(t - t_c(N))\,N^{-\theta}\,\Big),
$$

onde:

- $t_c(N)$ é o **pseudo-crítico** — a posição da janela crítica para tamanho $N$;
- $\theta = 1/\nu'$ é o expoente que controla a **largura** dessa janela;
- $\mathcal{F}$ é uma **função mestra** (universal, independente de $N$).

Operacionalmente, isso significa duas tarefas independentes:

1. Localizar $t_c(N)$ para cada $N$ a partir da curva $p(\cdot, N)$.
2. Encontrar o $\theta$ que faz todas as curvas, reescaladas por $(t-t_c)N^{-\theta}$, **colapsarem** sobre $\mathcal{F}$.

---

## 2. Estimativa de $t_c(N)$ — inflexão em escala logarítmica

### 2.1 Definição

Em sistemas onde $t$ varia por várias ordens de grandeza com $N$ (ataques cumulativos, processos multiplicativos), a derivada natural é com respeito a $\ln t$. Define-se

$$
g(t,N) \;\equiv\; \frac{\partial p(t,N)}{\partial \ln t},
$$

e o estimador

$$
\boxed{\,t_c(N) \;=\; \arg\max_{t}\; g(t,N).\,}
$$

Justificativa: na vizinhança de $t_c(N)$, a derivada $g$ tem forma quase-Gaussiana,

$$
g(t,N) \;\approx\; \frac{A_N}{\sqrt{2\pi}\,\sigma_N}\exp\!\left(-\frac{(\ln t-\ln t_c(N))^2}{2\sigma_N^2}\right),
\qquad
\sigma_N \xrightarrow{N\to\infty} 0,
$$

de modo que o argmax de $g$ coincide com $\ln t_c(N)$.

### 2.2 Discretização

Dadas amostras $\{(t_i, p_i)\}_{i=1}^n$ ordenadas em $t$, calcula-se a derivada por diferenças centrais em $\ln t$:

$$
g_i \;=\; \frac{p_{i+1} - p_{i-1}}{\ln t_{i+1} - \ln t_{i-1}}, \qquad i = 2,\dots, n-1.
$$

No código isto é feito por `np.gradient(p, np.log(t))`, que aceita um eixo não-uniforme e produz diferenças centrais ponderadas pelos espaçamentos reais. O erro de truncamento por diferença central é

$$
g_i \;=\; g(t_i) \;+\; \tfrac{1}{6}\,g'''(t_i)\,(\Delta\ln t)^2 \;+\;\mathcal{O}\!\left(\frac{\sigma_p}{\Delta\ln t}\right),
$$

onde $\sigma_p$ é o ruído de Monte-Carlo em $p$. O segundo termo domina quando $\Delta\ln t$ é pequeno — daí a necessidade de suavização.

### 2.3 Suavização por Savitzky–Golay

A suavização é feita por **filtro de Savitzky–Golay** de janela ímpar $w$ e ordem polinomial $r=3$:

$$
\tilde g_i \;=\; \sum_{k=-(w-1)/2}^{(w-1)/2} c_k\, g_{i+k},
$$

onde os coeficientes $\{c_k\}$ vêm do ajuste local, por mínimos quadrados, de um polinômio de grau $r$ em $w$ pontos adjacentes. SG **preserva a posição e a amplitude do pico** melhor que uma média móvel simples, porque uma média móvel achata o pico introduzindo um viés sistemático na posição do argmax da ordem de $2\,g''(t_c)(\Delta\ln t)^2$.

A janela é **adaptativa** em $N$:

$$
w(n) \;=\; \mathrm{nearest\_odd}\Big(\min\big(w_{\max},\,\max(w_{\min},\,\lceil 0.05\,n\rceil)\big)\Big),
\qquad w_{\min}=7,\; w_{\max}=51.
$$

A motivação é que $n$ cresce com $N$ — usar janela fixa suaviza demais para $N$ pequenos e de menos para $N$ grandes.

### 2.4 Recorte de borda

Os pontos $i < k$ e $i > n-k$ (com $k = \max(2,\,\lfloor 0.05\,n\rfloor)$) são **descartados** antes de tomar o argmax. Isto evita que artefatos de extremidade — onde Savitzky–Golay extrapola implicitamente — sejam selecionados como picos espúrios. O estimador final é então

$$
\boxed{\;t_c(N) \;=\; t_{i^*}, \quad i^* \;=\; \arg\max_{i \in [k,\, n-k]} \tilde g_i.\;}
$$

---

## 3. $t_{\max}(N)$ e ajustes em lei de potência

### 3.1 Definição direta

$$
t_{\max}(N) \;=\; \max_i t_i,
$$

o último tempo observado na simulação para tamanho $N$.

### 3.2 Lei de potência log-log

Para $y\in\{t_c, t_{\max}, t_c/t_{\max}\}$, ajusta-se $y(N) = A\,N^{\alpha}$ por mínimos quadrados em escala log:

$$
\hat\alpha \;=\; \frac{\sum_i (\ln N_i - \overline{\ln N})(\ln y_i - \overline{\ln y})}{\sum_i (\ln N_i - \overline{\ln N})^2},
\qquad
\ln \hat A \;=\; \overline{\ln y} - \hat\alpha\,\overline{\ln N}.
$$

O coeficiente de determinação é

$$
R^2 \;=\; 1 \;-\; \frac{\sum_i (\ln y_i - \hat\alpha\ln N_i - \ln\hat A)^2}{\sum_i (\ln y_i - \overline{\ln y})^2}.
$$

### 3.3 Ajuste com offset (sem assumir $t_c^\infty = 0$)

Quando se quer **testar** se existe $t_c^\infty$ finito, ajusta-se o modelo de três parâmetros

$$
t_c(N) \;=\; t_\infty \;+\; a\,N^{-\theta_{tc}},
$$

por mínimos quadrados não-lineares (Levenberg–Marquardt via `scipy.optimize.curve_fit`). Os chutes iniciais vêm do ajuste log-log puro:

$$
\theta_0 = -\hat\alpha,\qquad a_0 = e^{\ln\hat A},\qquad t_{\infty,0} = \tfrac{1}{2}\min_i t_c(N_i).
$$

**Diagnóstico.** Se $t_\infty$ retorna compatível com zero (dentro do erro), o sistema **não** tem ponto crítico termodinâmico — $t_c$ diverge como potência de $N$. Se $t_\infty > 0$ com $R^2$ alto, o ajuste log-log puro é enviesado e o expoente correto vem do ajuste não-linear.

---

## 4. Expoente $\theta$ por colapso de dados

### 4.1 Princípio

Dado $\{t_c(N_k)\}$ já estimado, para cada candidato $\theta$ reescala-se cada curva:

$$
x_k(t) \;=\; (t - t_c(N_k))\,N_k^{-\theta},\qquad y_k(t) = p(t,N_k).
$$

Se a hipótese de scaling vale e $\theta$ é o expoente correto, todas as curvas $\{(x_k, y_k)\}$ devem colapsar sobre $\mathcal{F}$. A qualidade do colapso é medida pelo **espalhamento vertical** das curvas reescaladas em um grid comum de $x$.

### 4.2 Funcional de score

Restringe-se a janela crítica em $p$ via máscara $p \in [p_{\min}, p_{\max}] = [0.05,\, 0.9]$, para ignorar os platôs ($p\approx 0$ e $p\approx 1$) onde a variância seria trivialmente baixa.

Constrói-se um grid uniforme $\{x_j\}_{j=1}^{M}$ em

$$
[L,\,R] \;=\; \left[\max_k \min_t x_k(t),\;\; \min_k \max_t x_k(t)\right],
$$

— a interseção dos domínios — e interpola-se cada curva linearmente:

$$
P_{kj} \;=\; \mathrm{interp}\big(x_j;\; x_k, y_k\big).
$$

A curva média e seu gradiente são

$$
\bar P_j \;=\; \frac{1}{K}\sum_{k=1}^{K} P_{kj},
\qquad
\bar P'_j \;=\; \frac{d\bar P}{dx}\bigg|_{x_j}.
$$

Define-se o **peso** que enfatiza a janela crítica (onde $\bar P$ varia) e ignora platôs:

$$
w_j \;=\; |\bar P'_j| \;+\; \varepsilon,\qquad \varepsilon = 10^{-7}.
$$

A função objetivo é a **variância vertical média ponderada**:

$$
\boxed{\;\mathcal{S}(\theta) \;=\; \frac{\displaystyle\sum_{j=1}^{M} w_j \cdot \mathrm{Var}_k[P_{kj}]}{\displaystyle\sum_{j=1}^{M} w_j},
\quad \mathrm{Var}_k[P_{kj}] = \frac{1}{K}\sum_{k=1}^{K}(P_{kj} - \bar P_j)^2.\;}
$$

### 4.3 Robustez: descarte da pior curva

Se $K \ge 4$, calcula-se o desvio quadrático médio ponderado de cada curva ao redor da média:

$$
D_k \;=\; \frac{\sum_j w_j (P_{kj} - \bar P_j)^2}{\sum_j w_j},
$$

e descarta-se a curva com $D_k$ máximo (ataque-falso-positivo de $t_c$ ruim). $\bar P$ é recalculada e $\mathcal{S}(\theta)$ é avaliado sobre as $K-1$ curvas restantes. Isto torna o estimador robusto a um único $t_c(N)$ ruidoso.

### 4.4 Busca em duas etapas

A minimização de $\mathcal{S}$ é não-convexa. Faz-se:

**Etapa coarse.** $\mathcal{S}$ é avaliado em $300$ pontos uniformes em $\theta \in [0, 6]$. Localiza-se $\theta_1 = \arg\min \mathcal{S}$ nesse grid grosso.

**Etapa refine.** Recoloca-se um grid de $300$ pontos em uma janela estreita

$$
\big[\max(0,\,\theta_1 - 3\Delta),\;\;\min(6,\,\theta_1 + 3\Delta)\big],
\qquad \Delta = \frac{6}{299},
$$

e localiza-se

$$
\boxed{\;\theta^* \;=\; \arg\min_{\theta} \mathcal{S}(\theta).\;}
$$

### 4.5 Interpretação geométrica

Em $\theta = \theta^*$ as curvas $\{P_{k\cdot}\}$ se sobrepõem; $\mathcal{S}(\theta^*)$ aproxima a variância residual atribuível a:

- correções de scaling ($N^{-\omega}$ não capturadas),
- erro estatístico das simulações,
- erro no estimador de $t_c(N)$.

O quanto $\mathcal{S}(\theta)$ é "afiado" em torno de $\theta^*$ é uma medida qualitativa da incerteza em $\theta$.

---

## 5. Separação entre $\theta_{tc}$ e $\theta_{\text{collapse}}$

São **dois objetos matemáticos distintos**:

- $\theta_{tc}$ vem do ajuste $t_c(N) \sim N^{\theta_{tc}}$ — descreve a **posição** da janela crítica.
- $\theta_{\text{collapse}}$ vem da minimização de $\mathcal{S}(\theta)$ — descreve a **largura** da janela crítica.

Em FSS padrão, os dois coincidem ($= 1/\nu'$) se e somente se:

1. $t_c^\infty$ é finito e o ajuste é feito sobre $t_c(N) - t_c^\infty$;
2. ou $t_c^\infty = 0$ e a transição é "concentrada na origem" em $t$.

No regime em que $t_c$ diverge polinomialmente em $N$ (ataques sem ponto crítico termodinâmico), os dois expoentes podem **legitimamente diferir** e ambos devem ser reportados. Por essa razão, o código guarda os dois sob nomes distintos.

---

## 6. Resumo das fórmulas implementadas

| Quantidade | Fórmula | Função no código |
|---|---|---|
| Derivada em log | $g_i = \dfrac{p_{i+1}-p_{i-1}}{\ln t_{i+1}-\ln t_{i-1}}$ | `np.gradient(p, np.log(t))` |
| Suavização | $\tilde g = \mathrm{SG}_{w,3}(g)$ | `savgol_filter(g, w, 3)` |
| Janela adaptativa | $w = \mathrm{odd}(\mathrm{clip}(0.05n,\,7,\,51))$ | `_sg_window` |
| $t_c(N)$ | $t_{i^*},\; i^* = \arg\max_{i\in[k,n-k]} \tilde g_i$ | `estimate_tc_by_max_slope` |
| Ajuste log-log | $\hat\alpha = \dfrac{\mathrm{cov}(\ln N,\ln y)}{\mathrm{var}(\ln N)}$ | `fit_power_law` |
| Ajuste offset | $t_c = t_\infty + a N^{-\theta}$ (LM) | `fit_tc_with_offset` |
| Score colapso | $\mathcal{S}(\theta) = \dfrac{\sum_j w_j\,\mathrm{Var}_k[P_{kj}]}{\sum_j w_j}$ | `collapse_score_physical` |
| Peso | $w_j = |\bar P'_j| + \varepsilon$ | idem |
| $\theta^*$ | $\arg\min_\theta \mathcal{S}$ (coarse + refine) | `estimate_theta_by_collapse` |

---

## 7. Referências conceituais

- A ideia de minimizar variância vertical ponderada sob reescala é o método de **Bhattacharjee & Seno** (J. Phys. A 34, 6375, 2001), com variantes em Houdayer & Hartmann (PRB 70, 014418, 2004).
- Savitzky–Golay como derivador robusto: Savitzky & Golay, Anal. Chem. 36, 1627 (1964).
- A discussão sobre por que $\theta_{tc} \ne \theta_{\text{collapse}}$ em geral está em qualquer texto moderno de FSS — e.g. Privman (ed.), *Finite Size Scaling and Numerical Simulation of Statistical Systems*, World Scientific (1990).
