# Janela de transição (janela crítica)

Documento conceitual: o que é a *janela de transição* da curva $p(t,N)$, como ela aparece no nosso estimador de $t_c$, e por que sua largura controla o expoente do colapso.

---

## 1. Definição operacional

A curva $p(t,N)$ representa a fração de nós sobreviventes (ou alguma outra observável de ordem) ao longo de um ataque com $t$ passos numa rede de tamanho $N$. Tipicamente ela tem três regiões:

1. **Plateau inicial** ($p \approx 1$): o ataque ainda não tocou no esqueleto da rede.
2. **Queda abrupta**: o sistema atravessa a transição de percolação.
3. **Plateau final** ($p \approx 0$): a rede está fragmentada, restam só ataques redundantes até $t_{\max}$.

A **janela de transição** $W(N)$ é a região (2): o intervalo de tempo em torno de $t_c(N)$ onde $p$ varia de $\approx 1$ para $\approx 0$. Formalmente, fixando um par de níveis $p_{\max}, p_{\min}$ (por exemplo $0.9$ e $0.1$):

$$
W(N) \;\equiv\; t(p_{\min};N) \;-\; t(p_{\max};N),
\qquad
\bar t(N) \;\equiv\; \tfrac{1}{2}\bigl[t(p_{\min};N) + t(p_{\max};N)\bigr].
$$

$\bar t(N)$ é uma estimativa do tempo crítico; $W(N)$ é a *largura* da transição. As duas grandezas escalam de formas potencialmente **independentes** com $N$:

$$
\bar t(N) \;\sim\; N^{\theta_{tc}}, \qquad W(N) \;\sim\; N^{\theta_{W}}.
$$

---

## 2. Relação com o estimador $\arg\max\,dp/d\ln t$

O estimador implementado em `estimate_tc_by_max_slope` é

$$
\hat t_c(N) \;=\; \arg\max_t\, \frac{dp(t,N)}{d\ln t}.
$$

Em palavras: $\hat t_c$ é o instante de **inclinação máxima** de $p$ em escala log de tempo. Esse instante vive **dentro** da janela de transição — é, por construção, o centro mais inflexionado da queda. Em particular:

- Se $W(N)$ for **estreita** ($\theta_W < 0$ em fração de $t_{\max}$), o máximo de $dp/d\ln t$ é um pico **alto e bem definido** $\Rightarrow$ $\hat t_c$ tem variância pequena.
- Se $W(N)$ for **larga** (transição suave), $dp/d\ln t$ é quase chato $\Rightarrow$ qualquer ruído desloca o `argmax` $\Rightarrow$ $\hat t_c$ tem variância grande.
- Se a curva tiver **dois cotovelos** (cross-over), há dois máximos competindo e o `argmax` pode saltar de um para o outro conforme $N$ cresce. **Esse é o modo de falha observado em $C=N$** (ver [analise_tc.md](analise_tc.md) §1 e §5).

A janela de transição é, portanto, a região onde o estimador *deve* operar. Mascarar a busca a $p \in [p_{\min}, p_{\max}]$ é a forma direta de **forçar** o `argmax` a viver dentro dela.

---

## 3. Janela absoluta vs janela relativa

Há duas escalas naturais:

**Janela absoluta** $W(N)$ — medida em unidades de $t$. Para um ataque aleatório, cresce com $N$ porque o próprio $t_{\max}\sim N^{\alpha}$ cresce ($\alpha\approx 2.16$ nos seus dados).

**Janela relativa** $W(N)/t_{\max}(N)$ — fração de tempo gasta atravessando a transição. É a grandeza fisicamente significativa: diz se, no limite $N\to\infty$, a transição converge para uma *singularidade* (janela relativa $\to 0$) ou se permanece *espalhada* (janela relativa $\to$ const).

A razão $t_c/t_{\max} \sim N^{\beta}$ que aparece em [analise_tc.md](analise_tc.md) é uma proxy direta dessa janela relativa: $\beta < 0$ significa que o centro da transição encolhe relativo a $t_{\max}$.

| Conjunto | $\beta$ | Janela relativa | Leitura |
|---|---|---|---|
| $C=1$ | $-0.16$ | quase platô | transição fica até o fim do ataque; não há janela crítica fina |
| $C=2$ | $-0.71$ | encolhe | janela intermediária emerge |
| $C=3$ | $-0.98$ | $\sim N^{-1}$ | janela crítica nítida, comportamento de percolação típico |
| $C=N$ | $-4.16$ | colapsa | **artefato**: o que escala é o pico espúrio, não a transição (ver §5 de [analise_tc.md](analise_tc.md)) |

---

## 4. Por que a largura $W(N)$ é o coração do FSS

A hipótese de *finite-size scaling* (FSS) afirma que existem expoentes $\theta_{tc}$ e $\theta_W$ tais que

$$
p(t,N) \;=\; \mathcal{F}\!\left(\frac{t - t_c(N)}{W(N)}\right),
$$

com $\mathcal{F}$ uma **função universal** independente de $N$. Reescrevendo $W(N) = b\,N^{-\theta_{\text{coll}}}$ (convenção do código com $\theta_{\text{coll}} > 0$ se a janela encolhe):

$$
p(t,N) \;=\; \mathcal{F}\bigl((t - t_c(N))\,N^{-\theta_{\text{coll}}}\bigr).
$$

Esta é exatamente a forma usada por `collapse_score_physical`. O algoritmo busca o $\theta_{\text{coll}}$ que **minimiza a dispersão vertical** entre as curvas reescaladas, ponderada por $|d\bar p/dx|$ — ou seja, dá mais peso à janela de transição (onde a curva se move) e ignora os platôs (onde nada acontece).

Por isso o $\theta$ do colapso e o $\theta$ do ajuste $t_c \sim N^{\theta_{tc}}$ **medem coisas diferentes**:

- $\theta_{tc}$ é o expoente de **posição** (deslocamento do centro);
- $\theta_{\text{coll}}$ é o expoente de **largura** (estiramento horizontal necessário para sobrepor).

Para um sistema de percolação clássico bem comportado eles coincidem; quando $t_c^\infty$ não existe, podem divergir.

---

## 5. Diagnóstico visual: como reconhecer uma janela patológica

Ao olhar `plot_tc_diagnostic`:

| Sintoma | Janela diz | Provável causa |
|---|---|---|
| Linha vermelha no centro da queda, pico nítido | Janela bem definida | Estimador OK |
| Linha vermelha encostada no plateau inicial ou final | Janela violada por borda | Pico espúrio: aplicar máscara $p\in[p_{\min},p_{\max}]$ |
| Dois máximos no painel de $dp/d\ln t$ | Cross-over | Há duas escalas físicas — separá-las ou restringir o domínio |
| Pico largo e baixo | Janela larga ($\theta_W$ pequeno) | Sistema longe da criticalidade ou $N$ insuficiente |
| Pico evolui de "centro" para "borda" conforme $N$ cresce | Janela colapsa nas bordas | Modo de falha do $C=N$: o transitório domina o `argmax` |

---

## 6. Resumo de uma frase

> A **janela de transição** é o intervalo de tempo onde $p(t,N)$ realmente cai. Tudo de física relevante para criticalidade — $t_c$, $\theta_{tc}$, $\theta_{\text{coll}}$, validade do colapso, possibilidade de cross-over — vive nela; o estimador só funciona quando o `argmax` é forçado a permanecer dentro dela.

---

## 7. Referências cruzadas

- Definição matemática de $t_c$ e do estimador: §2 de [metodologia_tc_theta.md](metodologia_tc_theta.md).
- Separação $\theta_{tc}$ vs $\theta_{\text{coll}}$: §5 de [metodologia_tc_theta.md](metodologia_tc_theta.md).
- Diagnóstico do caso $C=N$ (modo de falha por borda): §5 de [analise_tc.md](analise_tc.md).
- Implementação do colapso ponderado por $|d\bar p/dx|$: `collapse_score_physical` em [plots_random_attack_refined.ipynb](plots_random_attack_refined.ipynb).
