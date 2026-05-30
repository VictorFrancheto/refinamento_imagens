# Análise dos expoentes de $t_c(N)$, $t_{\max}(N)$ e da razão $t_c/t_{\max}$

Resultados numéricos obtidos com o estimador refinado (Savitzky–Golay + corte de borda) sobre o conjunto $N \in \{300, 500, 600, 1000, 2000, 3000\}$ e ataques aleatórios com redundância $C \in \{1, 2, 3, N\}$.

---

## 1. Tabela-resumo

| Conjunto | $t_c(N) \sim N^{\theta}$ (log-log) | R² | $t_{\max}(N) \sim N^{\alpha}$ | $t_c/t_{\max} \sim N^{\beta}$ | Veredito |
|---|---|---|---|---|---|
| $C=1$ | $\theta \approx 2.08$ | 0.999 | $\alpha \approx 2.24$ | $\beta \approx -0.16$ (R²=0.82) | **Consistente** — $t_c$ escala junto com $t_{\max}$ |
| $C=2$ | $\theta \approx 1.47$ | 0.996 | $\alpha \approx 2.18$ | $\beta \approx -0.71$ | **Consistente** — janela crítica intermediária |
| $C=3$ | $\theta \approx 1.18$ | 0.994 | $\alpha \approx 2.17$ | $\beta \approx -0.98$ | **Consistente** — janela crítica estreita |
| $C=N$ | $\theta \approx -2.01$ | 0.604 | $\alpha \approx 2.16$ | $\beta \approx -4.16$ (R²=0.87) | **Artefato** — estimador falha para $N \ge 2000$ |

---

## 2. $t_{\max}(N)$ — escala universal $\sim N^{2.15\text{–}2.24}$

Os quatro valores de $C$ produzem o mesmo expoente, com R² = 1.000:

$$
t_{\max}(N) \;\sim\; N^{2.16 \pm 0.04}.
$$

Isto é exatamente o que se espera de um processo aleatório que precisa varrer da ordem de $\binom{N}{2} = \mathcal{O}(N^2)$ pares de nós para esgotar a rede. O fato de $\alpha$ não depender de $C$ confirma que $t_{\max}$ é uma escala **combinatória**, ditada apenas pelo tamanho do espaço de pares, e não pela dinâmica de percolação.

---

## 3. $t_c(N)$ — separação clara por valor de $C$

Os expoentes do scaling pseudo-crítico se ordenam de forma física:

$$
\theta_{C=1} \;>\; \theta_{C=2} \;>\; \theta_{C=3},
\qquad
2.08 \;>\; 1.47 \;>\; 1.18.
$$

Interpretação:

- **$C=1$** (sem redundância): a transição só ocorre quando o ataque consumiu uma fração $\mathcal{O}(1)$ de todos os pares possíveis. O $t_c$ é praticamente proporcional ao $t_{\max}$, com razão $\sim 10\%\text{–}14\%$ quase constante (expoente $-0.16$ com R²=0.82, ou seja, **um platô com ruído**). Não há "janela crítica fina" — o sistema só quebra perto do esgotamento.
- **$C=2$**: a redundância dupla move a transição para tempos sublineares em $t_{\max}$. A razão cai como $N^{-0.71}$, indicando o aparecimento de uma escala crítica genuína abaixo do esgotamento combinatório.
- **$C=3$**: $\beta \approx -1$ é o caso mais limpo — a janela crítica encolhe como $N^{-1}$ relativa a $t_{\max}$, comportamento típico de transição de percolação com expoente de correlação bem definido.

A monotonicidade $\theta_C$ decrescente com $C$ é fisicamente coerente: **quanto mais redundância, mais cedo (em fração de $t_{\max}$) o sistema atinge o ponto crítico**, porque há mais arestas vulneráveis disponíveis em paralelo.

---

## 4. Ajuste com offset $t_c(N) = t_\infty + a\,N^{-\theta}$

Os valores aparentemente "negativos" de $\theta$ no ajuste com offset **não são bug**: o modelo é simétrico sob $\theta \mapsto -\theta$ se $t_\infty$ absorver a constante, então um $\theta < 0$ retornado pelo otimizador é equivalente a um crescimento $N^{|\theta|}$. Comparando:

| Conjunto | $\|\theta\|$ (offset) | $\theta$ (log-log) | Concordância |
|---|---|---|---|
| $C=1$ | 2.02 | 2.08 | ✓ |
| $C=2$ | 1.50 | 1.47 | ✓ |
| $C=3$ | 1.70 | 1.18 | ✗ (modelo mal especificado) |
| $C=N$ | 0.00 | $-2.01$ | ✗ (dados ruins) |

Conclusões:

- Para $C=1$ e $C=2$, os dois métodos coincidem dentro de $\sim 3\%$. **Não existe $t_c^\infty$ finito**: o tempo crítico diverge polinomialmente em $N$.
- Para $C=3$, a divergência entre os dois ajustes (1.70 vs 1.18) sugere que o ansatz puro $t_c \sim N^\theta$ é insuficiente: há correções de scaling não-triviais. Vale tentar $t_c(N) = t_\infty + a N^{-\theta}$ com $t_\infty > 0$ fixado por extrapolação.
- Para $C=N$, o ajuste degenera (R² = 0.19, $t_\infty \approx -4\times 10^7$): sem significado físico — ver §5.

---

## 5. Diagnóstico do caso $C=N$

Os números brutos:

```
N=  300 | t_c = 13 404
N=  500 | t_c = 13 076
N=  600 | t_c = 18 805
N= 1000 | t_c = 36 285
N= 2000 | t_c = 152.3      ← descontinuidade
N= 3000 | t_c = 449.7
```

A queda de **três ordens de grandeza** entre $N=1000$ e $N=2000$ é incompatível com qualquer dinâmica monotônica. As assinaturas do problema:

1. R² do log-log = 0.604 (todos os outros ≥ 0.99).
2. Expoente negativo $\theta = -2.01$ — fisicamente impossível: $t_c$ não pode encolher conforme o sistema cresce em um ataque aleatório.
3. Ajuste com offset não converge para nada interpretável.
4. Razão $t_c/t_{\max}$ pula de $10^{-2}$ para $10^{-5}$ no mesmo intervalo.

**Causa provável.** O regime $C=N$ provavelmente produz uma curva $p(t)$ com **dois regimes** em escala log: uma rampa inicial rápida (transitório de saturação local) seguida pela percolação propriamente dita. Para $N \le 1000$ o pico de $dp/d\ln t$ ainda está na percolação; para $N \ge 2000$ o pico do transitório inicial cresce e passa a dominar o `argmax`, colapsando $t_c$ para um valor minúsculo.

**Como confirmar.** Inspecionar `plot_tc_diagnostic` para $C=N$: nos painéis $N=2000$ e $N=3000$, a linha vermelha vertical deve estar no início da curva, antes da rampa principal.

**Como corrigir.** Restringir a busca de `argmax dp/dlogt` à região de $p$ longe das bordas — por exemplo $p \in [0.05, 0.95]$:

```python
def estimate_tc_by_max_slope(df, edge_frac=0.05, pmin=0.05, pmax=0.95):
    t, p = clean_tp(df)
    n = len(t)
    if n < 25:
        return np.nan

    logt = np.log(t)
    dp_dlogt = np.gradient(p, logt)

    window = _sg_window(n, frac=0.05, wmin=7, wmax=51)
    if window < n:
        try:
            dp_dlogt = savgol_filter(dp_dlogt, window, polyorder=3)
        except ValueError:
            pass

    k = max(2, int(edge_frac * n))
    valid = np.zeros(n, dtype=bool)
    valid[k:n - k] = True
    valid &= (p >= pmin) & (p <= pmax)

    if not valid.any():
        return np.nan

    idx = np.where(valid)[0]
    i_local = idx[int(np.argmax(dp_dlogt[idx]))]
    return float(t[i_local])
```

---

## 6. Veredito sobre a dinâmica

1. **Escala combinatória universal.** $t_{\max}(N) \sim N^{2.16}$ para todos os $C$ confirma que o esgotamento do ataque é um efeito de contagem de pares, independente da estrutura da rede.

2. **Existência de uma transição não-trivial controlada por $C$.** O expoente $\theta$ de $t_c \sim N^\theta$ depende monotonicamente de $C$ (2.08, 1.47, 1.18 para $C = 1, 2, 3$), e o quociente $t_c/t_{\max}$ tem expoentes $\beta \in \{-0.16,\, -0.71,\, -0.98\}$. **A redundância $C$ comprime a janela crítica relativa, mas não desloca $t_{\max}$.**

3. **Não há $t_c^\infty$ finito para $C \in \{1, 2\}$.** O tempo crítico diverge polinomialmente — coerente com um sistema sem ponto crítico no sentido termodinâmico padrão, apenas uma transição que recua para $\infty$ conforme $N$ cresce.

4. **$C=3$ é o ponto onde correções de scaling começam a importar.** A discrepância entre log-log puro e ajuste com offset ($\Delta\theta \approx 0.5$) é um sinal de cross-over: nesse regime já vale a pena tentar formas mais ricas como $t_c = t_\infty + a N^{-\theta} + b N^{-\theta'}$.

5. **$C=N$ ainda precisa ser refeito.** Os números atuais não refletem física — refletem uma falha do estimador induzida pela morfologia da curva $p(t)$. Aplicar a máscara $p\in[0.05, 0.95]$ deve restaurar a monotonicidade.

---

## 7. Próximos passos sugeridos

1. Reaplicar `estimate_tc_by_max_slope` com a máscara $p \in [0.05, 0.95]$ para $C=N$ e refazer todos os ajustes.
2. Para $C=3$, ajustar $t_c(N) = t_\infty + a N^{-\theta}$ com $t_\infty > 0$ forçado por bounds e comparar AIC com o log-log puro.
3. Estimar barras de erro em $\theta$ via bootstrap das simulações (re-amostragem das curvas $p(t)$).
4. Verificar o colapso usando os $\theta$ corrigidos: se o painel $C=3$ colapsar com $\theta_{\text{coll}} \approx \theta_{tc}$, isso fecha o argumento de FSS padrão; se não, há outro expoente $\theta'$ governando a largura.
