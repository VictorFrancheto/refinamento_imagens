# Estimativa de $t_c(N)$ e colapso de dados em ataques aleatórios

Este README descreve o pipeline implementado em [plots_random_attack_refined.ipynb](plots_random_attack_refined.ipynb): o que cada bloco faz, como o tempo pseudo-crítico $t_c(N)$ é extraído das curvas $p(t,N)$ e como o expoente de scaling $\theta$ é obtido por colapso de dados.

Para a derivação matemática completa, ver [metodologia_tc_theta.md](metodologia_tc_theta.md). Para os vereditos físicos sobre os números obtidos, ver [analise_tc.md](analise_tc.md).

---

## 1. Visão geral do pipeline

```
CSVs por (N, C)
      │
      ▼
load_csv_attack + clean_tp           ── leitura, ordenação, filtro de NaN/zeros
      │
      ▼
estimate_tc_by_max_slope             ── t_c(N) para cada curva
estimate_tmax                        ── t_max(N) para cada curva
      │
      ▼
fit_power_law / fit_tc_with_offset   ── scaling t_c ~ N^θ, t_max ~ N^α
      │
      ▼
estimate_theta_by_collapse           ── θ do colapso (min de variância)
      │
      ▼
plot_scaling_tc_tmax                 ── (e), (f) log-log
plot_tc_over_tmax                    ── razão t_c/t_max
test_tmax_scaling                    ── teste assintótico em α
plot_tc_diagnostic                   ── inspeção visual de t_c
plot_collapse_panel                  ── colapso (4 painéis: C=1,2,3,N)
```

A função `plot_final(mode="random")` orquestra tudo. Os tamanhos varridos são $N \in \{300, 500, 600, 1000, 2000, 3000\}$ e os regimes de redundância $C \in \{1, 2, 3, N\}$.

---

## 2. Leitura e limpeza dos dados

**`load_csv_attack(N, C, mode)`**

Carrega `./data_random_medias_random/{mode}_{N}_{C}.csv`. Cada arquivo é a média de várias realizações Monte-Carlo do ataque para o par $(N, C)$. Remove colunas auxiliares `node1`, `node2`; mantém as colunas físicas $t$ e $p$.

**`clean_tp(df)`**

Aplica três operações em sequência:

1. Filtro `np.isfinite(t) & np.isfinite(p) & (t > 0)` — garante que `np.log(t)` seja válido.
2. Ordenação por $t$ crescente — necessária para diferenças centrais e interpolação.
3. Retorna $(t, p)$ como `np.ndarray` de `float`.

---

## 3. Estimativa de $t_c(N)$

**Definição operacional**

$$
t_c(N) \;=\; \arg\max_t \; \frac{dp(t,N)}{d\ln t}.
$$

O ponto onde $p$ varia mais rápido em escala log é o centro da janela crítica. A escolha de $\ln t$ (em vez de $t$) é essencial porque $t$ varre cinco ordens de grandeza nas suas simulações — derivar em $t$ produziria um pico espúrio perto da origem.

**`estimate_tc_by_max_slope(df, edge_frac=0.05)`** faz, em ordem:

1. **Mudança de variável:** `logt = np.log(t)`.
2. **Derivada central em $\ln t$:** `dp_dlogt = np.gradient(p, logt)`. O `np.gradient` aceita eixo não-uniforme e usa diferenças centrais ponderadas pelo espaçamento real.
3. **Suavização Savitzky–Golay grau 3** com janela ímpar adaptativa via `_sg_window(n, frac=0.05, wmin=7, wmax=51)`. SG ajusta localmente um polinômio cúbico em $w$ pontos vizinhos e preserva a posição/amplitude do pico — uma média móvel simples viesaria o argmax em $\sim 2\,g''(t_c)(\Delta \ln t)^2$.
4. **Recorte de borda:** descarta os primeiros e últimos $k = \max(2,\lfloor 0.05 n\rfloor)$ pontos antes de tomar o argmax, eliminando artefatos de extremidade do filtro.
5. **Retorno:** $t_{i^*}$ onde $i^* = \arg\max$ da derivada suavizada no interior.

**`estimate_tmax(df)`** apenas devolve `np.max(t)` — o último tempo observado na simulação.

---

## 4. Ajustes de lei de potência

**`fit_power_law(x, y)`**

Mínimos quadrados em escala log para $y = A\,N^\alpha$:

- ajusta $\ln y = \alpha \ln N + \ln A$ por `np.polyfit`;
- devolve $(\alpha, \ln A, R^2)$.

Usado para os scalings $t_c(N) \sim N^{\theta_{tc}}$, $t_{\max}(N) \sim N^\alpha$ e $t_c/t_{\max} \sim N^\beta$.

**`fit_tc_with_offset(Ns, tcs)`**

Mínimos quadrados não-lineares (Levenberg–Marquardt via `scipy.optimize.curve_fit`) para o modelo de três parâmetros

$$
t_c(N) = t_\infty + a\,N^{-\theta_{tc}}.
$$

Chutes iniciais derivados do ajuste log-log puro. Serve como **diagnóstico**: se $t_\infty$ retorna compatível com zero, não há ponto crítico termodinâmico; se $t_\infty > 0$ com $R^2$ alto, o ajuste log-log puro está enviesado e o expoente real vem deste ajuste não-linear.

---

## 5. Estimativa de $\theta$ por colapso de dados

**Princípio.** A hipótese de finite-size scaling é

$$
p(t,N) \;=\; \mathcal{F}\!\big((t - t_c(N))\,N^{-\theta}\big),
$$

ou seja, todas as curvas $p$ vs $(t-t_c)N^{-\theta}$ devem colapsar sobre uma única função mestra $\mathcal{F}$ quando $\theta$ é o expoente correto.

**`collapse_score_physical(dfs, Ns, tcs, theta, ...)`** mede a qualidade do colapso para um $\theta$ candidato:

1. Para cada curva $k$, calcula $x_k = (t_k - t_c(N_k))\,N_k^{-\theta}$ e aplica máscara $p \in [0.05,\,0.9]$ para focar na janela crítica e ignorar platôs.
2. Constrói grid uniforme $\{x_j\}$ de 400 pontos na **interseção** dos domínios das $K$ curvas e interpola linearmente cada $p_k$ nesse grid, formando a matriz $P_{kj}$.
3. Calcula a curva média $\bar P_j$ e o peso $w_j = |d\bar P/dx|_j + \varepsilon$, que enfatiza onde $\bar P$ varia (janela crítica) e ignora platôs.
4. **Descarte da pior curva** (se $K \ge 4$): elimina a curva com maior desvio quadrático ponderado em relação a $\bar P$, tornando o score robusto a um $t_c$ ruidoso.
5. Retorna $\mathcal{S}(\theta) = \dfrac{\sum_j w_j \cdot \mathrm{Var}_k[P_{kj}]}{\sum_j w_j}$.

**`estimate_theta_by_collapse(dfs, Ns, tcs, ...)`** minimiza $\mathcal{S}(\theta)$ em duas etapas:

1. **Coarse:** avalia $\mathcal{S}$ em 300 pontos uniformes em $\theta \in [0, 6]$; localiza $\theta_1$.
2. **Refine:** avalia $\mathcal{S}$ em 300 pontos numa janela estreita ao redor de $\theta_1$ e devolve $\theta^* = \arg\min \mathcal{S}$.

---

## 6. Diagnósticos e plots

| Função | O que mostra |
|---|---|
| `plot_scaling_tc_tmax` | Painel (e): $t_c(N)$ log-log para os 4 valores de $C$. Painel (f): $t_{\max}(N)$ log-log. Imprime expoentes e R². |
| `plot_tc_over_tmax` | Razão $t_c/t_{\max}$ vs $N$ em log-log. Indica se a janela crítica encolhe (regime FSS) ou se mantém fração constante de $t_{\max}$. |
| `test_tmax_scaling` | Para cada $C$: ajuste direto $t_{\max} \sim N^\alpha$; tabela de $t_{\max}/N^2$; varredura de $\alpha\in[1.8, 2.4]$ minimizando $\mathrm{Var}[\ln(t_{\max}/N^\alpha)]$ — confirma se $\alpha = 2$ exato é compatível. |
| `plot_tc_diagnostic` | Para cada $C$, 6 painéis $p(t)$ em escala log com linha vertical em $t_c(N)$ — inspeção visual de que o estimador acertou a inflexão. |
| `plot_collapse_panel` | Curvas $p$ vs $(t-t_c)N^{-\theta^*}$ com escala `symlog` em $x$ — verificação visual final do colapso. |

---

## 7. Robustez e melhorias em relação à versão original

| Aspecto | Original | Versão refinada |
|---|---|---|
| Suavização da derivada | média móvel fixa de janela 7 | **Savitzky–Golay grau 3** com janela $w = \mathrm{odd}(\mathrm{clip}(0.05n, 7, 51))$ adaptativa em $N$ |
| Bordas | nenhum recorte (`mode="same"` no `convolve` distorce extremidades) | descarte de $\max(2, \lfloor 0.05n\rfloor)$ pontos em cada lado antes do argmax |
| Ajuste de $t_c$ vs $N$ | só log-log puro (assume implicitamente $t_c^\infty = 0$) | log-log puro **+** modelo com offset $t_c = t_\infty + aN^{-\theta}$ via `curve_fit` |
| Score de colapso | sensível a um único $t_c(N)$ ruim | descarte automático da curva mais discordante quando $K \ge 4$ |
| Janelas $[p_{\min}, p_{\max}]$ | inconsistente entre `collapse_score_physical` e `estimate_theta_by_collapse` | unificadas em $[0.05,\,0.9]$ |
| API de colormap | `plt.cm.get_cmap` (depreciado) | `plt.get_cmap` |

---

## 8. Resumo das funções

| Função | Entrada | Saída | Papel no pipeline |
|---|---|---|---|
| `load_csv_attack(N, C, mode)` | $(N, C)$ | DataFrame com colunas $t, p$ | leitura |
| `clean_tp(df)` | DataFrame | arrays $(t, p)$ ordenados | sanitização |
| `_sg_window(n, frac, wmin, wmax)` | $n$ | janela ímpar | helper SG |
| `estimate_tc_by_max_slope(df)` | DataFrame | $t_c(N)$ | estimador principal |
| `estimate_tmax(df)` | DataFrame | $t_{\max}(N)$ | esgotamento |
| `fit_power_law(x, y)` | arrays | $(\alpha, \ln A, R^2)$ | scaling log-log |
| `fit_tc_with_offset(Ns, tcs)` | arrays | $\{t_\infty, a, \theta, R^2\}$ | diagnóstico $t_c^\infty$ |
| `collapse_score_physical(...)` | curvas, $\theta$ | $\mathcal{S}(\theta)$ | qualidade do colapso |
| `estimate_theta_by_collapse(...)` | curvas | $\theta^*$ | expoente do colapso |
| `plot_*` | dicts de curvas | figuras | visualização |
| `print_*` | dicts | stdout | tabelas de scaling |
| `plot_final(mode)` | — | tudo | orquestrador |

---

## 9. Como rodar

1. Garanta os dados em `./data_random_medias_random/random_{N}_{C}.csv` para cada $(N, C)$.
2. Abra [plots_random_attack_refined.ipynb](plots_random_attack_refined.ipynb).
3. Execute as células em ordem (a última chama `plot_final(mode="random")`).
4. Saídas geradas no diretório atual:
   - `scaling_tc_tmax.png`
   - `tc_over_tmax.png`
   - `collapse_corrected.png`
   - tabelas de $t_c(N)$ e expoentes impressas no stdout.

**Dependências:** `numpy`, `pandas`, `matplotlib`, `scipy` (para `savgol_filter` e `curve_fit`). LaTeX no sistema (para `mpl.rcParams["text.usetex"] = True`); se não houver, desligue essa linha.

---

## 10. Limitações conhecidas

- **6 pontos $N$** é o mínimo para distinguir expoentes próximos (e.g. $\alpha = 2.00$ vs $2.10$). Barras de erro por bootstrap das simulações Monte-Carlo seriam o próximo passo natural.
- O estimador de $t_c$ assume **unimodalidade** de $dp/d\ln t$ na janela crítica. Falha quando a curva $p(t)$ tem múltiplos regimes — foi o que ocorreu em $C=N$ para $N \ge 2000$ (ver `analise_tc.md`, §5). Mitigação sugerida: restringir o argmax à região $p \in [0.05, 0.95]$.
- O ajuste com offset tem **3 parâmetros e 6 pontos**: convergência depende de bons chutes iniciais e pode falhar quando os dados não preferem essa forma funcional (e.g. $C=N$, com $R^2 = 0.19$).
