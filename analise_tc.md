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

### O que significa "o estimador falhar"

O estimador de $t_c$ usado aqui é puramente geométrico:

$$
\hat{t}_c(N) \;=\; \arg\max_{t}\; \frac{dp(t,N)}{d\ln t},
$$

ou seja, **localizamos o ponto de inclinação máxima da curva $p(t,N)$ em escala log de tempo**. Esse procedimento só recupera o tempo crítico físico se a curva tiver **um único cotovelo dominante**, correspondente à transição de percolação. "O estimador falha" significa, concretamente, que essa hipótese é violada e o `argmax` passa a apontar para algo que **não é** a transição de interesse. Três modos de falha típicos:

1. **Pico espúrio nas bordas.** Quando $p(t)$ está saturando perto de $0$ ou $1$, ruído numérico ou um transitório inicial podem gerar uma derivada localmente maior que a do verdadeiro cotovelo crítico. O `argmax` escolhe esse pico falso, e $\hat{t}_c$ é "atraído" para a borda da curva. É exatamente o que parece ocorrer em $C=N$ a partir de $N \ge 2000$: $\hat{t}_c$ desaba três ordens de grandeza porque o `argmax` migra do cotovelo de percolação para um pico inicial de saturação local (ver §5).

2. **Curva com dois cotovelos (cross-over).** Se houver duas escalas de tempo competindo — por exemplo, um transitório rápido seguido da percolação propriamente dita — cada uma gera um máximo local em $dp/d\ln t$. Para $N$ pequeno um deles é maior, para $N$ grande o outro, e $\hat{t}_c(N)$ "salta" descontinuamente entre regimes. Esse salto aparece na tabela como $\theta$ negativo e R² baixo, sintomas inequívocos de que o estimador não está medindo uma única escala física.

3. **Janela de suavização inadequada.** Uma janela Savitzky–Golay larga demais alarga e desloca o pico; estreita demais deixa ruído passar e cria máximos espúrios. O código usa janela adaptativa $w = \mathrm{odd}(\mathrm{clip}(0.05\,n,\,7,\,51))$ justamente para minimizar isso, mas em curvas patológicas (como $C=N$ com poucos pontos no platô central) nenhuma escolha de janela resolve sozinha.

Em todos os três casos o sintoma observável é o mesmo: $\hat{t}_c(N)$ deixa de ser uma função monótona e suave de $N$, o R² do ajuste $\ln \hat{t}_c$ vs $\ln N$ cai abruptamente, e o sinal do expoente pode até inverter. **Não é a física que mudou — é o estimador que parou de medir o que se propôs a medir.** O remédio é restringir o domínio do `argmax` à janela onde a transição física certamente vive (máscara $p \in [p_{\min}, p_{\max}]$, ver §5).

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

---

## 8. Conclusão: conexão com *shortest path percolation* e redes de internet quântica

Os números obtidos não são apenas expoentes em ajustes log-log — eles têm uma leitura física direta no contexto da dinâmica de **shortest path percolation (SPP)** que governa a operação de uma rede de internet quântica sob consumo repetido de recursos de emaranhamento.

### 8.1 O mapeamento físico

Em uma rede de internet quântica, cada par de nós $(u,v)$ que precisa estabelecer um canal quântico consome o **caminho mais curto** disponível entre eles — o *entanglement swapping* ao longo desse caminho esgota um "estoque" de pares de Bell (estados maximamente emaranhados de dois qubits, historicamente chamados de pares EPR em referência ao artigo de Einstein–Podolsky–Rosen de 1935) em cada aresta intermediária. Isso define a dinâmica de SPP:

- a cada passo $t$, sorteia-se uma demanda $(u,v)$;
- o caminho mais curto $\pi_t(u,v)$ é computado na rede residual $G_t$;
- as $|\pi_t|$ arestas ao longo dele têm sua **capacidade** $C$ decrementada em 1;
- arestas com capacidade exaurida são removidas, $G_{t+1} = G_t \setminus \{e : C_e = 0\}$.

O parâmetro $C$ — que na nossa varredura assume valores $\{1, 2, 3, N\}$ — é exatamente a **redundância de pares de Bell por enlace**, ou equivalentemente a taxa de geração de emaranhamento por unidade de tempo de uso. O observável $p(t,N) = $ fração de pares ainda conectados na componente gigante é o **rendimento operacional** da rede.

### 8.2 Leitura física dos expoentes

| Resultado numérico | Significado em SPP | Implicação para internet quântica |
|---|---|---|
| $t_{\max} \sim N^{2.16}$ universal | Esgotamento combinatório: $\binom{N}{2}$ pares possíveis | Tempo de vida absoluto da rede é ditado pela contagem de demandas, não pela topologia |
| $\theta_{C=1} = 2.08 \approx \alpha$ | $t_c/t_{\max} \to $ const. | **Sem proteção:** rede colapsa só perto do esgotamento; não há aviso prévio |
| $\theta_{C=2} = 1.47$, $\beta = -0.71$ | Janela crítica sublinear em $t_{\max}$ | Redundância dupla cria escala crítica genuína $\ll t_{\max}$ |
| $\theta_{C=3} = 1.18$, $\beta \approx -1$ | $t_c/t_{\max} \sim 1/N$ | Transição abrupta com janela vanishing — comportamento de **percolação clássica** |
| Monotonicidade $\theta_C \downarrow$ com $C$ | Mais redundância $\Rightarrow$ transição mais precoce em fração de $t_{\max}$ | Trade-off: aumentar $C$ baixa o $t_c$ relativo mas estreita a janela de degradação |

### 8.3 Por que $\theta$ decresce com $C$ — interpretação de redes complexas

A intuição ingênua diria o oposto: mais Bell pairs por aresta deveriam **adiar** o colapso. O fato observado — $\theta_C$ decresce — revela que estamos medindo a **largura relativa** da transição, não sua posição absoluta. Em termos de teoria de percolação:

- Com $C=1$, cada aresta é frágil; o sistema se comporta como **percolação de arestas independentes**, onde a remoção é local e descorrelacionada. A correlação $\xi(t)$ diverge lentamente porque o dano se espalha por difusão geométrica. Isso produz uma janela crítica larga, com $\xi \sim N$ (corte por tamanho finito), e portanto $\theta \to d_{\text{eff}}$.

- Com $C \ge 2$, a dinâmica adquire **memória de caminhos**: arestas em caminhos centrais (alto *edge betweenness*) são reutilizadas sistematicamente, e o consumo deixa de ser aleatório sobre arestas — torna-se **correlacionado pela topologia**. As arestas centrais saturam primeiro, e a rede passa por uma **reconfiguração dos caminhos mais curtos** (o que em redes complexas se conhece como *rerouting cascade*). Essa cascata é o que estreita a janela crítica: uma vez que o backbone de baixo-diâmetro quebra, o sistema desmorona em $\Delta t = o(t_c)$.

- O limite $C=3$ com $\beta \approx -1$ corresponde precisamente ao regime onde a janela crítica escala como o **diâmetro inverso** de uma rede de mundo pequeno, $t_c/t_{\max} \sim 1/\log N$ ou $\sim 1/N$ dependendo da topologia subjacente. Esse é o **regime de percolação explosiva controlada**, característico de SPP em grafos com distribuição de betweenness heavy-tailed.

### 8.4 Implicações para o projeto de redes de internet quântica

1. **Capacidade não é a métrica que importa — é a curvatura de $p(t)$.** O expoente $\alpha \approx 2.16$ universal mostra que duplicar $C$ não duplica o tempo de vida operacional; o ganho está em $\theta$, ou seja, em **quão abrupta** é a degradação. Redes com $C$ alto são *quase indistinguíveis* de redes saudáveis até $t \sim t_c$, e então colapsam em uma janela $\sim N^{\theta-\alpha}$ — o que do ponto de vista operacional significa **ausência de degradação gradual**: o operador não tem sinais precoces de fadiga.

2. **A redundância $C$ é um trade-off entre robustez e detectabilidade.** $C=1$ é frágil mas previsível (a degradação acompanha o uso); $C=3$ é robusto mas perigoso (sem aviso). Em internet quântica, isso impõe que **políticas de monitoramento** sejam adaptadas ao $C$ provisionado: para $C$ alto, é obrigatório monitorar **derivadas** de $p(t)$, não níveis.

3. **A universalidade de $\alpha$ atravessa a topologia.** Como $t_{\max} \sim N^{2.16}$ não depende de $C$ — e, em particular, está acima do limite trivial $\binom{N}{2} \sim N^2$ por um fator logarítmico-em-$N$ disfarçado — temos evidência de que o **comprimento médio de caminho** $\langle \ell \rangle$ entra na contagem: cada demanda consome $\langle \ell \rangle$ unidades de capacidade. Para redes de mundo pequeno, $\langle \ell \rangle \sim \log N$, o que reproduz $\alpha \approx 2 + \epsilon$ observado.

4. **O caso $C=N$ é o limite de "capacidade infinita" — e é onde mora a física assintótica do SPP.** O fato de o estimador falhar nesse regime (§5) é em si um achado: significa que a curva $p(t,N)$ desenvolve **dois cotovelos** — um transitório topológico (perda das arestas de menor betweenness, sem afetar a conectividade) e a percolação genuína. Em redes de internet quântica com capacidade essencialmente ilimitada, o operador veria **duas escalas de tempo distintas**: uma de "envelhecimento topológico" e outra de "morte por percolação", separadas por um platô. Identificar e separar essas duas escalas — não meramente "corrigir o estimador" — é o próximo problema teórico relevante.

### 8.5 Síntese

Os resultados confirmam que a transição observada em $p(t,N)$ é uma **transição de percolação de tamanho finito sobre o ensemble de caminhos mais curtos**, não sobre o ensemble de arestas. O expoente $\theta$ não mede uma dimensão geométrica do grafo — mede como a **distribuição de carga sobre esses caminhos** se concentra à medida que arestas são removidas. É essa concentração que governa a janela operacional de uma internet quântica real: redes com alta redundância são objetivamente mais duradouras, porém estatisticamente mais frágeis no sentido de que sua falha é **menos anunciada**. Em última análise, o que a tabela do §1 está nos dizendo é que **a métrica de projeto correta não é $t_c$ nem $t_{\max}$ isoladamente, mas o expoente $\beta = \theta - \alpha$** — a única quantidade que captura simultaneamente quanto a rede dura *e* quanto tempo o operador tem para reagir.
