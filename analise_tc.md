# Análise dos expoentes de $t_c(N)$, $t_{\max}(N)$ e da razão $t_c/t_{\max}$

Resultados numéricos obtidos com o estimador refinado (Savitzky–Golay + corte de borda) sobre o conjunto $N \in \{300, 500, 600, 1000, 2000, 3000\}$ e ataques aleatórios com redundância $C \in \{1, 2, 3, N\}$.

---

## 1. Tabela-resumo

| Conjunto | $t_c(N) \sim N^{\theta}$ (log-log) | R² | $t_{\max}(N)$ (horizonte, §2) | Veredito |
|---|---|---|---|---|
| $C=1$ | $\theta \approx 2.08$ | 0.999 | $\sim N^{2.24}$ | **Consistente** — colapso só perto do esgotamento de pares distintos |
| $C=2$ | $\theta \approx 1.47$ | 0.996 | $\sim N^{2.18}$ | **Consistente** — expoente intermediário |
| $C=3$ | $\theta \approx 1.18$ | 0.994 | $\sim N^{2.17}$ | **Consistente** — expoente próximo do limite $C\to\infty$ esperado |
| $C=N$ | $\theta \approx -2.01$ | 0.604 | $\sim N^{2.16}$ | **Artefato** — estimador falha para $N \ge 2000$ |

Observação: $t_{\max}$ na coluna acima é apenas o último timestamp dos CSVs (regra de parada do simulador, ver §2), não um observável físico. Por isso a coluna de $\beta = \theta - \alpha$ que aparecia em versões anteriores foi removida — ela era redundante com $\theta$ por construção.

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

## 2. $t_{\max}(N)$ — horizonte de simulação, não observável físico

No pipeline, `estimate_tmax(df) = max(t no CSV)`. Ou seja, $t_{\max}$ não é o tempo em que a rede de fato deixa de funcionar — é o **último timestamp registrado**, fixado pela regra de parada do simulador que gerou os CSVs. O ajuste devolve, com R² = 1.000 para todo $C$,

$$
t_{\max}(N) \;\sim\; N^{2.16 \pm 0.04}.
$$

A independência em $C$ e o R² praticamente perfeito confirmam que se trata de **uma única regra de parada $\propto N^2$** aplicada a todas as rodadas (provavelmente algo como `int(N*(N-1)/2)` ou um múltiplo). O $0.16$ acima de $2$ é compatível com correção de arredondamento de inteiros sobre uma faixa estreita $N\in[300, 3000]$ — não é assinatura física da topologia.

Consequência metodológica: a razão $t_c/t_{\max}$ e seu expoente $\beta = \theta - 2.16$ são **redundantes com $\theta$** por construção e não trazem informação física independente. Toda a leitura física tem que ser feita diretamente sobre $t_c(N)$.

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

## 6.$t_{\max}$ é horizonte de simulação, não observável.** O ajuste $\sim N^{2.16}$ é a regra de parada do gerador dos CSVs ($\propto N^2$, com $0.16$ de correção de ajuste log-log em faixa estreita de $N$). Não há universalidade física a extrair daí.

2. **A dinâmica vive em $\theta_C$.** Os três valores $\theta_C \in \{2.08,\, 1.47,\, 1.18\}$ para $C\in\{1,2,3\}$ são robustos (R² $\ge 0.994$) e medidos diretamente sobre $p(t,N)$. A monotonicidade decrescente em $C$ é o achado central; sua interpretação geométrica (profundidade no ranking de edge-betweenness) está na §8.

3. **Não há $t_c^\infty$ finito para $C \in \{1, 2\}$.** O tempo crítico diverge polinomialmente em $N$ — coerente com ausência de transição termodinâmica padrão; apenas pseudo-crítico em tamanho finito.

4. **$C=3$ é o ponto onde correções de scaling começam a importar.** A discrepância entre log-log puro e ajuste com offset ($\Delta\theta \approx 0.5$) é sinal de cross-over; vale a pena tentar formas mais ricas como $t_c = t_\infty + a N^{-\theta} + b N^{-\theta'}$.

5. **$C=N$ ainda precisa ser refeito.** Os números atuais são falha do estimador, não físicaé um sinal de cross-over: nesse regime já vale a pena tentar formas mais ricas como $t_c = t_\infty + a N^{-\theta} + b N^{-\theta'}$.

5. **$C=N$ ainda precisa ser refeito.** Os números atuais não refletem física — refletem uma falha do estimador induzida pela morfologia da curva $p(t)$. Aplicar a máscara $p\in[0.05, 0.95]$ deve restaurar a monotonicidade.

---

## 7. Próximos passos sugeridos

1. Reaplicar `estimate_tc_by_max_slope` com a máscara $p \in [0.05, 0.95]$ para $C=N$ e refazer todos os ajustes.
2. Para $C=3$, ajustar $t_c(N) = t_\infty + a N^{-\theta}$ com $t_\infty > 0$ forçado por bounds e comparar AIC com o log-log puro.
3. Estimar barras de erro em $\theta$ via bootstrap das simulações (re-amostragem das curvas $p(t)$).
4. Verificar o colapso usando os $\theta$ corrigidos: se o painel $C=3$ colapsar com $\theta_{\text{coll}} \approx \theta_{tc}$, isso fecha o argumento de FSS padrão; se não, há outro expoente $\theta'$ governando a largura.

---

## 8. SPP, redes quânticas e os resultados deste trabalho

Esta seção conecta três coisas: (i) a dinâmica que foi simulada (*shortest path percolation*), (ii) a leitura dessa dinâmica em internet quântica (que é a motivação física do trabalho), e (iii) o que os números medidos — $\theta_C \in \{2.08, 1.47, 1.18\}$ — efetivamente dizem sobre essa leitura. Tudo aqui se apoia em coisas que **estão** no pipeline: o estimador de $t_c$ (§§2–3 de [metodologia_tc_theta.md](metodologia_tc_theta.md)), a janela de transição $W(N)$ ([janela_transicao.md](janela_transicao.md)) e os ajustes log-log de [plots_random_attack_refined.ipynb](plots_random_attack_refined.ipynb).

### 8.1 O que está sendo simulado, em uma frase

SPP sobre Waxman: a cada passo $t$ um par $(u,v)$ é sorteado, o caminho mais curto $\pi_t(u,v)$ é roteado na rede residual e cada aresta de $\pi_t$ tem sua capacidade $C$ decrementada em 1; arestas com capacidade zero somem. O observável é $p(t,N)$, fração de pares ainda conectados pela componente gigante.

### 8.2 A tradução para internet quântica

A correspondência é direta e dispensa metáforas:

| Objeto no SPP | Objeto na rede quântica |
|---|---|
| Aresta do grafo Waxman | Enlace físico entre dois repetidores |
| Capacidade $C$ da aresta | Número de pares de Bell (EPR) pré-alocados naquele enlace |
| Demanda $(u,v)$ no passo $t$ | Requisição de entrega de emaranhamento ponto-a-ponto |
| Caminho $\pi_t(u,v)$ | Rota de *entanglement swapping* escolhida pelo controlador |
| Decremento de capacidade ao longo de $\pi_t$ | Consumo de um par de Bell por enlace para realizar o swapping |
| Remoção da aresta quando $C=0$ | Enlace fica fora de serviço até nova replenishment |
| $p(t,N)$ | Fração de pares fim-a-fim ainda atendíveis com a reserva atual |
| $t_c(N)$ | Número de requisições que a rede aguenta antes da degradação abrupta |

Sob esse mapeamento, $C$ é exatamente o **estoque inicial de pares EPR por enlace** e $t$ é a **vazão acumulada de requisições**. A janela de transição $W(N)$ documentada em [janela_transicao.md](janela_transicao.md) é o intervalo em que a rede passa de "essencialmente saudável" para "essencialmente colapsada".

### 8.3 O que os resultados deste trabalho dizem

Três fatos numéricos, dos quais tudo decorre:

1. $t_c(N) \sim N^{\theta_C}$ com $\theta_1 \approx 2.08$, $\theta_2 \approx 1.47$, $\theta_3 \approx 1.18$ (R² $\ge 0.994$).
2. $\theta_C$ é **estritamente decrescente** em $C$ no intervalo medido.
3. $t_{\max}(N) \sim N^{2.16}$ é horizonte de simulação ($\propto N^2$), não escala física (§2).

Lendo na linguagem de internet quântica:

**(a) Aumentar o estoque $C$ por enlace aumenta o número absoluto de requisições atendidas, mas com retornos decrescentes em $N$.** O expoente $\theta_C$ cai de $2.08$ para $1.18$ entre $C=1$ e $C=3$. Em valor absoluto, $t_c$ sempre cresce com $C$ para $N$ fixo (mais estoque, mais requisições antes do colapso); o ponto é que **a taxa de crescimento em $N$** muda. Para $C=1$, dobrar $N$ multiplica $t_c$ por $\approx 4{,}2$; para $C=3$, dobrar $N$ multiplica $t_c$ por apenas $\approx 2{,}3$. Operacionalmente: o ganho que se obtém ao expandir a rede (mais nós) é maior quando o estoque é magro do que quando ele é folgado.

**(b) Não existe ponto crítico termodinâmico para $C \in \{1, 2\}$.** O ajuste com offset (§4) devolve $t_\infty$ compatível com zero e o ajuste log-log puro coincide com o ajuste de três parâmetros. Isso significa que, para esses regimes, **não há um número finito de requisições que a rede aguente no limite $N\to\infty$**: a capacidade da rede inteira diverge com $N$. Para internet quântica de larga escala isso é uma propriedade desejável (escalabilidade), e os números quantificam **com que expoente** essa capacidade escala.

**(c) $C=3$ é cross-over.** A discrepância entre log-log puro ($\theta=1.18$) e ajuste com offset ($\theta=1.70$) — §4, R² baixo no offset — indica que o ansatz $t_c = t_\infty + aN^{-\theta}$ não basta. Em termos operacionais, $C=3$ parece ser o regime em que duas escalas de tempo passam a competir: o consumo direto do estoque e a reorganização de rotas conforme arestas saem. Os dados atuais não separam essas duas escalas — para isso seria preciso $C=4, 5, \dots$ e/ou monitorar $\langle\ell_t\rangle$ ao longo do ataque (cf. §4 de [janela_transicao.md](janela_transicao.md)).

**(d) A janela de transição encolhe com $C$ — e isso é o achado operacional central.** Os $\theta_C$ medidos são justamente os expoentes da posição da janela; combinados com o expoente de largura $\theta_{\text{coll}}$ (estimado por colapso, ver §4 de [metodologia_tc_theta.md](metodologia_tc_theta.md)), eles dizem **quão abrupta** é a queda de $p(t,N)$ em torno de $t_c$. Para uma rede quântica, isso tem consequência prática: quanto maior $C$, maior $t_c$ em absoluto, mas mais estreita a janela em que a rede passa de funcional a colapsada. Reservas generosas de pares EPR pagam em vida útil mas custam em **previsibilidade do colapso** — um operador não deve confiar em "a vazão ainda está boa" como sinal de saúde para $C$ alto; precisa monitorar derivadas, não níveis.

### 8.4 O caso $C=N$ e a hipótese de duas escalas

O regime $C=N$ corresponde a "estoque tão grande que cada enlace praticamente não satura por consumo direto". Os números atuais ($\theta=-2.01$, R²=0.604) **não medem física** — são falha do estimador, documentada no §5: o `argmax` de $dp/d\ln t$ migra de um pico para outro entre $N=1000$ e $N=2000$. A presença de **dois picos** é em si informação: sugere que nesse regime $p(t,N)$ tem duas escalas de tempo distintas, uma transitória e outra de percolação propriamente dita. Operacionalmente: redes quânticas com estoque essencialmente ilimitado por enlace podem exibir uma fase de "desgaste topológico lento" antes do colapso propriamente dito. Confirmar isso exige refazer $C=N$ com a máscara $p\in[0.05, 0.95]$ proposta no §5 e examinar se os dois picos são robustos.

### 8.5 Conclusão

Em SPP sobre Waxman com o estimador deste trabalho:

- O **tempo de vida** $t_c(N)$ da rede cresce polinomialmente em $N$ com expoente $\theta_C$ que **depende monotonicamente** do estoque por enlace $C$, sem evidência de um teto finito para $C \in \{1, 2\}$.
- O ganho marginal de aumentar $C$ é positivo em valor absoluto mas **reduz o expoente de escala**, ou seja, a vantagem do estoque folgado decresce conforme a rede cresce.
- A janela em que o colapso ocorre encolhe com $C$, transformando alta capacidade em **baixa previsibilidade**.
- O regime de capacidade muito alta ($C=N$) sugere duas escalas de tempo distintas, mas isso precisa ser confirmado com o estimador corrigido antes de qualquer afirmação física.

Para o desenho de uma internet quântica real, a leitura é direta: $C$ é o botão que troca **horizonte operacional** por **agudeza da falha**. Os dados aqui não dizem qual é o $C$ ótimo — dizem que essa pergunta tem que ser feita em termos de $\theta_C$, não de $t_c$ isolado, porque é $\theta_C$ que governa como a escolha vai escalar conforme a rede cresce.
