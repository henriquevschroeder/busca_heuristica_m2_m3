# Atividade de Busca Heurística — Distribuição de Tarefas em Máquinas Paralelas

Implementação em C de dois programas de simulação de busca local para o
problema de minimização do *makespan* (tempo de uso de todas as máquinas)
na distribuição de `n` tarefas entre `m` máquinas paralelas.

## Arquivos

| Arquivo                     | Descrição                                              |
|------------------------------|--------------------------------------------------------|
| `busca_monotona.c`           | Busca Local Monótona Randomizada                        |
| `busca_nao_monotona.c`       | Busca Local Não-Monótona (Têmpera Simulada)             |
| `resultados_monotona.csv`    | Saída de exemplo já gerada (540 execuções)              |
| `resultados_nao_monotona.csv`| Saída de exemplo já gerada (300 execuções)              |

## Compilação e execução

```bash
gcc -O2 -o busca_monotona busca_monotona.c -lm
gcc -O2 -o busca_nao_monotona busca_nao_monotona.c -lm

./busca_monotona          # gera resultados_monotona.csv
./busca_nao_monotona      # gera resultados_nao_monotona.csv
```

## Modelagem do problema

- Solução: vetor `assign[i]` = máquina alocada para a tarefa `i`.
- Vizinhança: mover **uma** tarefa para outra máquina (move simples).
- Solução inicial: heurística gulosa **LPT** (Longest Processing Time
  first) — ordena as tarefas da maior para a menor e aloca cada uma na
  máquina menos carregada no momento. Serve apenas como ponto de partida
  razoável para as buscas locais.
- Objetivo: minimizar `makespan = max_j( soma dos p_i alocados em j )`.

## As duas buscas

### 1. Busca Local Monótona Randomizada (`busca_monotona.c`)

Variação sugerida no enunciado: parâmetro `alpha ∈ {0.1,...,0.9}` controla a
frequência de passos de **caminhada aleatória**.

- Com probabilidade `alpha`: move uma tarefa sorteada para uma máquina
  sorteada e aceita o movimento **incondicionalmente** (diversificação).
- Com probabilidade `1-alpha`: só aceita o movimento se ele **reduzir** o
  makespan corrente (busca local tradicional).
- O melhor valor encontrado é guardado à parte e **nunca piora** — é isso
  que garante o caráter monótono do resultado reportado, mesmo a solução de
  trabalho oscilando por causa da caminhada aleatória.
- Parada: 1000 iterações consecutivas sem melhora do melhor valor
  (conforme o enunciado).

### 2. Busca Local Não-Monótona / Têmpera Simulada (`busca_nao_monotona.c`)

Variação sugerida no enunciado: resfriamento geométrico `T = T * alpha`,
com `alpha ∈ {0.80, 0.85, 0.90, 0.95, 0.99}`.

- A cada iteração sorteia-se um vizinho (mover uma tarefa para outra
  máquina).
- Se o vizinho é igual ou melhor, o movimento é sempre aceito.
- Se o vizinho é pior, é aceito com probabilidade `exp(-delta / T)` — é
  isso que torna a busca **não-monótona**: soluções piores podem ser
  aceitas, principalmente quando `T` ainda está alta.
- Parada: 1000 iterações sem melhora do melhor valor **ou** `T` cair
  abaixo de `1e-3` (o que ocorrer primeiro).

## Premissas assumidas (não especificadas no enunciado)

O enunciado (arquivo `enunciadoHeurísticas.pdf`) fica truncado exatamente
na fórmula de resfriamento da Têmpera Simulada e não define, por exemplo, a
temperatura inicial nem o critério de parada específico dela. Foram
assumidos:

- **Temperatura inicial** `T0` = makespan da solução inicial (escala o
  algoritmo automaticamente ao tamanho/dificuldade da instância).
- **Temperatura mínima** `Tmin = 1e-3`.
- **Critério de parada da Têmpera Simulada**: o mesmo das demais buscas
  (1000 iterações sem melhora), complementado por `T < Tmin`.
- **Instância por combinação (m, r)**: gerada **uma única vez** e reutilizada
  em todas as 10 replicações e em todos os valores do parâmetro, para que a
  comparação entre parâmetros seja justa (mesma instância, heurística
  randomizada).
- **Semente aleatória**: `srand(time(NULL))`, uma vez no início do
  programa — cada execução do binário produz resultados diferentes.

## Formato de saída

Ambos os programas gravam um CSV no formato definido no enunciado:

```
heuristica,n,m,replicacao,tempo,iteracoes,valor,parametro
monotona,32,10,1,0.0001,1000,150,0.1
temperasimulada,32,10,1,0.0187,54,167,0.80
```

- `heuristica`: nome da heurística (`monotona` ou `temperasimulada`).
- `n`, `m`: tamanho da instância.
- `replicacao`: número da execução (1 a 10) para aquela instância/parâmetro.
- `tempo`: tempo de CPU gasto na execução, em segundos.
- `iteracoes`: total de iterações realizadas.
- `valor`: melhor makespan encontrado.
- `parametro`: valor de `alpha` usado naquela execução.

## Perguntas do enunciado

Com os dois CSVs gerados, é possível responder diretamente (por exemplo
com Excel/pandas) às perguntas pedidas na atividade:

- Qual heurística demandou mais iterações?
- Qual heurística demandou mais tempo?
- Qual encontrou resultados de maior qualidade (menor `valor`)?
- Quais parâmetros garantem maior qualidade?
- Quais parâmetros são mais rápidos?

Basta agrupar cada CSV por `parametro` (e por `heuristica` ao juntar os
dois arquivos) e comparar médias de `tempo`, `iteracoes` e `valor`.
