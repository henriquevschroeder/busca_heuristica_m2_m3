/* ============================================================================
 * busca_nao_monotona.c
 *
 * Universidade do Vale do Itajaí - Complexidade de Algoritmos
 * Atividade de Busca Heurística
 *
 * BUSCA LOCAL NÃO-MONÓTONA (TÊMPERA SIMULADA / SIMULATED ANNEALING)
 * ---------------------------------------------------------------------------
 * Problema: distribuição de n tarefas entre m máquinas paralelas, onde a
 * tarefa i possui tempo de processamento p_i. Objetivo: minimizar o tempo de
 * uso de todas as máquinas (makespan = maior carga entre as máquinas).
 *
 * Variação usada (sugerida no enunciado):
 *   A temperatura varia entre iterações adjacentes segundo t = t * alpha,
 *   com alpha ∈ {0.80, 0.85, 0.90, 0.95, 0.99}.
 *
 *   A cada iteração é sorteado um vizinho (move-se uma tarefa sorteada para
 *   uma máquina de destino sorteada, diferente da atual):
 *     - se o vizinho melhora (ou empata) o makespan corrente, o movimento é
 *       SEMPRE aceito;
 *     - caso contrário (o vizinho piora a solução), o movimento é aceito com
 *       probabilidade exp(-delta / T) -- é isso que torna a busca
 *       NÃO-MONÓTONA: soluções piores podem ser aceitas, especialmente
 *       quando a temperatura T ainda está alta, permitindo escapar de
 *       ótimos locais.
 *   A cada iteração a temperatura é resfriada geometricamente: T = T * alpha.
 *
 * ATENÇÃO - premissas assumidas (o enunciado não define esses detalhes para
 * a Têmpera Simulada, apenas a fórmula de resfriamento):
 *   - Temperatura inicial T0: makespan da solução inicial (proporcional à
 *     escala do problema, para que, no início, praticamente qualquer
 *     movimento tenha boa chance de ser aceito).
 *   - Temperatura mínima Tmin = 1e-3 (abaixo disso, T é tratada como muito
 *     fria e a busca para).
 *   - Critério de parada: 1000 iterações consecutivas sem melhora do melhor
 *     valor encontrado (mesmo critério usado na busca monótona) OU
 *     T < Tmin, o que ocorrer primeiro.
 *
 * Instâncias simuladas: m ∈ {10, 20, 50} máquinas; n = m^r tarefas, com
 * r ∈ {1.5, 2.0}; tempo de cada tarefa sorteado uniformemente em [1,100].
 * Cada instância (m, r) é executada 10 vezes (replicações) para cada valor
 * de alpha.
 *
 * Saída: resultados_nao_monotona.csv, no formato
 *   heuristica,n,m,replicacao,tempo,iteracoes,valor,parametro
 *
 * Compilar: gcc -O2 -o busca_nao_monotona busca_nao_monotona.c -lm
 * Executar: ./busca_nao_monotona
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MAX_ITER_SEM_MELHORA 1000
#define NUM_REPLICACOES 10
#define T_MIN 1e-3

/* ---------------------------------------------------------------------- */
/* Estrutura da solução (igual a busca_monotona.c)                         */
/* ---------------------------------------------------------------------- */
typedef struct {
    int n;
    int m;
    int  *p;
    int  *assign;
    long *load;
} Solucao;

/* ---------------------------------------------------------------------- */
/* Utilidades de numeros aleatorios                                        */
/* ---------------------------------------------------------------------- */
static double rand01(void) {
    return (double) rand() / ((double) RAND_MAX + 1.0);
}

static int rand_int(int loInclusive, int hiExclusive) {
    return loInclusive + rand() % (hiExclusive - loInclusive);
}

/* ---------------------------------------------------------------------- */
/* Funcoes basicas sobre a solucao                                         */
/* ---------------------------------------------------------------------- */
static long makespan(const Solucao *s) {
    long mx = 0;
    for (int j = 0; j < s->m; j++)
        if (s->load[j] > mx) mx = s->load[j];
    return mx;
}

static Solucao *criar_solucao(int n, int m) {
    Solucao *s = malloc(sizeof(Solucao));
    s->n = n;
    s->m = m;
    s->p = malloc(sizeof(int) * n);
    s->assign = malloc(sizeof(int) * n);
    s->load = calloc(m, sizeof(long));
    return s;
}

static void liberar_solucao(Solucao *s) {
    free(s->p);
    free(s->assign);
    free(s->load);
    free(s);
}

static void gerar_tarefas(Solucao *s) {
    for (int i = 0; i < s->n; i++)
        s->p[i] = 1 + rand() % 100;
}

typedef struct { int idx; int val; } Par;

static int cmp_par_desc(const void *a, const void *b) {
    return ((const Par *) b)->val - ((const Par *) a)->val;
}

static void solucao_inicial_gulosa(Solucao *s) {
    Par *ord = malloc(sizeof(Par) * s->n);
    for (int i = 0; i < s->n; i++) { ord[i].idx = i; ord[i].val = s->p[i]; }
    qsort(ord, s->n, sizeof(Par), cmp_par_desc);

    for (int j = 0; j < s->m; j++) s->load[j] = 0;

    for (int k = 0; k < s->n; k++) {
        int i = ord[k].idx;
        int melhor = 0;
        for (int j = 1; j < s->m; j++)
            if (s->load[j] < s->load[melhor]) melhor = j;
        s->assign[i] = melhor;
        s->load[melhor] += s->p[i];
    }
    free(ord);
}

static void aplicar_movimento(Solucao *s, int tarefa, int destino) {
    int origem = s->assign[tarefa];
    s->load[origem] -= s->p[tarefa];
    s->load[destino] += s->p[tarefa];
    s->assign[tarefa] = destino;
}

/* ---------------------------------------------------------------------- */
/* Têmpera Simulada (busca local não-monótona)                             */
/* ---------------------------------------------------------------------- */
typedef struct {
    long valor;
    long iteracoes;
} ResultadoBusca;

static ResultadoBusca tempera_simulada(Solucao *s, double alpha) {
    long valorCorrente = makespan(s);
    long melhorValor = valorCorrente;

    double T = (double) valorCorrente > 0 ? (double) valorCorrente : 1.0;

    long iterSemMelhora = 0;
    long totalIter = 0;

    while (iterSemMelhora < MAX_ITER_SEM_MELHORA && T > T_MIN) {
        totalIter++;

        int tarefa = rand_int(0, s->n);
        int destino;
        do {
            destino = rand_int(0, s->m);
        } while (destino == s->assign[tarefa]);

        int origem = s->assign[tarefa];
        long loadOrigemAntiga = s->load[origem];
        long loadDestinoAntiga = s->load[destino];

        aplicar_movimento(s, tarefa, destino);
        long valorVizinho = makespan(s);
        long delta = valorVizinho - valorCorrente;

        int aceitar;
        if (delta <= 0) {
            aceitar = 1; /* vizinho igual ou melhor: sempre aceita */
        } else {
            /* vizinho pior: aceita com probabilidade exp(-delta / T)
             * (isso torna a busca nao-monotona) */
            double prob = exp(-(double) delta / T);
            aceitar = (rand01() < prob);
        }

        if (aceitar) {
            valorCorrente = valorVizinho;
        } else {
            /* desfaz o movimento */
            s->load[origem] = loadOrigemAntiga;
            s->load[destino] = loadDestinoAntiga;
            s->assign[tarefa] = origem;
        }

        if (valorCorrente < melhorValor) {
            melhorValor = valorCorrente;
            iterSemMelhora = 0;
        } else {
            iterSemMelhora++;
        }

        /* resfriamento geometrico: t = t * alpha */
        T = T * alpha;
    }

    ResultadoBusca r;
    r.valor = melhorValor;
    r.iteracoes = totalIter;
    return r;
}

/* ---------------------------------------------------------------------- */
/* Programa principal: gera as instancias e executa os experimentos        */
/* ---------------------------------------------------------------------- */
int main(void) {
    srand((unsigned int) time(NULL));

    int maquinas[] = {10, 20, 50};
    double rs[] = {1.5, 2.0};
    double alphas[] = {0.80, 0.85, 0.90, 0.95, 0.99};

    int numM = sizeof(maquinas) / sizeof(maquinas[0]);
    int numR = sizeof(rs) / sizeof(rs[0]);
    int numAlpha = sizeof(alphas) / sizeof(alphas[0]);

    FILE *out = fopen("resultados_nao_monotona.csv", "w");
    if (!out) {
        perror("erro ao criar resultados_nao_monotona.csv");
        return 1;
    }
    fprintf(out, "heuristica,n,m,replicacao,tempo,iteracoes,valor,parametro\n");

    for (int im = 0; im < numM; im++) {
        int m = maquinas[im];
        for (int ir = 0; ir < numR; ir++) {
            double r = rs[ir];
            int n = (int) llround(pow((double) m, r));

            Solucao *base = criar_solucao(n, m);
            gerar_tarefas(base);

            printf("Instancia: m=%d, r=%.1f -> n=%d\n", m, r, n);

            for (int ia = 0; ia < numAlpha; ia++) {
                double alpha = alphas[ia];

                for (int rep = 1; rep <= NUM_REPLICACOES; rep++) {
                    Solucao *s = criar_solucao(n, m);
                    memcpy(s->p, base->p, sizeof(int) * n);
                    solucao_inicial_gulosa(s);

                    clock_t inicio = clock();
                    ResultadoBusca res = tempera_simulada(s, alpha);
                    clock_t fim = clock();
                    double tempo = (double) (fim - inicio) / CLOCKS_PER_SEC;

                    fprintf(out, "temperasimulada,%d,%d,%d,%.4f,%ld,%ld,%.2f\n",
                            n, m, rep, tempo, res.iteracoes, res.valor, alpha);

                    liberar_solucao(s);
                }
            }
            liberar_solucao(base);
        }
    }

    fclose(out);
    printf("Concluido. Resultados gravados em resultados_nao_monotona.csv\n");
    return 0;
}
