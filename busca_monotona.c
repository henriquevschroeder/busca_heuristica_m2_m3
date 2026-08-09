/* ============================================================================
 * busca_monotona.c
 *
 * Universidade do Vale do Itajaí - Complexidade de Algoritmos
 * Atividade de Busca Heurística
 * 
 * **Alunos:**
 * - Henrique Schroeder
 * - Nilson Andrade
 * - Matheus Barbiéri
 * - Vinícius Grisa
 *
 * BUSCA LOCAL MONÓTONA RANDOMIZADA
 * ---------------------------------------------------------------------------
 * Problema: distribuição de n tarefas entre m máquinas paralelas, onde a
 * tarefa i possui tempo de processamento p_i. Objetivo: minimizar o tempo de
 * uso de todas as máquinas (makespan = maior carga entre as máquinas).
 *
 * Variação usada (sugerida no enunciado):
 *   Parâmetro alpha ∈ {0.1,...,0.9} indica a frequência com que é realizada
 *   uma caminhada aleatória no espaço de busca.
 *
 *   A cada iteração:
 *     - com probabilidade alpha: dá um passo de CAMINHADA ALEATÓRIA, ou
 *       seja, move uma tarefa sorteada para uma máquina sorteada e aceita o
 *       movimento incondicionalmente (serve para diversificar a busca e
 *       escapar de ótimos locais);
 *     - com probabilidade (1 - alpha): dá um passo de BUSCA LOCAL comum
 *       (first-improvement), só aceitando o movimento se ele reduzir o
 *       makespan corrente.
 *
 *   A busca é "monótona" no sentido de que o MELHOR valor encontrado
 *   (registrado separadamente da solução de trabalho) nunca piora ao longo
 *   da execução, mesmo quando passos de caminhada aleatória pioram a
 *   solução corrente.
 *
 * Critério de parada: 1000 iterações consecutivas sem melhora do melhor
 * valor encontrado (conforme enunciado).
 *
 * Instâncias simuladas: m ∈ {10, 20, 50} máquinas; n = m^r tarefas, com
 * r ∈ {1.5, 2.0}; tempo de cada tarefa sorteado uniformemente em [1,100].
 * Cada instância (m, r) é executada 10 vezes (replicações) para cada valor
 * de alpha.
 *
 * Saída: resultados_monotona.csv, no formato
 *   heuristica,n,m,replicacao,tempo,iteracoes,valor,parametro
 *
 * Compilar: gcc -O2 -o busca_monotona busca_monotona.c -lm
 * Executar: ./busca_monotona
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MAX_ITER_SEM_MELHORA 1000
#define NUM_REPLICACOES 10

/* ---------------------------------------------------------------------- */
/* Estrutura da solução                                                    */
/* ---------------------------------------------------------------------- */
typedef struct {
    int n;          /* numero de tarefas */
    int m;          /* numero de maquinas */
    int  *p;        /* p[i]      = tempo de processamento da tarefa i */
    int  *assign;   /* assign[i] = maquina alocada para a tarefa i */
    long *load;     /* load[j]   = soma dos p_i das tarefas na maquina j */
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

/* gera os tempos de processamento aleatorios das tarefas, entre 1 e 100 */
static void gerar_tarefas(Solucao *s) {
    for (int i = 0; i < s->n; i++)
        s->p[i] = 1 + rand() % 100;
}

typedef struct { int idx; int val; } Par;

static int cmp_par_desc(const void *a, const void *b) {
    return ((const Par *) b)->val - ((const Par *) a)->val;
}

/* solucao inicial gulosa (LPT - Longest Processing Time first): ordena as
 * tarefas da maior para a menor e aloca cada uma na maquina menos carregada
 * no momento. Usada como ponto de partida para a busca local. */
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

/* aplica o movimento: move a tarefa `tarefa` para a maquina `destino` */
static void aplicar_movimento(Solucao *s, int tarefa, int destino) {
    int origem = s->assign[tarefa];
    s->load[origem] -= s->p[tarefa];
    s->load[destino] += s->p[tarefa];
    s->assign[tarefa] = destino;
}

/* ---------------------------------------------------------------------- */
/* Busca Local Monótona Randomizada                                        */
/* ---------------------------------------------------------------------- */
typedef struct {
    long valor;
    long iteracoes;
} ResultadoBusca;

static ResultadoBusca busca_monotona_randomizada(Solucao *s, double alpha) {
    long melhorValor = makespan(s);
    long iterSemMelhora = 0;
    long totalIter = 0;

    while (iterSemMelhora < MAX_ITER_SEM_MELHORA) {
        totalIter++;

        int tarefa = rand_int(0, s->n);
        int destino;
        do {
            destino = rand_int(0, s->m);
        } while (destino == s->assign[tarefa]);

        if (rand01() < alpha) {
            /* passo de caminhada aleatoria: aceita incondicionalmente */
            aplicar_movimento(s, tarefa, destino);
        } else {
            /* busca local monotona: so aceita se melhora o makespan corrente */
            int origem = s->assign[tarefa];
            long loadOrigemAntiga = s->load[origem];
            long loadDestinoAntiga = s->load[destino];
            long makespanAntigo = makespan(s);

            aplicar_movimento(s, tarefa, destino);
            long makespanNovo = makespan(s);

            if (makespanNovo >= makespanAntigo) {
                /* nao melhorou: desfaz o movimento */
                s->load[origem] = loadOrigemAntiga;
                s->load[destino] = loadDestinoAntiga;
                s->assign[tarefa] = origem;
            }
        }

        long valorCorrente = makespan(s);
        if (valorCorrente < melhorValor) {
            melhorValor = valorCorrente;
            iterSemMelhora = 0;
        } else {
            iterSemMelhora++;
        }
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
    double alphas[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    int numM = sizeof(maquinas) / sizeof(maquinas[0]);
    int numR = sizeof(rs) / sizeof(rs[0]);
    int numAlpha = sizeof(alphas) / sizeof(alphas[0]);

    FILE *out = fopen("resultados_monotona.csv", "w");
    if (!out) {
        perror("erro ao criar resultados_monotona.csv");
        return 1;
    }
    fprintf(out, "heuristica,n,m,replicacao,tempo,iteracoes,valor,parametro\n");

    for (int im = 0; im < numM; im++) {
        int m = maquinas[im];
        for (int ir = 0; ir < numR; ir++) {
            double r = rs[ir];
            int n = (int) llround(pow((double) m, r));

            /* gera a instancia (tempos das tarefas) uma vez por (m, r) */
            Solucao *base = criar_solucao(n, m);
            gerar_tarefas(base);

            printf("Instancia: m=%d, r=%.1f -> n=%d\n", m, r, n);

            for (int ia = 0; ia < numAlpha; ia++) {
                double alpha = alphas[ia];

                for (int rep = 1; rep <= NUM_REPLICACOES; rep++) {
                    /* copia a instancia (mesmos p_i) e gera solucao inicial */
                    Solucao *s = criar_solucao(n, m);
                    memcpy(s->p, base->p, sizeof(int) * n);
                    solucao_inicial_gulosa(s);

                    clock_t inicio = clock();
                    ResultadoBusca res = busca_monotona_randomizada(s, alpha);
                    clock_t fim = clock();
                    double tempo = (double) (fim - inicio) / CLOCKS_PER_SEC;

                    fprintf(out, "monotona,%d,%d,%d,%.4f,%ld,%ld,%.1f\n",
                            n, m, rep, tempo, res.iteracoes, res.valor, alpha);

                    liberar_solucao(s);
                }
            }
            liberar_solucao(base);
        }
    }

    fclose(out);
    printf("Concluido. Resultados gravados em resultados_monotona.csv\n");
    return 0;
}
