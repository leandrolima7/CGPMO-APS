#include <math.h>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <float.h>

///Propagation delay(ns)
#define tdAND 1.7
#define tdNAND 1.8
#define tdOR 1.7
#define tdNOR 1.7
#define tdXOR 1.9
#define tdXNOR 72.5 ///mean value

///Probability that the output of a gate is in high level
#define P1AND 0.25
#define P1NAND 0.75
#define P1OR 0.75
#define P1NOR 0.25
#define P1XOR 0.5
#define P1XNOR 0.5

///Operational parameters
#define freq 50000000
#define Capa_Load 0.00000000005
#define Vcc 5

///CGP parameters
#define Max_Execucoes_Inde 1
#define NCOL 160  ///Number of columns
#define NLIN 1 ///Number of lines
#define LB NCOL /// LEVEL-BACK
#define PORTAS 6 ///Number of gates in Gamma set

#define Entradas 16 ///Number of Inputs
#define SAIDAS 9 ///Number of Outputs
#define NMUTACOES 0.05*NCOL ///Point Mutation
#define PArent_POPa 50 ///Number of parents
#define Limit 1000 ///Maximum in the front

///Seed
#define seed_do_circuito 31///a cada execução mudar isso


typedef struct Gene
{
    int lig1;
    int lig2;
    int porta; ///Gate
    int checked;
    int avaliar;
    int logical_out;
    int output[SAIDAS];
    double distancia;
    double atraso;
    double error;
    double power;
    double delay;
    struct Gene *prox1, *prox2;
}Gene;

typedef struct lista
{
    int val;
    struct lista * prox;
}lista , *vetor;

typedef struct lista_Euc_nearest
{
    int nearest;
    double val;
    struct lista_Euc_nearest * prox;
}lista_Euc_nearest;

void remover_R(lista ** lst);
void remover_R_Nearest(lista_Euc_nearest ** lst);
void inserir(lista** lst, int val);
void remover(lista ** lst, int val);
void insere_ord(lista** lst, int val);
void print_results(Gene **Rt, int linhaTabela, int colunaTabela,int number_solutions);
void insere_decrescebte_nearest(lista_Euc_nearest** lst, double error, int ind);
int busca(lista** lst, int val);

void FNDS(Gene **pop, lista** Front,int size_pop); ///Fast non dominated sort
void CrowDis(Gene **ind,lista **Front, double *Vet_CRD); ///Crowd Distance

int saida_porta(int r1, int r2, int porta);
int retornaValor(Gene Rt[NCOL], int linha[], int saida);
double dist_hamming(Gene ind[NCOL], int **tabela,int linhaTabela, int colunaTabela);
double delay(Gene ind[NCOL]);
double power(Gene ind[NCOL]);
void fitness(Gene **ind,int **tabela,int linhaTabela, int colunaTabela,int size_pop);
void Nu_Gates(Gene ind[NCOL]);

void criaInd(Gene ind[NCOL],int portas[], int saidas[]); ///Generates the first individual
void mutacao(Gene son[NCOL],Gene father[NCOL],int portas[PORTAS]); ///Mutation

int numAleatorio(int a, int b); ///Random number
unsigned long long int bin_to_dec(unsigned long long int bin); ///Binary to decimal
void sorteiaSaidas(int saidas[]); ///Defines output
void leTxt(int **tabela,int nL,int nC); ///Reads truth table
void imprimeTabela(int **vet, int nL,int nC); ///Print truth table
void selecionaPortas(int vet[PORTAS]);///Select gates

double cont_evaluation, percent, n_gates;
double Max_EVALUATION = 1000; /// Maximum circuit evaluation

int main()
{
    percent = pow(2,Entradas)*SAIDAS;

    int saidas[SAIDAS]; ///Output array of each individual
    int portas[PORTAS]; ///Available gates
    int linhaTabela,colunaTabela,number_solutions,actual_parent,size_pop;
    int **tabela = NULL; ///Truth Table

    int i1 = 0;
    int j1= 0;
    int trash;

    while(i1 < 1)
    {
        j1 = 0;
        while(j1 < 3)
        {
            if(j1 == 0)
                scanf("%d ",&linhaTabela);
            else if(j1 == 1)
            {
                scanf("%d",&colunaTabela);
                colunaTabela = colunaTabela + 1;
            }
            else
                scanf("%d",&trash);
            j1++;
        }
        i1++;
    }

    tabela = (int**)calloc(linhaTabela,sizeof(int*));
    for(int j = 0; j<linhaTabela;j++)
        tabela[j] = (int*)calloc(colunaTabela,sizeof(int));

    leTxt(tabela,linhaTabela,colunaTabela);
    selecionaPortas(portas); /// Seleciona as portas que farão parte do conjunto Gamma

    srand((unsigned)seed_do_circuito);

    printf("seed: %d Max_EVALUATION: %f\n",seed_do_circuito,Max_EVALUATION);

    size_pop = 1;
    cont_evaluation = 0;

    Gene **Rt = NULL;
    Gene **Pt = NULL;
    lista *tmp = NULL;

    vetor Front[1];
    Front[0] = NULL;

    Rt = (Gene**)calloc(size_pop,sizeof(Gene*)); ///Parent

    for(int z=0; z<size_pop;z++)
        Rt[z] = (Gene*)calloc(NCOL,sizeof(Gene));

    for(int j=0; j< NCOL; j++)
    {
        Rt[0][j].lig1 = -1;
        Rt[0][j].lig2 = -1;
        Rt[0][j].logical_out = -1;
        Rt[0][j].checked = -1;
        Rt[0][j].avaliar = 1;
        Rt[0][j].porta = -1;
        Rt[0][j].distancia = -1;
        Rt[0][j].error = -1;
        Rt[0][j].power = -1;
        Rt[0][j].delay = -1;
        Rt[0][j].atraso = -1;
        Rt[0][j].prox1 = NULL;
        Rt[0][j].prox2 = NULL;
    }

    for(int j=0; j< SAIDAS; j++)
        Rt[0][0].output[j] = 0;

    ///Creates the First individual **************************************************************************************************

    for(int j=0; j< size_pop; j++)
        criaInd(Rt[j],portas,saidas);

    ///**************************************************************************************************************************


    while(cont_evaluation < Max_EVALUATION)
    {
        fitness(Rt,tabela,linhaTabela,colunaTabela,size_pop);
        FNDS(Rt,Front,size_pop);

        tmp = Front[0];
        number_solutions = 0;

        while(tmp != NULL)
        {
            number_solutions++;
            tmp=tmp->prox;
        }

        cont_evaluation = cont_evaluation + size_pop;
        tmp = Front[0];

        Pt = (Gene**)calloc(number_solutions,sizeof(Gene*));///Offspring

        for(int z=0; z<number_solutions;z++)
            Pt[z] = (Gene*)calloc(NCOL,sizeof(Gene));
        for(int z = 0; z<number_solutions; z++)
        {
            for(int j=0; j< NCOL; j++)
            {
                Pt[z][j].lig1 = -1;
                Pt[z][j].lig2 = -1;
                Pt[z][j].logical_out = -1;
                Pt[z][j].checked = -1;
                Pt[z][j].avaliar = 0;
                Pt[z][j].porta = -1;
                Pt[z][j].atraso = -1;
                Pt[z][j].distancia = -1;
                Pt[z][j].power = -1;
                Pt[z][j].error = -1;
                Pt[z][j].delay = -1;
                Pt[z][j].prox1 = NULL;
                Pt[z][j].prox2 = NULL;
            }
            for(int j=0; j< SAIDAS; j++)
                Pt[z][0].output[j] = 0;
        }
        actual_parent = 0;

        for(int z = 0; z < number_solutions;z++)
        {
            for(int j=0; j< NCOL; j++)
            {
                Pt[actual_parent][j].lig1 = Rt[tmp->val][j].lig1;
                Pt[actual_parent][j].lig2 = Rt[tmp->val][j].lig2;
                Pt[actual_parent][j].checked = Rt[tmp->val][j].checked;
                Pt[actual_parent][j].logical_out = Rt[tmp->val][j].logical_out;
                Pt[actual_parent][j].porta = Rt[tmp->val][j].porta;
                Pt[actual_parent][j].atraso = Rt[tmp->val][j].atraso;
                Pt[actual_parent][j].error = Rt[tmp->val][j].error;
                Pt[actual_parent][j].distancia = Rt[tmp->val][j].distancia;
                Pt[actual_parent][j].delay = Rt[tmp->val][j].delay;
                Pt[actual_parent][j].power = Rt[tmp->val][j].power;
                if(Pt[actual_parent][j].lig1 < Entradas)
                    Pt[actual_parent][j].prox1 = NULL;
                else
                    Pt[actual_parent][j].prox1 = &Pt[actual_parent][Pt[actual_parent][j].lig1 - Entradas];

                if(Pt[actual_parent][j].lig2 < Entradas)
                    Pt[actual_parent][j].prox2 = NULL;
                else
                    Pt[actual_parent][j].prox2 = &Pt[actual_parent][Pt[actual_parent][j].lig2 - Entradas];
            }

            for(int u=0; u< SAIDAS; u++)
                Pt[actual_parent][0].output[u] = Rt[tmp->val][0].output[u];

            tmp = tmp->prox;
            actual_parent++;
        }

        if(cont_evaluation >= Max_EVALUATION)
        {
            print_results(Pt,linhaTabela,colunaTabela,number_solutions);
            break;
        }

        while(Front[0] != NULL)
            remover_R(&Front[0]);

        Front[0] = NULL;
        tmp = NULL;

        for(int z = 0;z< size_pop; z++)
            free(Rt[z]);

        free(Rt);
        Rt = NULL;

        Gene **Qt = (Gene**)calloc(number_solutions,sizeof(Gene*));
        for(int z=0; z<number_solutions;z++)
            Qt[z] = (Gene*)calloc(NCOL,sizeof(Gene));

        for(int z = 0; z<number_solutions; z++)
        {
            for(int j=0; j< NCOL; j++)
            {
                Qt[z][j].lig1 = -1;
                Qt[z][j].lig2 = -1;
                Qt[z][j].checked = -1;
                Qt[z][j].avaliar = 1;
                Qt[z][j].logical_out = -1;
                Qt[z][j].porta = -1;
                Qt[z][j].atraso = -1;
                Qt[z][j].distancia = -1;
                Qt[z][j].power = -1;
                Qt[z][j].error = -1;
                Qt[z][j].delay = -1;
                Qt[z][j].prox1 = NULL;
                Qt[z][j].prox2 = NULL;
            }
            for(int j=0; j< SAIDAS; j++)
                Qt[z][0].output[j] = 0;
        }

        size_pop = 2*number_solutions;

        for(int z = 0; z< number_solutions; z++)
            mutacao(Qt[z],Pt[z],portas);

        Rt = (Gene**)calloc(size_pop,sizeof(Gene*));

        for(int z=0; z<size_pop;z++)
            Rt[z] = (Gene*)calloc(NCOL,sizeof(Gene));
        for(int i = 0; i<number_solutions;i++)
        {
            for(int z = 0; z<NCOL;z++)
            {
                Rt[i][z].lig1 = Pt[i][z].lig1;
                Rt[i][z].lig2 = Pt[i][z].lig2;
                Rt[i][z].checked = Pt[i][z].checked;
                Rt[i][z].avaliar = Pt[i][z].avaliar;
                Rt[i][z].logical_out = Pt[i][z].logical_out;
                Rt[i][z].porta = Pt[i][z].porta;
                Rt[i][z].atraso = Pt[i][z].atraso;
                Rt[i][z].delay = Pt[i][z].delay;
                Rt[i][z].power = Pt[i][z].power;
                Rt[i][z].distancia = Pt[i][z].distancia;
                Rt[i][z].error = Pt[i][z].error;
                if(Rt[i][z].lig1 < Entradas)
                    Rt[i][z].prox1 = NULL;
                else
                    Rt[i][z].prox1 = &Rt[i][Rt[i][z].lig1 - Entradas];

                if(Rt[i][z].lig2 < Entradas)
                    Rt[i][z].prox2 = NULL;
                else
                    Rt[i][z].prox2 = &Rt[i][Rt[i][z].lig2 - Entradas];
            }

            for(int j=0; j< SAIDAS; j++)
                Rt[i][0].output[j] = Pt[i][0].output[j];
        }

        for(int i = number_solutions; i<size_pop;i++)
        {
            for(int z = 0; z<NCOL;z++)
            {
                Rt[i][z].lig1 = Qt[i-number_solutions][z].lig1;
                Rt[i][z].lig2 = Qt[i-number_solutions][z].lig2;
                Rt[i][z].checked = Qt[i-number_solutions][z].checked;
                Rt[i][z].avaliar = Qt[i-number_solutions][z].avaliar;
                Rt[i][z].logical_out = Qt[i-number_solutions][z].logical_out;
                Rt[i][z].porta = Qt[i-number_solutions][z].porta;
                Rt[i][z].atraso = Qt[i-number_solutions][z].atraso;
                Rt[i][z].delay = Qt[i-number_solutions][z].delay;
                Rt[i][z].power = Qt[i-number_solutions][z].power;
                Rt[i][z].distancia = Qt[i-number_solutions][z].distancia;
                Rt[i][z].error = Qt[i-number_solutions][z].error;
                if(Rt[i][z].lig1 < Entradas)
                    Rt[i][z].prox1 = NULL;
                else
                    Rt[i][z].prox1 = &Rt[i][Rt[i][z].lig1 - Entradas];

                if(Rt[i][z].lig2 < Entradas)
                    Rt[i][z].prox2 = NULL;
                else
                    Rt[i][z].prox2 = &Rt[i][Rt[i][z].lig2 - Entradas];
            }

            for(int j=0; j< SAIDAS; j++)
                Rt[i][0].output[j] = Qt[i-number_solutions][0].output[j];
        }

        for(int z = 0;z< number_solutions; z++)
            free(Qt[z]);
        free(Qt);
        Qt = NULL;

        for(int z = 0;z< number_solutions; z++)
            free(Pt[z]);
        free(Pt);
        Pt = NULL;
    }

    for(int z = 0;z< number_solutions; z++)
        free(Pt[z]);
    free(Pt);
    Pt = NULL;

    for(int z = 0;z< size_pop; z++)
        free(Rt[z]);
    free(Rt);
    Rt = NULL;

    for(int z = 0;z<linhaTabela; z++)
        free(tabela[z]);
    free(tabela);
    tabela = NULL;

    return 0;
}

void print_results(Gene **Rt, int linhaTabela, int colunaTabela,int number_solutions)
{
    printf("\n");
    printf("Numero avaliacoes: %f\n",cont_evaluation);
    printf("\n");


    printf("Num Ind: %d\n",number_solutions);
    printf("\nD_H  | Delay | Potencia\n\n");

    for(int i = 0; i < number_solutions; i++)
        printf("%f %f %f\n",Rt[i][0].error/percent,Rt[i][0].delay,Rt[i][0].power);

    printf("\n");

    printf("-------------------------------------------------------------------------------------------------------------------------\n");
}

void FNDS(Gene **Rt, lista** Front, int size_pop)
{
    int front_size = 0;
    int np[size_pop];

    for(int p = 0;p<size_pop;p++)
        np[p] = 0;

    for(int p=0;p<size_pop;p++)
    {
        for(int q=p+1;q<size_pop;q++)
        {
            if(((Rt[p][0].power <= Rt[q][0].power) && (Rt[p][0].error <= Rt[q][0].error) && (Rt[p][0].delay <= Rt[q][0].delay)) && ((Rt[p][0].power < Rt[q][0].power) || (Rt[p][0].error < Rt[q][0].error) || (Rt[p][0].delay < Rt[q][0].delay)))
                np[q] = np[q] + 1;
            else if (((Rt[q][0].power <= Rt[p][0].power) && (Rt[q][0].error <= Rt[p][0].error) && (Rt[q][0].delay <= Rt[p][0].delay)) && ((Rt[q][0].power < Rt[p][0].power) || (Rt[q][0].error < Rt[p][0].error) || (Rt[q][0].delay < Rt[p][0].delay)))
                np[p] = np[p] + 1;
        }
        if(np[p] == 0)
        {
            front_size++;
            inserir(&Front[0],p);
        }
    }

    if(size_pop > Limit)
    {
        lista *tmp = Front[0];

        double Vet_CRD[size_pop+1];

        for(int j = 0;j<size_pop+1;j++)
            Vet_CRD[j] = 0;

        CrowDis(Rt,&Front[0],Vet_CRD);
        int contador = 0;
        lista *avaliar = NULL;

        lista_Euc_nearest *Maiores_distancias = NULL;;
        lista_Euc_nearest *aux_Maiores_distancias = NULL;

        for(int z = 0; z < (int)Vet_CRD[0];z++)
        {
            Rt[tmp->val][0].distancia = Vet_CRD[tmp->val+1];
            insere_decrescebte_nearest(&Maiores_distancias,Rt[tmp->val][0].distancia,tmp->val);
            tmp = tmp->prox;
        }

        aux_Maiores_distancias = Maiores_distancias;

        while((contador < Limit) && (aux_Maiores_distancias != NULL))
        {
            inserir(&avaliar,aux_Maiores_distancias->nearest);
            contador++;
            aux_Maiores_distancias = aux_Maiores_distancias-> prox;
        }

        tmp = avaliar;

        while(Front[0] != NULL)
            remover_R(&Front[0]);
        Front[0] = NULL;

        front_size = 0;

        while(tmp != NULL)
        {
            inserir(&Front[0],tmp->val);
            front_size++;
            tmp = tmp->prox;
        }

        while(avaliar!= NULL)
            remover_R(&avaliar);
        free(avaliar);


        while(Maiores_distancias!= NULL)
            remover_R_Nearest(&Maiores_distancias);
        free(Maiores_distancias);
    }
}

void CrowDis(Gene **ind,lista **Front, double *Vet_CRD)
{
    int objectives = 3;
    lista *tmp_number_solutions = *Front;
    int number_solutions = 0;

    while(tmp_number_solutions != NULL)
    {
        number_solutions++;
        tmp_number_solutions=tmp_number_solutions->prox;
    }
    tmp_number_solutions = NULL;

    Vet_CRD[0] = (double)number_solutions;

    for(int solution  = 0; solution < objectives; solution++)
    {
        if(solution == 0)
        {
            int i,j,aux1;
            lista *tmp = *Front;
            lista *tmp2 = tmp->prox;
            lista *tmp_min = NULL;

            for (i = 0; i < (number_solutions-1); i++)
            {
                tmp_min = tmp;
                for (j = (i+1); j < number_solutions; j++)
                {
                    if(ind[tmp2->val][0].error < ind[tmp_min->val][0].error)
                       tmp_min = tmp2;
                    tmp2=tmp2->prox;
                }
                if(ind[tmp->val][0].error != ind[tmp_min->val][0].error)
                {
                    aux1 = tmp->val;
                    tmp->val = tmp_min->val;
                    tmp_min->val = aux1;
                }
                tmp=tmp->prox;
                tmp2 = tmp->prox;
            }

            tmp2 = *Front;

            if(number_solutions == 1)
            {
                Vet_CRD[tmp2->val+1] = DBL_MAX;

                continue;
            }
            else if(number_solutions == 2)
            {
                Vet_CRD[tmp2->val+1] = DBL_MAX;
                Vet_CRD[tmp2->prox->val+1] = DBL_MAX;
                continue;
            }
            else
            {
                tmp = *Front;
                Vet_CRD[tmp->val+1] = DBL_MAX;
                while(tmp->prox != NULL)
                {
                    tmp = tmp -> prox;
                }
                Vet_CRD[tmp->val+1] = DBL_MAX;
                tmp = *Front;
                lista *tmp_before = tmp;
                tmp = tmp -> prox;
                lista *tmp_after = tmp->prox;
                int r;
                for(r = 2;r<number_solutions;r++)
                {
                    Vet_CRD[tmp->val+1] = Vet_CRD[tmp->val+1] + (ind[tmp_after->val][0].error - ind[tmp_before->val][0].error);
                    tmp_before = tmp;
                    tmp = tmp->prox;
                    tmp_after = tmp->prox;
                }
                tmp_before = NULL;
                tmp_after = NULL;
            }
            tmp = NULL;
            tmp2 = NULL;
            tmp_min = NULL;
        }

        if(solution == 1)
        {
            int i,j,aux2;
            lista *tmp = *Front;
            lista *tmp2 = tmp->prox;
            lista *tmp_min = NULL;;
            for (i = 0; i < (number_solutions-1); i++)
            {
                tmp_min = tmp;
                for (j = (i+1); j < number_solutions; j++)
                {
                    if(ind[tmp2->val][0].delay < ind[tmp_min->val][0].delay)
                       tmp_min = tmp2;
                    tmp2=tmp2->prox;
                }
                if (ind[tmp->val][0].delay != ind[tmp_min->val][0].delay)
                {
                    aux2 = tmp->val;
                    tmp->val = tmp_min->val;
                    tmp_min->val = aux2;
                }
                tmp=tmp->prox;
                tmp2 = tmp->prox;
            }

            tmp2 = *Front;

            if(number_solutions == 1)
            {
                Vet_CRD[tmp2->val+1] = DBL_MAX;
                continue;
            }
            else if(number_solutions == 2)
            {
                Vet_CRD[tmp2->val+1] = DBL_MAX;
                Vet_CRD[tmp2->prox->val+1] = DBL_MAX;
                continue;
            }
            else
            {
                tmp = *Front;
                Vet_CRD[tmp->val+1] = DBL_MAX;
                while(tmp->prox != NULL)
                {
                    tmp = tmp -> prox;
                }
                Vet_CRD[tmp->val+1] = DBL_MAX;
                tmp = *Front;
                tmp = tmp -> prox;
                lista *tmp_after = tmp->prox;
                lista *tmp_before = *Front;
                int r;
                for(r = 2;r<number_solutions;r++)
                {
                    Vet_CRD[tmp->val+1] = Vet_CRD[tmp->val+1] + (ind[tmp_after->val][0].delay - ind[tmp_before->val][0].delay);
                    tmp_before = tmp;
                    tmp = tmp->prox;
                    tmp_after = tmp->prox;
                }
                tmp_before = NULL;
                tmp_after = NULL;
            }
            tmp = NULL;
            tmp2 = NULL;
            tmp_min = NULL;
        }
        if(solution == 2)
        {
            int i,j,aux2;
            lista *tmp = *Front;
            lista *tmp2 = tmp->prox;
            lista *tmp_min = NULL;;
            for (i = 0; i < (number_solutions-1); i++)
            {
                tmp_min = tmp;
                for (j = (i+1); j < number_solutions; j++)
                {
                    if(ind[tmp2->val][0].power < ind[tmp_min->val][0].power)
                       tmp_min = tmp2;
                    tmp2=tmp2->prox;
                }
                if (ind[tmp->val][0].power != ind[tmp_min->val][0].power)
                {
                    aux2 = tmp->val;
                    tmp->val = tmp_min->val;
                    tmp_min->val = aux2;
                }
                tmp=tmp->prox;
                tmp2 = tmp->prox;
            }

            tmp2 = *Front;

            if(number_solutions == 1)
            {
                Vet_CRD[tmp2->val+1] = DBL_MAX;
                continue;
            }
            else if(number_solutions == 2)
            {
                Vet_CRD[tmp2->val+1] = DBL_MAX;
                Vet_CRD[tmp2->prox->val+1] = DBL_MAX;
                continue;
            }
            else
            {
                tmp = *Front;
                Vet_CRD[tmp->val+1] = DBL_MAX;
                while(tmp->prox != NULL)
                {
                    tmp = tmp -> prox;
                }
                Vet_CRD[tmp->val+1] = DBL_MAX;
                tmp = *Front;
                tmp = tmp -> prox;
                lista *tmp_after = tmp->prox;
                lista *tmp_before = *Front;
                int r;
                for(r = 2;r<number_solutions;r++)
                {
                    Vet_CRD[tmp->val+1] = Vet_CRD[tmp->val+1] + (ind[tmp_after->val][0].power - ind[tmp_before->val][0].power);
                    tmp_before = tmp;
                    tmp = tmp->prox;
                    tmp_after = tmp->prox;
                }
                tmp_before = NULL;
                tmp_after = NULL;
            }
            tmp = NULL;
            tmp2 = NULL;
            tmp_min = NULL;
        }
    }
}

int numAleatorio(int a, int b)
{
    return (a + rand()%(b - a + 1));
}

void leTxt(int **tabela,int nL,int nC)
{
   for(int i = 0 ; i < nL; i++)
    {
        for(int j = 0 ; j < nC; j++)
            scanf("%d ",&tabela[i][j]);
    }
}

void imprimeTabela(int **vet, int nL,int nC)
{
    int i = 0;
    int j = 0;

    for(i = 0; i < nL; i++)
    {
        for(j = 0; j < nC; j++)
        {
            printf("%d  ",vet[i][j]);
        }
        printf("\n");
    }
}

void selecionaPortas(int vet[])
{
    int contador = 0;
    int porta;

    while(contador < PORTAS)
    {
        scanf("%d",&porta);
        vet[contador] = porta;
        contador++;
    }

    contador = 0;
    printf("Portas escolhidas: ");
    while(contador < PORTAS)
    {
        printf("%d ",vet[contador]);
        contador++;
    }
    printf("\n");

}

void sorteiaSaidas(int saidas[])
{
    int i;
    for(i = 0; i < SAIDAS; i++)
        saidas[i] = numAleatorio(0, (NLIN*NCOL)+(Entradas) - 1);
}

unsigned long long int bin_to_dec(unsigned long long int bin)
{
    unsigned long long int total  = 0;
    unsigned long long int potenc = 1;

    while(bin > 0)
    {
        total += bin % 10 * potenc;
        bin = bin / 10;
        potenc = potenc * 2;
    }
    return total;
}

void criaInd(Gene ind[NCOL],int portas[], int saidas[])
{
    int auxiliar, aleatorio;

    sorteiaSaidas(saidas);

    for(int o = 0; o < SAIDAS; o++)
        ind[0].output[o] = saidas[o];

    for(int j = 0; j < NCOL; j++)
    {
        if(j == 0 || LB == 0)
        {
            ind[j].lig1 = numAleatorio(0, Entradas - 1);
            ind[j].lig2 = numAleatorio(0, Entradas - 1);
            ind[j].prox1 = NULL;
            ind[j].prox2 = NULL;
        }
        else
        {
            if(LB >= j || LB == -1)
            {
                ind[j].lig1 = numAleatorio(0, (Entradas) + (NLIN*j) - 1);
                ind[j].lig2 = numAleatorio(0, (Entradas) + (NLIN*j) - 1);
                if(ind[j].lig1 < Entradas)
                {
                    ind[j].prox1 = NULL;
                }
                else
                {
                    ind[j].prox1 = &ind[ind[j].lig1 - Entradas];
                }
                if(ind[j].lig2 < Entradas)
                {
                    ind[j].prox2 = NULL;
                }
                else
                {
                    ind[j].prox2 = &ind[ind[j].lig2 - Entradas];
                }
            }
            else
            {
                if(LB < j)
                {
                    aleatorio = numAleatorio(0,(Entradas)+ LB*NLIN - 1);
                    if(aleatorio >= 0  && aleatorio < (Entradas-1))
                    {
                        ind[j].lig1 = aleatorio;
                        if(ind[j].lig1 < Entradas)
                        {
                            ind[j].prox1 = NULL;
                        }
                        else
                        {
                            ind[j].prox1 = &ind[ind[j].lig1 - Entradas];
                        }
                    }
                    else
                    {
                        ind[j].lig1 = numAleatorio((Entradas) + (j - LB)*NLIN, (Entradas) + (NLIN*j) - 1);
                        if(ind[j].lig1 < Entradas)
                        {
                            ind[j].prox1 = NULL;
                        }
                        else
                        {
                            ind[j].prox1 = &ind[ind[j].lig1 - Entradas];
                        }
                    }

                    aleatorio = numAleatorio(0,(Entradas)+ LB*NLIN - 1);

                    if(aleatorio >= 0  && aleatorio < (Entradas-1))
                    {
                        ind[j].lig2 = aleatorio;
                        if(ind[j].lig2 < Entradas)
                        {
                            ind[j].prox2 = NULL;
                        }
                        else
                        {
                            ind[j].prox2 = &ind[ind[j].lig2 - Entradas];
                        }
                    }
                    else
                    {
                        ind[j].lig2 = numAleatorio((Entradas) + (j - LB)*NLIN, (Entradas) + (NLIN*j) - 1);
                        if(ind[j].lig2 < Entradas)
                        {
                            ind[j].prox2 = NULL;
                        }
                        else
                        {
                            ind[j].prox2 = &ind[ind[j].lig2 - Entradas];
                        }
                    }
                }
            }
        }
        auxiliar = numAleatorio(0, (PORTAS-1));
        ind[j].porta = portas[auxiliar];

        if(ind[j].porta == 0)
            ind[j].atraso = tdAND;
        else if(ind[j].porta == 1)
            ind[j].atraso = tdOR;
        else if(ind[j].porta == 2)
            ind[j].atraso = tdNOR;
        else if(ind[j].porta == 3)
            ind[j].atraso = tdNAND;
        else if(ind[j].porta == 4)
            ind[j].atraso = tdXOR;
        else if(ind[j].porta == 5)
            ind[j].atraso = tdXNOR;
    }
}

void mutacao(Gene son[NCOL],Gene father[NCOL],int portas[PORTAS])
{
    int teste = (int)NMUTACOES;
    int auxiliar;
    int gen, coluna, aleat;

    if(teste < 1)
        teste = 1;

    for(int z = 0; z < NCOL; z++)
    {
        son[z].lig1 = father[z].lig1;
        son[z].lig2 = father[z].lig2;
        son[z].porta = father[z].porta;
        son[z].atraso = father[z].atraso;
        son[z].distancia = -1;
        son[z].power = father[z].power;
        son[z].error = father[z].error;
        son[z].delay = father[z].delay;

        if(son[z].lig1 < Entradas)
            son[z].prox1 = NULL;
        else
            son[z].prox1 = &son[son[z].lig1 - Entradas];

        if(son[z].lig2 < Entradas)
            son[z].prox2 = NULL;
        else
            son[z].prox2 = &son[son[z].lig2 - Entradas];
    }

    for(int j=0; j< SAIDAS; j++)
        son[0].output[j] = father[0].output[j];

    int u = 0;

    lista *lst_nos = NULL;
    lista *tmp_nos = NULL;

    for(int i=0;i<SAIDAS;i++)
    {
        if(son[0].output[i] < Entradas)
            continue;

        u = son[0].output[i] - Entradas;

        inserir(&lst_nos,u);
        tmp_nos=lst_nos;

        do
        {
            if(son[(int)tmp_nos->val].prox1 != NULL)
                inserir(&lst_nos,son[(int)tmp_nos->val].lig1 - Entradas);

            if(son[(int)tmp_nos->val].prox2 != NULL)
                inserir(&lst_nos,son[(int)tmp_nos->val].lig2 - Entradas);

            tmp_nos = tmp_nos -> prox;

        }while(tmp_nos!=NULL);
    }

    tmp_nos = lst_nos;

    int checa = 0;

    for(int z = 0; z < teste; z++)
    {

        gen = numAleatorio((Entradas), (Entradas) + NCOL*NLIN - 1);
        coluna  = (gen - (Entradas))/(NLIN);

        if(checa == 0)
        {
            checa = busca(&lst_nos,coluna);
        }

        aleat = numAleatorio(0, 3);

        if(aleat == 0)
        {
            if(coluna == 0 || LB == 0)
            {
                son[coluna].lig1 = numAleatorio(0, Entradas - 1) ;
                if(son[coluna].lig1 < Entradas)
                {
                    son[coluna].prox1 = NULL;
                }
                else
                {
                    son[coluna].prox1 = &son[son[coluna].lig1 - Entradas];
                }
            }
            else
            {
                if(LB >= coluna || LB == -1)
                {
                    son[coluna].lig1 = numAleatorio(0, (Entradas) + (NLIN*coluna) - 1);
                    if(son[coluna].lig1 < Entradas)
                    {
                        son[coluna].prox1 = NULL;
                    }
                    else
                    {
                        son[coluna].prox1 = &son[son[coluna].lig1 - Entradas];
                    }
                }
                else
                {
                    if(LB < coluna && LB > 0)
                    {
                        int aleatorio = numAleatorio(0,(Entradas)+ LB*NLIN - 1);
                        if(aleatorio >= 0  && aleatorio < (Entradas-1))
                        {
                            son[coluna].lig1 = aleatorio;
                            if(son[coluna].lig1 < Entradas)
                            {
                                son[coluna].prox1 = NULL;
                            }
                            else
                            {
                                son[coluna].prox1 = &son[son[coluna].lig1 - Entradas];
                            }
                        }
                        else
                        {
                            son[coluna].lig1 = numAleatorio((Entradas) + (coluna - LB)*NLIN, (Entradas) + (NLIN*coluna) - 1);
                            if(son[coluna].lig1 < Entradas)
                            {
                                son[coluna].prox1 = NULL;
                            }
                            else
                            {
                                son[coluna].prox1 = &son[son[coluna].lig1 - Entradas];
                            }
                        }
                    }
                }
            }
        }
        else
        {
            if(aleat == 1)
            {
                if(coluna == 0 || LB == 0)
                {
                    son[coluna].lig2 = numAleatorio(0, Entradas - 1) ;
                    if(son[coluna].lig2 < Entradas)
                    {
                        son[coluna].prox2 = NULL;
                    }
                    else
                    {
                        son[coluna].prox2 = &son[son[coluna].lig2 - Entradas];
                    }
                }
                else
                {
                    if(LB >= coluna || LB == -1)
                    {
                        son[coluna].lig2 = numAleatorio(0, (Entradas) + (NLIN*coluna) - 1);
                        if(son[coluna].lig2 < Entradas)
                        {
                            son[coluna].prox2 = NULL;
                        }
                        else
                        {
                            son[coluna].prox2 = &son[son[coluna].lig2 - Entradas];
                        }
                    }
                    else
                    {
                        if(LB < coluna && LB > 0)
                        {
                            int aleatorio = numAleatorio(0,(Entradas)+ LB*NLIN - 1);
                            if(aleatorio >= 0  && aleatorio < (Entradas-1))
                            {
                                son[coluna].lig2 = aleatorio;
                                if(son[coluna].lig2 < Entradas)
                                {
                                    son[coluna].prox2 = NULL;
                                }
                                else
                                {
                                    son[coluna].prox2 = &son[son[coluna].lig2 - Entradas];
                                }
                            }
                            else
                            {
                                son[coluna].lig2 = numAleatorio((Entradas) + (coluna - LB)*NLIN, (Entradas) + (NLIN*coluna) - 1);
                                if(son[coluna].lig2 < Entradas)
                                {
                                    son[coluna].prox2 = NULL;
                                }
                                else
                                {
                                    son[coluna].prox2 = &son[son[coluna].lig2 - Entradas];
                                }
                            }
                        }
                    }
                }
            }
            else if(aleat == 2)
            {
                auxiliar = numAleatorio(0, (PORTAS-1));
                son[coluna].porta = portas[auxiliar];

                if(son[coluna].porta == 0)
                    son[coluna].atraso = tdAND;
                else if(son[coluna].porta == 1)
                    son[coluna].atraso = tdOR;
                else if(son[coluna].porta == 2)
                    son[coluna].atraso = tdNOR;
                else if(son[coluna].porta == 3)
                    son[coluna].atraso = tdNAND;
                else if(son[coluna].porta == 4)
                    son[coluna].atraso = tdXOR;
                else if(son[coluna].porta == 5)
                    son[coluna].atraso = tdXNOR;
            }
            else if(aleat == 3)
            {
                checa = 1;
                auxiliar = numAleatorio(0, (SAIDAS-1));
                son[0].output[auxiliar] = numAleatorio(Entradas, (NLIN*NCOL)+(Entradas) - 1);;
            }
        }
    }

    son[0].avaliar = checa;

    tmp_nos = NULL;

    while(lst_nos != NULL)
        remover_R(&lst_nos);
    free(lst_nos);
    lst_nos = NULL;
}

int saida_porta(int r1, int r2, int porta)
{
    /// XOR variables
    int r1n;
    int r2n;
    int op1;
    int op2;
    ///
    switch(porta)
    {
        case 0: /// AND
            if(r1 == 1 && r2 == 1)
                return 1;
            else
                return 0;
            break;

        case 1: /// OR
            if(r1 == 1 || r2 == 1)
                return 1;
            else
                return 0;
            break;

        case 2: /// NOR
            if(r1 == 1 || r2 == 1)
                return 0;
            else
                return 1;
            break;

        case 3: /// NAND
            if(r1 == 1 && r2 == 1)
                return 0;
            else
                return 1;
            break;

        case 4: /// XOR
            if(r1 == 1)
                r1n = 0;
            else
                r1n = 1;

            if(r2 == 1)
                r2n = 0;
            else
                r2n = 1;

            if(r1n == 1 && r2 == 1)
                op1 = 1;
            else
                op1 = 0;

            if(r1 == 1 && r2n == 1)
                op2 = 1;
            else
                op2 = 0;

            if(op1 == 1 || op2 == 1)
                return 1;
            else
                return 0;
            break;

        case 5: /// XNOR
            if(r1 == r2)
                return 1;
            else
                return 0;
            break;

        case 6: /// A AND NOT B
            if(r1 == 1 && r2 == 0)
                return 1;
            else
                return 0;
            break;

        case 7: /// A
            return r1;
            break;

        case 8: /// B AND NOT A
            if(r1 == 0 && r2 == 1)
                return 1;
            else
                return 0;
            break;

        case 9: /// B
            return r2;
            break;

        case 10: /// NOT B
            if(r2 == 1)
                return 0;
            else
                return 1;
            break;

        case 11: /// A OR NOT B
            if(r1 == 0 && r2 == 1)
                return 0;
            else
                return 1;
            break;

        case 12: /// NOT A
            if(r1 == 1)
                return 0;
            else
                return 1;
            break;
        case 13: /// B OR NOT A
            if(r1 == 1 && r2 == 0)
                return 0;
            else
                return 1;
            break;
    }
    return 0;
}

int retornaValor(Gene ind[NCOL], int linha[], int saida)
{
    if(saida >= 0 && saida < (Entradas))
    {
        return linha[saida];
    }
    else
    {
        if(saida >= (Entradas))
        {
            int numeroColuna = (saida - (Entradas))/(NLIN);

            if(ind[numeroColuna].checked == -1)
            {
                int r1 = retornaValor(ind, linha, ind[numeroColuna].lig1);
                int r2 = retornaValor(ind, linha, ind[numeroColuna].lig2);

                ind[numeroColuna].logical_out = saida_porta(r1, r2, ind[numeroColuna].porta);

                ind[numeroColuna].checked = 1;

                return ind[numeroColuna].logical_out;
            }
            else if(ind[numeroColuna].checked == 1)
                return ind[numeroColuna].logical_out;
            else
            {
                cont_evaluation = DBL_MAX;
                printf("Erro de execucao! \n");
            }
        }
    }
    return 0;
}

double dist_hamming(Gene ind[NCOL], int **tabela,int linhaTabela, int colunaTabela)
{
    double Number_Incorrect_bits = 0;
    double resultados[SAIDAS];

    for(int z_a = 0; z_a < linhaTabela; z_a++)
    {
        for(int j_a = 0; j_a < NCOL; j_a++)
            ind[j_a].checked = -1;

        for(int j_a = 0; j_a < SAIDAS; j_a++)
            resultados[j_a] = retornaValor(ind, tabela[z_a], ind[0].output[j_a]);

        double cont = pow(10, SAIDAS-1);
        double pontuacao = 0;
        for(int j_b = 0; j_b < SAIDAS; j_b++)
        {
            pontuacao += cont*resultados[j_b];
            cont = cont/10;
        }

        pontuacao = (double)bin_to_dec((unsigned long long int)pontuacao);
        double dist = 0;
        int val = tabela[z_a][colunaTabela-1]^(int)pontuacao;
        while(val)
        {
            ++dist;
            val &= val - 1;
        }
        Number_Incorrect_bits = Number_Incorrect_bits + dist;
    }
    return Number_Incorrect_bits;
}

double delay(Gene ind[NCOL])
{
    int u, adj[2];
    double sum_d[SAIDAS];
    double numero_max = 0;

    u = 0;

    for(int j=0;j<2;j++)
       adj[j] = 0;

    for(int j=0;j<SAIDAS;j++)
       sum_d[j] = 0;

    for(int i=0;i<SAIDAS;i++)
    {
        if(ind[0].output[i] < Entradas)
            continue;

        int tam = (ind[0].output[i] - Entradas + 2);

        double d[tam];

        for(int j=0;j<tam;j++)
            d[j] = DBL_MAX;

        u = ind[0].output[i] - Entradas;
        d[u+1] = 0;

        lista *lst = NULL;
        lista *tmp = NULL;

        inserir(&lst,u);
        tmp = lst;

        do
        {
            if(ind[(int)tmp->val].prox1 != NULL)
                inserir(&lst,ind[(int)tmp->val].lig1 - Entradas);

            if(ind[(int)tmp->val].prox2 != NULL)
                inserir(&lst,ind[(int)tmp->val].lig2 - Entradas);

            tmp = tmp -> prox;
        }while(tmp!=NULL);

        int flag = 0;
        tmp = lst;

        while(lst != NULL)
        {
            u = lst->val;

            if(ind[u].lig1 < Entradas)
                adj[0] = 0;
            else
                adj[0] = ind[u].lig1-Entradas+1;

            if(ind[u].lig2 < Entradas)
                adj[1] = 0;
            else
                adj[1] = ind[u].lig2-Entradas+1;

            for(int j=0;j<2;j++)
            {
                if(flag == 0)
                {
                    if(d[adj[j]] > d[u+1] + ind[u].atraso)
                        d[adj[j]] = d[u+1] + ind[u].atraso;
                }
                else
                {
                    if(d[adj[j]] < d[u+1] + ind[u].atraso)
                        d[adj[j]] = d[u+1] + ind[u].atraso;
                }
            }
            flag = 1;
            if(flag == 1)
            {
                for(int j=0;j<tam;j++)
                {
                    if(d[j] == DBL_MAX)
                        d[j] = -1;
                }
            }

            remover(&lst,u);
        }

        numero_max = d[0];
        for (int w = 0;w < tam; w++)
            if (d[w] > numero_max)
                numero_max = d[w];
        sum_d[i] = numero_max;

        while(lst != NULL)
            remover_R(&lst);
        free(lst);

        lst = NULL;
        tmp = lst;
    }

    numero_max = sum_d[0];
    for (int i = 0;i < SAIDAS; ++i)
        if (sum_d[i] > numero_max)
            numero_max = sum_d[i];
    return numero_max;
}

double power(Gene ind[NCOL])
{
    int u = 0;
    double switching_component_of_power = 0;
    lista *g_lst = NULL;
    lista *g_tmp = NULL;

    for(int i=0;i<SAIDAS;i++)
    {

        if(ind[0].output[i] < Entradas)
            continue;

        u = ind[0].output[i] - Entradas;

        lista *lst = NULL;
        lista *tmp1 = NULL;


        inserir(&lst,u);
        tmp1 = lst;


        do
        {
            if(ind[(int)tmp1->val].prox1 != NULL)
                inserir(&lst,ind[(int)tmp1->val].lig1 - Entradas);

            if(ind[(int)tmp1->val].prox2 != NULL)
                inserir(&lst,ind[(int)tmp1->val].lig2 - Entradas);

            tmp1 = tmp1 -> prox;
        }while(tmp1!=NULL);


        g_tmp = lst;

        while(g_tmp != NULL)
        {
            insere_ord(&g_lst,g_tmp->val);
            g_tmp = g_tmp->prox;
        }

        tmp1 = NULL;

        while(lst != NULL)
            remover_R(&lst);
        free(lst);

        lst = NULL;
    }


    if(g_lst == NULL)
        return switching_component_of_power;
    else
    {
        g_tmp = g_lst;

        while(g_tmp->prox != NULL)
            g_tmp = g_tmp->prox;

        double *aux_vet = NULL;
        int size_aux_vet = g_tmp->val+1;

        aux_vet = (double*)calloc(size_aux_vet,sizeof(double));

        for(int i = 0 ;i< size_aux_vet; i++)
            aux_vet[i] = -1;

        lista *tmp = g_lst;

        do
        {

            if(ind[(int)tmp->val].prox1 == NULL && ind[(int)tmp->val].prox2 == NULL)
            {
                if(ind[(int)tmp->val].porta == 0)
                    aux_vet[(int)tmp->val] = P1AND;
                else if(ind[(int)tmp->val].porta == 1)
                    aux_vet[(int)tmp->val] = P1OR;
                else if(ind[(int)tmp->val].porta == 2)
                    aux_vet[(int)tmp->val] = P1NOR;
                else if(ind[(int)tmp->val].porta == 3)
                    aux_vet[(int)tmp->val] = P1NAND;
                else if(ind[(int)tmp->val].porta == 4)
                    aux_vet[(int)tmp->val] = P1XOR;
                else if(ind[(int)tmp->val].porta == 5)
                    aux_vet[(int)tmp->val] = P1XNOR;
            }
            else if(ind[(int)tmp->val].prox1 != NULL && ind[(int)tmp->val].prox2 != NULL)
            {
                double prob_lig1,prob_lig2;
                prob_lig1 = aux_vet[ind[(int)tmp->val].lig1 - Entradas];
                prob_lig2 = aux_vet[ind[(int)tmp->val].lig2 - Entradas];

                if(ind[(int)tmp->val].porta == 0) ///AND
                    aux_vet[(int)tmp->val] = prob_lig1*prob_lig2;
                else if(ind[(int)tmp->val].porta == 1) ///OR
                    aux_vet[(int)tmp->val] = (1-(1-prob_lig1)*(1-prob_lig2));
                else if(ind[(int)tmp->val].porta == 2) ///NOR
                    aux_vet[(int)tmp->val] = (1-prob_lig1)*(1-prob_lig2);
                else if(ind[(int)tmp->val].porta == 3) ///NAND
                    aux_vet[(int)tmp->val] = 1-(prob_lig1*prob_lig2);
                else if(ind[(int)tmp->val].porta == 4) ///XOR
                    aux_vet[(int)tmp->val] = 1 - ((1-prob_lig1)*(1-prob_lig2) + (prob_lig1*prob_lig2));
                else if(ind[(int)tmp->val].porta == 5) ///XNOR
                    aux_vet[(int)tmp->val] = (1-prob_lig1)*(1-prob_lig2) + (prob_lig1*prob_lig2);

            }
            else if(ind[(int)tmp->val].prox1 == NULL && ind[(int)tmp->val].prox2 != NULL)
            {
                double prob_lig1,prob_lig2;
                prob_lig1 = 0.5;
                prob_lig2 = aux_vet[ind[(int)tmp->val].lig2 - Entradas];

                if(ind[(int)tmp->val].porta == 0) ///AND
                    aux_vet[(int)tmp->val] = prob_lig1*prob_lig2;
                else if(ind[(int)tmp->val].porta == 1) ///OR
                    aux_vet[(int)tmp->val] = (1-(1-prob_lig1)*(1-prob_lig2));
                else if(ind[(int)tmp->val].porta == 2) ///NOR
                    aux_vet[(int)tmp->val] = (1-prob_lig1)*(1-prob_lig2);
                else if(ind[(int)tmp->val].porta == 3) ///NAND
                    aux_vet[(int)tmp->val] = 1-(prob_lig1*prob_lig2);
                else if(ind[(int)tmp->val].porta == 4) ///XOR
                    aux_vet[(int)tmp->val] = 1 - ((1-prob_lig1)*(1-prob_lig2) + (prob_lig1*prob_lig2));
                else if(ind[(int)tmp->val].porta == 5) ///XNOR
                    aux_vet[(int)tmp->val] = (1-prob_lig1)*(1-prob_lig2) + (prob_lig1*prob_lig2);

            }
            else if(ind[(int)tmp->val].prox1 != NULL && ind[(int)tmp->val].prox2 == NULL)
            {
                double prob_lig1,prob_lig2;
                prob_lig1 = aux_vet[ind[(int)tmp->val].lig1 - Entradas];
                prob_lig2 = 0.5;

                if(ind[(int)tmp->val].porta == 0) ///AND
                    aux_vet[(int)tmp->val] = prob_lig1*prob_lig2;
                else if(ind[(int)tmp->val].porta == 1) ///OR
                    aux_vet[(int)tmp->val] = (1-(1-prob_lig1)*(1-prob_lig2));
                else if(ind[(int)tmp->val].porta == 2) ///NOR
                    aux_vet[(int)tmp->val] = (1-prob_lig1)*(1-prob_lig2);
                else if(ind[(int)tmp->val].porta == 3) ///NAND
                    aux_vet[(int)tmp->val] = 1-(prob_lig1*prob_lig2);
                else if(ind[(int)tmp->val].porta == 4) ///XOR
                    aux_vet[(int)tmp->val] = 1 - ((1-prob_lig1)*(1-prob_lig2) + (prob_lig1*prob_lig2));
                else if(ind[(int)tmp->val].porta == 5) ///XNOR
                    aux_vet[(int)tmp->val] = (1-prob_lig1)*(1-prob_lig2) + (prob_lig1*prob_lig2);
            }

            tmp = tmp -> prox;
        }while(tmp!=NULL);


        for(int i = 0 ;i< size_aux_vet; i++)
            if(aux_vet[i] != -1)
                aux_vet[i] =  aux_vet[i]*(1 - aux_vet[i]);


        double sum_node_transition_activity_factor = 0;

        for(int i = 0 ;i< size_aux_vet; i++)
        {
            if(aux_vet[i] != -1)
                sum_node_transition_activity_factor = sum_node_transition_activity_factor + aux_vet[i];
        }

        switching_component_of_power = 0;

        switching_component_of_power = sum_node_transition_activity_factor*Capa_Load*freq*Vcc*Vcc;

        free(aux_vet);
        aux_vet = NULL;

        g_tmp = NULL;
        while(g_lst != NULL)
            remover_R(&g_lst);
        free(g_lst);
        g_lst = NULL;

        return switching_component_of_power;
    }
}

void Nu_Gates(Gene ind[NCOL])
{
    int u = 0;
    double num_portas = 0;
    n_gates = num_portas;

    lista *lst = NULL;
    lista *tmp = NULL;

    for(int i=0;i<SAIDAS;i++)
    {
        if(ind[0].output[i] < Entradas)
            continue;

        u = ind[0].output[i] - Entradas;

        inserir(&lst,u);
        tmp=lst;

        do
        {
            if(ind[(int)tmp->val].prox1 != NULL)
                inserir(&lst,ind[(int)tmp->val].lig1 - Entradas);

            if(ind[(int)tmp->val].prox2 != NULL)
                inserir(&lst,ind[(int)tmp->val].lig2 - Entradas);

            tmp = tmp -> prox;
        }while(tmp!=NULL);
    }

    tmp=lst;
    if(tmp == NULL)
    {
        n_gates = 0;
        return;
    }
    ///Number of gates
    do
    {
        tmp = tmp -> prox;
        num_portas++;
    }while(tmp!=NULL);

    tmp = NULL;

    while(lst != NULL)
        remover_R(&lst);
    free(lst);
    lst = NULL;

    n_gates = num_portas;
}

void fitness(Gene **ind,int **tabela,int linhaTabela, int colunaTabela,int size_pop)
{
    for(int i=0;i<size_pop;i++)
    {
        if(ind[i][0].avaliar == 1)
        {
            ind[i][0].delay = delay(ind[i]);
            ind[i][0].power = power(ind[i]);
            ind[i][0].error = dist_hamming(ind[i],tabela,linhaTabela,colunaTabela);
        }
    }

}

void inserir(lista** lst, int dado)
{
    lista *aux = (lista *)malloc(sizeof(lista));
    aux->val = dado;
    aux -> prox = NULL;
    lista *tmp = NULL;

    if(*lst == NULL)
    {
        *lst = aux;
        aux = NULL;
        return;
    }
    else if(aux->val == (*lst)->val)
    {
        free(aux);
        aux = NULL;
        return;
    }

    else if(aux->val > (*lst)->val)
    {
        aux->prox = *lst;
        *lst = aux;
        aux = NULL;
        return;
    }
    else
    {
        tmp = *lst;
        while(tmp->prox != NULL)
        {
            if(aux->val == tmp->prox->val)
            {
                free(aux);
                aux = NULL;
                return;
            }
            else if(aux->val > tmp->prox->val)
            {
                aux->prox = tmp->prox;
                tmp->prox = aux;
                aux = NULL;
                return;
            }
            tmp = tmp -> prox;
        }
        tmp->prox = aux;
        aux = NULL;
        return;
    }
}

void remover_R(lista ** lst)
{
    lista *ptr;
    if(*lst == NULL)
         return;
    else
    {
        ptr = *lst;
        *lst = (*lst)->prox;
        free(ptr);
        return;
    }
}

void remover(lista ** lst, int dado)
{
    lista *ptr, *antes;
    if(*lst == NULL)
        return;
  else
  {
      ptr = *lst;
      antes = *lst;
      while (ptr !=NULL)
      {
         if (ptr->val == dado)
         {
            if (ptr == *lst)
            {
               *lst = (*lst)->prox;
               free(ptr);
               return;
            }
            else
            {
              antes->prox = ptr->prox;
              free(ptr);
              return;
            }
         }
         else
         {
            antes = ptr;
            ptr = ptr->prox;
         }
      }
      return;
  }

}

void insere_ord(lista** lst, int val)
{
    lista *aux = malloc(sizeof(lista));
    aux->val = val;
    aux -> prox = NULL;
    lista *tmp = NULL;

    if(*lst == NULL)
    {
        *lst = aux;
        aux = NULL;
        return;
    }
    else if(aux->val == (*lst)->val)
    {
        free(aux);
        aux = NULL;
        return;
    }
    else if(aux->val < (*lst)->val)
    {
        aux->prox = *lst;
        *lst = aux;
        aux = NULL;
        return;
    }
    else
    {
        tmp = *lst;
        while(tmp->prox != NULL)
        {
            if(aux->val == tmp->prox->val)
            {
                free(aux);
                aux = NULL;
                return;
            }
            else if(aux->val < tmp->prox->val)
            {
                aux->prox = tmp->prox;
                tmp->prox = aux;
                aux = NULL;
                return;
            }
            tmp = tmp -> prox;
        }
        tmp->prox = aux;
        aux = NULL;
        return;
    }
}

void remover_R_Nearest(lista_Euc_nearest ** lst)
{
    lista_Euc_nearest *ptr;
    if(*lst == NULL)
         return;
    else
    {
        ptr = *lst;
        *lst = (*lst)->prox;
        free(ptr);
        return;
    }
}

void insere_decrescebte_nearest(lista_Euc_nearest** lst, double error, int ind) ///jafoi
{
    lista_Euc_nearest *atual, *novo, *anterior;

    novo = (lista_Euc_nearest *) malloc(sizeof(lista_Euc_nearest));

    atual = *lst;
    anterior = NULL;

    novo->val = error;
    novo->nearest = ind;

    if(atual == NULL)
    {
        novo->prox = NULL;
        *lst = novo;
    }
    else
    {
        while(atual != NULL && atual->val > error)
        {
            anterior = atual;
            atual = atual->prox;
        }

        novo->prox = atual;

        if(anterior == NULL)
        {
            *lst = novo;
        } else
        {
            anterior->prox = novo;
        }
    }
}

int busca(lista** lst, int val)
{
    if(*lst == NULL)
        return 0;
    else
    {
        lista *tmp = *lst;
        while(tmp != NULL)
        {
            if(tmp->val == val)
                return 1;
            tmp = tmp -> prox;
        }
        return 0;
    }
}
