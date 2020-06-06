#include <math.h>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <time.h>

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
#define PArent_POP 50 ///Number of parents
#define NPOPULACAO 2*PArent_POP  ///Parents + offspring

///Seed
#define seed_do_circuito 31

typedef struct Gene
{
    int lig1;
    int lig2;
    int porta; ///Gate
    int checked;
    int logical_out;
    int frente;
    int output[SAIDAS];
    int avaliar;
    double distance;
    double atraso;
    double error;
    double power;
    double delay;
    struct Gene *prox1, *prox2;
}Gene;

typedef struct lista
{
    double val;
    struct lista * prox;
}lista , *vetor[NPOPULACAO];

typedef struct lista_Euc_nearest
{
    int nearest;
    double val;
    struct lista_Euc_nearest * prox;
}lista_Euc_nearest;

typedef struct lista_Euc
{
    int anterior;
    int posterior;
    double val;
    struct lista_Euc * prox;
}lista_Euc;

void remover_R(lista ** lst);
void inserir(lista** lst, double val);
void exibe(lista**lst);
void remover(lista ** lst, double val);
void insere_ord(lista** lst, double val);
void remover_R_EUC(lista_Euc ** lst);
void print_results(Gene **Rt, int linhaTabela, int colunaTabela);
void insere_ord_Crowded_Comparison_Operator(lista_Euc** lst, double distance, int frente, int ind);

int busca(lista** lst, int val);
void FNDS(Gene **Rt, lista** Front); /// Fast non dominated sort
void CrowDis(Gene **ind,lista **Front, double *Vet_CRD);///Crowd Distance
void Crowded_Comparison_Operator_Sorting(Gene **Pt,int actual_parent);
void copyPt_to_Rt(Gene** Rt,Gene** Pt);

int saida_porta(int r1, int r2, int porta);
int retornaValor(Gene Rt[NCOL], int linha[], int saida);
double dist_hamming(Gene ind[NCOL], int **tabela,int linhaTabela, int colunaTabela);
double delay(Gene ind[NCOL]);
double power(Gene ind[NCOL]);
void fitness(Gene **ind,int **tabela,int linhaTabela, int colunaTabela);
void Nu_Gates(Gene ind[NCOL]);
void create_Pt(Gene **Pt, Gene **Rt, int * actual_parent,lista** Front);

void criaInd(Gene ind[NCOL],int portas[], int saidas[]);
void mutacao(Gene** idividual_mutated,int portas[PORTAS], int comeco, int fim);

int numAleatorio(int a, int b);
unsigned long long int bin_to_dec(unsigned long long int bin);
void sorteiaSaidas(int saidas[]);
void leTxt(int **tabela,int nL,int nC);
void imprimeTabela(int **vet, int nL,int nC);
void selecionaPortas(int vet[PORTAS]);

double cont_evaluation,percent,n_gates;
double Max_EVALUATION = 1000;

int main()
{
    percent = pow(2,Entradas)*SAIDAS;
    int saidas[SAIDAS]; /// Vetor de Saidas de cada indivíduo
    int portas[PORTAS]; /// Vetor que armazena as portas que podem fazer parte do circuito
    int linhaTabela,colunaTabela;// actual_parent, current_front;

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

    int **tabela = NULL; /// Truth Table
    tabela = (int**)calloc(linhaTabela,sizeof(int*));

    for(int j = 0; j<linhaTabela;j++)
        tabela[j] = (int*)calloc(colunaTabela,sizeof(int));

    leTxt(tabela,linhaTabela,colunaTabela);

    selecionaPortas(portas);

    Gene ** Rt = NULL;
    Gene ** Pt = NULL;

    Rt = (Gene**)calloc(NPOPULACAO,sizeof(Gene*));
    Pt = (Gene**)calloc(NPOPULACAO,sizeof(Gene*));


    for(int j=0; j< NPOPULACAO; j++)
        Rt[j] = (Gene*)calloc(NCOL,sizeof(Gene));

    for(int j=0; j< NPOPULACAO; j++)
        Pt[j] = (Gene*)calloc(NCOL,sizeof(Gene));

    srand((unsigned)seed_do_circuito);

    printf("seed: %d Max_EVALUATION: %f Parent Pop: %d Populacao: %d\n",seed_do_circuito,Max_EVALUATION,PArent_POP,NPOPULACAO);


    ///Rt = Pt + Qt *******************************************************************************************************************************
    for(int z = 0; z<NPOPULACAO; z++)
    {
        for(int j=0; j< NCOL; j++)
        {
            Rt[z][j].lig1 = -1;
            Rt[z][j].lig2 = -1;
            Rt[z][j].checked = -1;
            Rt[z][j].logical_out = -1;
            Rt[z][j].porta = -1;
            Rt[z][j].error = -1;
            Rt[z][j].avaliar = 1;
            Rt[z][j].delay = -1;
            Rt[z][j].power = -1;
            Rt[z][j].atraso = -1;
            Rt[z][j].frente = 2*NPOPULACAO;
            Rt[z][j].distance = -1;
            Rt[z][j].prox1 = NULL;
            Rt[z][j].prox2 = NULL;
        }

        for(int j=0; j< SAIDAS; j++)
        {
            Rt[z][0].output[j] = 0;
        }
    }

    ///***Creates Pop**************************************************************************************************

    for(int j = 0; j< NPOPULACAO; j++)
        criaInd(Rt[j],portas,saidas);

    ///Parents**************************************************************************************************************************

    for(int z = 0; z<NPOPULACAO; z++)
    {
        for(int j=0; j< NCOL; j++)
        {
            Pt[z][j].lig1 = -1;
            Pt[z][j].lig2 = -1;
            Pt[z][j].porta = -1;
            Pt[z][j].checked = -1;
            Pt[z][j].logical_out = -1;
            Pt[z][j].atraso = -1;
            Pt[z][j].power = -1;
            Pt[z][j].avaliar = -1;
            Pt[z][j].error = -1;
            Pt[z][j].delay = -1;
            Pt[z][j].frente = -1;
            Pt[z][j].distance = -1;
            Pt[z][j].prox1 = NULL;
            Pt[z][j].prox2 = NULL;
        }
        for(int j=0; j< SAIDAS; j++)
        {
            Pt[z][0].output[j] = 0;
        }
    }

    ///**************************************************************************************************************************
    cont_evaluation = 0;
    vetor Front;

    while(cont_evaluation < Max_EVALUATION)
    {
        fitness(Rt,tabela,linhaTabela,colunaTabela);
        cont_evaluation = cont_evaluation + NPOPULACAO;

        for(int z=0; z<NPOPULACAO;z++)
            Front[z] = NULL;

        FNDS(Rt,Front);

        int actual_parent = 0;
        create_Pt(Pt,Rt,&actual_parent,Front);

        for(int z = 0;z< NPOPULACAO; z++)
        {
            while(Front[z] != NULL)
                remover_R(&Front[z]);
            free(Front[z]);
        }

        if(actual_parent > PArent_POP)
            Crowded_Comparison_Operator_Sorting(Pt,actual_parent);

        if(cont_evaluation >= Max_EVALUATION)
        {
            print_results(Pt,linhaTabela,colunaTabela);
            break;
        }

        copyPt_to_Rt(Rt,Pt);
        mutacao(Rt,portas,PArent_POP,NPOPULACAO);
    }

    for(int j = 0; j<NPOPULACAO;j++)
        free(Rt[j]);
    free(Rt);

    for(int j = 0; j<NPOPULACAO;j++)
        free(Pt[j]);
    free(Pt);

    for(int j = 0; j<linhaTabela;j++)
        free(tabela[j]);
    free(tabela);

    return 0;
}

void create_Pt(Gene **Pt, Gene **Rt, int * actual_parent,lista** Front)
{
    int current_front = 0;

    while(*actual_parent<PArent_POP)
    {
        lista *tmp = Front[current_front];

        double Vet_CRD[NPOPULACAO+1];

        for(int j = 0;j<NPOPULACAO+1;j++)
            Vet_CRD[j] = 0;

        CrowDis(Rt,&Front[current_front],Vet_CRD);

        for(int z = 0; z < (int)Vet_CRD[0];z++)
        {
            for(int j=0; j< NCOL; j++)
            {
                Pt[*actual_parent][j].lig1 = Rt[(int)tmp->val][j].lig1;
                Pt[*actual_parent][j].lig2 = Rt[(int)tmp->val][j].lig2;
                Pt[*actual_parent][j].checked = Rt[(int)tmp->val][j].checked;
                Pt[*actual_parent][j].logical_out = Rt[(int)tmp->val][j].logical_out;
                Pt[*actual_parent][j].porta = Rt[(int)tmp->val][j].porta;
                Pt[*actual_parent][j].atraso = Rt[(int)tmp->val][j].atraso;
                Pt[*actual_parent][j].power = Rt[(int)tmp->val][j].power;
                Pt[*actual_parent][j].error = Rt[(int)tmp->val][j].error;
                Pt[*actual_parent][j].delay = Rt[(int)tmp->val][j].delay;
                Pt[*actual_parent][j].frente = current_front;
                Pt[*actual_parent][j].distance = Vet_CRD[(int)tmp->val+1];
                if(Pt[*actual_parent][j].lig1 < Entradas)
                    Pt[*actual_parent][j].prox1 = NULL;
                else
                    Pt[*actual_parent][j].prox1 = &Pt[*actual_parent][Pt[*actual_parent][j].lig1 - Entradas];

                if(Pt[*actual_parent][j].lig2 < Entradas)
                    Pt[*actual_parent][j].prox2 = NULL;
                else
                    Pt[*actual_parent][j].prox2 = &Pt[*actual_parent][Pt[*actual_parent][j].lig2 - Entradas];

            }

            for(int j=0; j< SAIDAS; j++)
            {
                Pt[*actual_parent][0].output[j] = Rt[(int)tmp->val][0].output[j];
            }

            *actual_parent = *actual_parent + 1;
            tmp = tmp->prox;
        }

        current_front++;
        tmp = NULL;
    }
}

void Crowded_Comparison_Operator_Sorting(Gene **Pt,int actual_parent)
{
    lista_Euc *lst = NULL;
    lista_Euc *aux = NULL;

    for(int y = 0; y < actual_parent; y++)
        insere_ord_Crowded_Comparison_Operator(&lst,Pt[y][0].distance,Pt[y][0].frente,y);

    Gene ** Temporary = NULL;

    Temporary = (Gene**)calloc(actual_parent,sizeof(Gene*));

    for(int j=0; j< actual_parent; j++)
        Temporary[j] = (Gene*)calloc(NCOL,sizeof(Gene));

    ///*****************************************************************************************************************************
    for(int z = 0; z<actual_parent; z++)
    {
        for(int j=0; j< NCOL; j++)
        {
            Temporary[z][j].lig1 = -1;
            Temporary[z][j].lig2 = -1;
            Temporary[z][j].checked = -1;
            Temporary[z][j].logical_out = -1;
            Temporary[z][j].porta = -1;
            Temporary[z][j].atraso = -1;
            Temporary[z][j].power = -1;
            Temporary[z][j].error = -1;
            Temporary[z][j].delay = -1;
            Temporary[z][j].frente = -1;
            Temporary[z][j].distance = -1;
            Temporary[z][j].prox1 = NULL;
            Temporary[z][j].prox2 = NULL;
        }
    }
    for(int z = 0; z<actual_parent; z++)
    {
        for(int j=0; j< SAIDAS; j++)
        {
            Temporary[z][0].output[j] = 0;
        }
    }

    ///*****************************************************************************************************************************
    for(int z = 0; z<actual_parent; z++)
    {
        for(int j=0; j< NCOL; j++)
        {
            Temporary[z][j].lig1 = Pt[z][j].lig1;
            Temporary[z][j].lig2 = Pt[z][j].lig2;
            Temporary[z][j].checked = Pt[z][j].checked;
            Temporary[z][j].logical_out = Pt[z][j].logical_out;
            Temporary[z][j].porta = Pt[z][j].porta;
            Temporary[z][j].atraso = Pt[z][j].atraso;
            Temporary[z][j].power = Pt[z][j].power;
            Temporary[z][j].error = Pt[z][j].error;
            Temporary[z][j].delay = Pt[z][j].delay;
            Temporary[z][j].frente = Pt[z][j].frente;
            Temporary[z][j].distance = Pt[z][j].distance;

            if(Temporary[z][j].lig1 < Entradas)
                Temporary[z][j].prox1 = NULL;
            else
                Temporary[z][j].prox1 = &Temporary[z][Temporary[z][j].lig1 - Entradas];

            if(Temporary[z][j].lig2 < Entradas)
                Temporary[z][j].prox2 = NULL;
            else
                Temporary[z][j].prox2 = &Temporary[z][Temporary[z][j].lig2 - Entradas];
        }
    }
    for(int z = 0; z<actual_parent; z++)
    {
        for(int j=0; j< SAIDAS; j++)
        {
            Temporary[z][0].output[j] = Pt[z][0].output[j];
        }
    }

    ///*****************************************************************************************************************************

    aux = lst;

    for(int new_index = 0; new_index < PArent_POP; new_index++) ///
    {
        for(int j=0; j< NCOL; j++)
        {
            Pt[new_index][j].lig1 = Temporary[aux->posterior][j].lig1;
            Pt[new_index][j].lig2 = Temporary[aux->posterior][j].lig2;
            Pt[new_index][j].checked = Temporary[aux->posterior][j].checked;
            Pt[new_index][j].logical_out = Temporary[aux->posterior][j].logical_out;
            Pt[new_index][j].porta = Temporary[aux->posterior][j].porta;
            Pt[new_index][j].error = Temporary[aux->posterior][j].error;
            Pt[new_index][j].delay = Temporary[aux->posterior][j].delay;
            Pt[new_index][j].power = Temporary[aux->posterior][j].power;
            Pt[new_index][j].frente = Temporary[aux->posterior][j].frente;
            Pt[new_index][j].distance = Temporary[aux->posterior][j].distance;
            Pt[new_index][j].atraso = Temporary[aux->posterior][j].atraso;

            if(Pt[new_index][j].lig1 < Entradas)
                Pt[new_index][j].prox1 = NULL;
            else
                Pt[new_index][j].prox1 = &Pt[new_index][Pt[new_index][j].lig1 - Entradas];

            if(Pt[new_index][j].lig2 < Entradas)
                Pt[new_index][j].prox2 = NULL;
            else
                Pt[new_index][j].prox2 = &Pt[new_index][Pt[new_index][j].lig2 - Entradas];
        }
        aux = aux->prox;
    }

    aux = lst;

    for(int new_index = 0; new_index < PArent_POP; new_index++)
    {
        for(int j=0; j< SAIDAS; j++)
        {
            Pt[new_index][0].output[j] = Temporary[aux->posterior][0].output[j];
        }
        aux = aux->prox;
    }

    aux = NULL;
    while(lst != NULL)
        remover_R_EUC(&lst);

    for(int j = 0; j<actual_parent;j++)
        free(Temporary[j]);
    free(Temporary);
}

void copyPt_to_Rt(Gene** Rt,Gene** Pt)
{
    for(int i = 0; i<PArent_POP;i++)
    {
        for(int z = 0; z<NCOL;z++)
        {
            Rt[i][z].lig1 =     Pt[i][z].lig1;
            Rt[i][z].lig2 =     Pt[i][z].lig2;
            Rt[i][z].checked =    Pt[i][z].checked;
            Rt[i][z].logical_out =    Pt[i][z].logical_out;
            Rt[i][z].porta =    Pt[i][z].porta;
            Rt[i][z].avaliar =     0;
            Rt[i][z].atraso =   Pt[i][z].atraso;
            Rt[i][z].delay =    Pt[i][z].delay;
            Rt[i][z].error =    Pt[i][z].error;
            Rt[i][z].power =    Pt[i][z].power;
            Rt[i][z].frente =   Pt[i][z].frente;
            Rt[i][z].distance = Pt[i][z].distance;
            if(Rt[i][z].lig1 < Entradas)
                Rt[i][z].prox1 = NULL;
            else
                Rt[i][z].prox1 = &Rt[i][Rt[i][z].lig1 - Entradas];

            if(Rt[i][z].lig2 < Entradas)
                Rt[i][z].prox2 = NULL;
            else
                Rt[i][z].prox2 = &Rt[i][Rt[i][z].lig2 - Entradas];
        }
    }
    for(int z = 0; z<PArent_POP; z++)
    {
        for(int j=0; j< SAIDAS; j++)
        {
            Rt[z][0].output[j] = Pt[z][0].output[j];
        }
    }

    ///****************************************************************************************************************************
    for(int i = PArent_POP; i<NPOPULACAO;i++)
    {
        for(int z = 0; z<NCOL;z++)
        {
            Rt[i][z].lig1 =     Pt[i - PArent_POP][z].lig1;
            Rt[i][z].lig2 =     Pt[i - PArent_POP][z].lig2;
            Rt[i][z].checked =     Pt[i - PArent_POP][z].checked;
            Rt[i][z].logical_out =     Pt[i - PArent_POP][z].logical_out;
            Rt[i][z].porta =    Pt[i - PArent_POP][z].porta;
            Rt[i][z].avaliar =     Pt[i - PArent_POP][z].avaliar;
            Rt[i][z].atraso =   Pt[i - PArent_POP][z].atraso;
            Rt[i][z].delay =    Pt[i - PArent_POP][z].delay;
            Rt[i][z].error =    Pt[i - PArent_POP][z].error;
            Rt[i][z].power =    Pt[i - PArent_POP][z].power;
            Rt[i][z].frente =   Pt[i - PArent_POP][z].frente;
            Rt[i][z].distance = Pt[i - PArent_POP][z].distance;
            if(Rt[i][z].lig1 < Entradas)
                Rt[i][z].prox1 = NULL;
            else
                Rt[i][z].prox1 = &Rt[i][Rt[i][z].lig1 - Entradas];

            if(Rt[i][z].lig2 < Entradas)
                Rt[i][z].prox2 = NULL;
            else
                Rt[i][z].prox2 = &Rt[i][Rt[i][z].lig2 - Entradas];
        }
    }
    for(int z = PArent_POP; z<NPOPULACAO; z++)
    {
        for(int j=0; j< SAIDAS; j++)
        {
            Rt[z][0].output[j] = Pt[z - PArent_POP][0].output[j];
        }
    }
}


void print_results(Gene **Rt, int linhaTabela, int colunaTabela)
{
    printf("\n");
    printf("Numero avaliacoes: %f\n",cont_evaluation);
    printf("\n");


    printf("Num Ind: %d\n",PArent_POP);
    printf("\nD_H  | Delay | Potencia\n\n");

    for(int i = 0; i < PArent_POP; i++)
        printf("%f %f %f\n",Rt[i][0].error/percent,Rt[i][0].delay,Rt[i][0].power);
    printf("\n");

    ///*------------------------------------------------------------------------------------------------------

    printf("-------------------------------------------------------------------------------------------------------------------------\n");
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

        ind[j].frente = -1;
        ind[j].distance = -1;
        ind[j].power = -1;
        ind[j].error = -1;
        ind[j].delay = -1;
    }
}

void mutacao(Gene** idividual_mutated,int portas[PORTAS], int comeco, int fim)
{
    int teste = (int)NMUTACOES;

    if(teste < 1)
        teste = 1;

    for(int s = comeco ; s < fim; s++)
    {
        int auxiliar;
        int gen, coluna, aleat;

        idividual_mutated[s][0].frente = -1;
        idividual_mutated[s][0].distance = 1;

        int u = 0;

        lista *lst_nos = NULL;
        lista *tmp_nos = NULL;

        for(int i=0;i<SAIDAS;i++)
        {
            if(idividual_mutated[s][0].output[i] < Entradas)
                continue;

            u = idividual_mutated[s][0].output[i] - Entradas;

            inserir(&lst_nos,u);
            tmp_nos=lst_nos;

            do
            {
                if(idividual_mutated[s][(int)tmp_nos->val].prox1 != NULL)
                    inserir(&lst_nos,idividual_mutated[s][(int)tmp_nos->val].lig1 - Entradas);

                if(idividual_mutated[s][(int)tmp_nos->val].prox2 != NULL)
                    inserir(&lst_nos,idividual_mutated[s][(int)tmp_nos->val].lig2 - Entradas);

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
                    idividual_mutated[s][coluna].lig1 = numAleatorio(0, Entradas - 1) ;
                    if(idividual_mutated[s][coluna].lig1 < Entradas)
                    {
                        idividual_mutated[s][coluna].prox1 = NULL;
                    }
                    else
                    {
                        idividual_mutated[s][coluna].prox1 = &idividual_mutated[s][idividual_mutated[s][coluna].lig1 - Entradas];
                    }
                }
                else
                {
                    if(LB >= coluna || LB == -1)
                    {
                        idividual_mutated[s][coluna].lig1 = numAleatorio(0, (Entradas) + (NLIN*coluna) - 1);
                        if(idividual_mutated[s][coluna].lig1 < Entradas)
                        {
                            idividual_mutated[s][coluna].prox1 = NULL;
                        }
                        else
                        {
                            idividual_mutated[s][coluna].prox1 = &idividual_mutated[s][idividual_mutated[s][coluna].lig1 - Entradas];
                        }
                    }
                    else
                    {
                        if(LB < coluna && LB > 0)
                        {
                            int aleatorio = numAleatorio(0,(Entradas)+ LB*NLIN - 1);
                            if(aleatorio >= 0  && aleatorio < (Entradas-1))
                            {
                                idividual_mutated[s][coluna].lig1 = aleatorio;
                                if(idividual_mutated[s][coluna].lig1 < Entradas)
                                {
                                    idividual_mutated[s][coluna].prox1 = NULL;
                                }
                                else
                                {
                                    idividual_mutated[s][coluna].prox1 = &idividual_mutated[s][idividual_mutated[s][coluna].lig1 - Entradas];
                                }
                            }
                            else
                            {
                                idividual_mutated[s][coluna].lig1 = numAleatorio((Entradas) + (coluna - LB)*NLIN, (Entradas) + (NLIN*coluna) - 1);
                                if(idividual_mutated[s][coluna].lig1 < Entradas)
                                {
                                    idividual_mutated[s][coluna].prox1 = NULL;
                                }
                                else
                                {
                                    idividual_mutated[s][coluna].prox1 = &idividual_mutated[s][idividual_mutated[s][coluna].lig1 - Entradas];
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
                        idividual_mutated[s][coluna].lig2 = numAleatorio(0, Entradas - 1) ;
                        if(idividual_mutated[s][coluna].lig2 < Entradas)
                        {
                            idividual_mutated[s][coluna].prox2 = NULL;
                        }
                        else
                        {
                            idividual_mutated[s][coluna].prox2 = &idividual_mutated[s][idividual_mutated[s][coluna].lig2 - Entradas];
                        }
                    }
                    else
                    {
                        if(LB >= coluna || LB == -1)
                        {
                            idividual_mutated[s][coluna].lig2 = numAleatorio(0, (Entradas) + (NLIN*coluna) - 1);
                            if(idividual_mutated[s][coluna].lig2 < Entradas)
                            {
                                idividual_mutated[s][coluna].prox2 = NULL;
                            }
                            else
                            {
                                idividual_mutated[s][coluna].prox2 = &idividual_mutated[s][idividual_mutated[s][coluna].lig2 - Entradas];
                            }
                        }
                        else
                        {
                            if(LB < coluna && LB > 0)
                            {
                                int aleatorio = numAleatorio(0,(Entradas)+ LB*NLIN - 1);
                                if(aleatorio >= 0  && aleatorio < (Entradas-1))
                                {
                                    idividual_mutated[s][coluna].lig2 = aleatorio;
                                    if(idividual_mutated[s][coluna].lig2 < Entradas)
                                    {
                                        idividual_mutated[s][coluna].prox2 = NULL;
                                    }
                                    else
                                    {
                                        idividual_mutated[s][coluna].prox2 = &idividual_mutated[s][idividual_mutated[s][coluna].lig2 - Entradas];
                                    }
                                }
                                else
                                {
                                    idividual_mutated[s][coluna].lig2 = numAleatorio((Entradas) + (coluna - LB)*NLIN, (Entradas) + (NLIN*coluna) - 1);
                                    if(idividual_mutated[s][coluna].lig2 < Entradas)
                                    {
                                        idividual_mutated[s][coluna].prox2 = NULL;
                                    }
                                    else
                                    {
                                        idividual_mutated[s][coluna].prox2 = &idividual_mutated[s][idividual_mutated[s][coluna].lig2 - Entradas];
                                    }
                                }
                            }
                        }
                    }
                }
                else if(aleat == 2)
                {
                    auxiliar = numAleatorio(0, (PORTAS-1));
                    idividual_mutated[s][coluna].porta = portas[auxiliar];

                    if(idividual_mutated[s][coluna].porta == 0)
                        idividual_mutated[s][coluna].atraso = tdAND;
                    else if(idividual_mutated[s][coluna].porta == 1)
                        idividual_mutated[s][coluna].atraso = tdOR;
                    else if(idividual_mutated[s][coluna].porta == 2)
                        idividual_mutated[s][coluna].atraso = tdNOR;
                    else if(idividual_mutated[s][coluna].porta == 3)
                        idividual_mutated[s][coluna].atraso = tdNAND;
                    else if(idividual_mutated[s][coluna].porta == 4)
                        idividual_mutated[s][coluna].atraso = tdXOR;
                    else if(idividual_mutated[s][coluna].porta == 5)
                        idividual_mutated[s][coluna].atraso = tdXNOR;
                }
                else if(aleat == 3)
                {
                    checa = 1;
                    auxiliar = numAleatorio(0, (SAIDAS-1));
                    idividual_mutated[s][0].output[auxiliar] = numAleatorio(Entradas, (NLIN*NCOL)+(Entradas) - 1);
                }
            }
        }

        idividual_mutated[s][0].avaliar = checa;
        tmp_nos = NULL;

        while(lst_nos != NULL)
            remover_R(&lst_nos);
        free(lst_nos);
        lst_nos = NULL;
    }
}

int saida_porta(int r1, int r2, int porta)
{
    /// Variables for XOR
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

void fitness(Gene **ind,int **tabela,int linhaTabela, int colunaTabela)
{
    for(int i=0;i<NPOPULACAO;i++)
    {
        if(ind[i][0].avaliar == 1)
        {
            ind[i][0].delay = delay(ind[i]);
            ind[i][0].power = power(ind[i]);
            ind[i][0].error = dist_hamming(ind[i],tabela,linhaTabela,colunaTabela);
        }
    }
}

void FNDS(Gene **Rt, lista** Front)
{
    int number_front = 0;
    int front_size = 0;
    int New_front_size = 0;
    int np[NPOPULACAO];
    vetor Sp;

    for(int p=0; p<NPOPULACAO;p++)
         Sp[p] = NULL;

    for(int p=0; p<NPOPULACAO;p++)
        np[p] = 0;

    for(int p=0;p<NPOPULACAO;p++)
    {
        for(int q=0;q<NPOPULACAO;q++)
        {
            if(p != q)
            {
                if(((Rt[p][0].error <= Rt[q][0].error) && (Rt[p][0].power <= Rt[q][0].power) && (Rt[p][0].delay <= Rt[q][0].delay)) && ((Rt[p][0].power < Rt[q][0].power) || (Rt[p][0].error < Rt[q][0].error) || (Rt[p][0].delay < Rt[q][0].delay)))
                    inserir(&Sp[p],q);
                else if (((Rt[q][0].power <= Rt[p][0].power) && (Rt[q][0].error <= Rt[p][0].error) && (Rt[q][0].delay <= Rt[p][0].delay)) && ((Rt[q][0].power < Rt[p][0].power) || (Rt[q][0].error < Rt[p][0].error) || (Rt[q][0].delay < Rt[p][0].delay)))
                    np[p] = np[p] + 1;
            }
        }
        if(np[p] == 0)
        {
            front_size++;
            inserir(&Front[0],p);
        }
    }

    lista *tmp = NULL;
    lista *tmp_f = NULL;

    while(Front[number_front] != NULL)
    {
        tmp_f = Front[number_front];
        for(int p=0;p<front_size;p++)
        {
            tmp = Sp[(int)tmp_f->val];
            tmp_f = tmp_f->prox;

            while(tmp != NULL)
            {
                np[(int)tmp->val] = np[(int)tmp->val] - 1;
                if(np[(int)tmp->val] == 0)
                {
                    inserir(&Front[number_front+1],tmp->val);
                    New_front_size++;
                }
                tmp = tmp -> prox;
            }
        }

        front_size =  New_front_size;
        New_front_size = 0;
        number_front++;
    }

    tmp = NULL;
    tmp_f = NULL;

    for(int p = 0;p< NPOPULACAO; p++)
    {
        while(Sp[p] != NULL)
            remover_R(&Sp[p]);
        free(Sp[p]);
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
                    if(ind[(int)tmp2->val][0].error < ind[(int)tmp_min->val][0].error)
                       tmp_min = tmp2;
                    tmp2=tmp2->prox;
                }
                if(ind[(int)tmp->val][0].error != ind[(int)tmp_min->val][0].error)
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
                Vet_CRD[(int)tmp2->val+1] = DBL_MAX;

                continue;
            }
            else if(number_solutions == 2)
            {
                Vet_CRD[(int)tmp2->val+1] = DBL_MAX;
                Vet_CRD[(int)tmp2->prox->val+1] = DBL_MAX;
                continue;
            }
            else
            {
                tmp = *Front;
                Vet_CRD[(int)tmp->val+1] = DBL_MAX;
                while(tmp->prox != NULL)
                {
                    tmp = tmp -> prox;
                }
                Vet_CRD[(int)tmp->val+1] = DBL_MAX;
                tmp = *Front;
                lista *tmp_before = tmp;
                tmp = tmp -> prox;
                lista *tmp_after = tmp->prox;
                int r;
                for(r = 2;r<number_solutions;r++)
                {
                    Vet_CRD[(int)tmp->val+1] = Vet_CRD[(int)tmp->val+1] + (ind[(int)tmp_after->val][0].error - ind[(int)tmp_before->val][0].error);
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
                    if(ind[(int)tmp2->val][0].delay < ind[(int)tmp_min->val][0].delay)
                       tmp_min = tmp2;
                    tmp2=tmp2->prox;
                }
                if (ind[(int)tmp->val][0].delay != ind[(int)tmp_min->val][0].delay)
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
                Vet_CRD[(int)tmp2->val+1] = DBL_MAX;
                continue;
            }
            else if(number_solutions == 2)
            {
                Vet_CRD[(int)tmp2->val+1] = DBL_MAX;
                Vet_CRD[(int)tmp2->prox->val+1] = DBL_MAX;
                continue;
            }
            else
            {
                tmp = *Front;
                Vet_CRD[(int)tmp->val+1] = DBL_MAX;
                while(tmp->prox != NULL)
                {
                    tmp = tmp -> prox;
                }
                Vet_CRD[(int)tmp->val+1] = DBL_MAX;
                tmp = *Front;
                tmp = tmp -> prox;
                lista *tmp_after = tmp->prox;
                lista *tmp_before = *Front;
                int r;
                for(r = 2;r<number_solutions;r++)
                {
                    Vet_CRD[(int)tmp->val+1] = Vet_CRD[(int)tmp->val+1] + (ind[(int)tmp_after->val][0].delay - ind[(int)tmp_before->val][0].delay);
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
                    if(ind[(int)tmp2->val][0].power < ind[(int)tmp_min->val][0].power)
                       tmp_min = tmp2;
                    tmp2=tmp2->prox;
                }
                if (ind[(int)tmp->val][0].power != ind[(int)tmp_min->val][0].power)
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
                Vet_CRD[(int)tmp2->val+1] = DBL_MAX;
                continue;
            }
            else if(number_solutions == 2)
            {
                Vet_CRD[(int)tmp2->val+1] = DBL_MAX;
                Vet_CRD[(int)tmp2->prox->val+1] = DBL_MAX;
                continue;
            }
            else
            {
                tmp = *Front;
                Vet_CRD[(int)tmp->val+1] = DBL_MAX;
                while(tmp->prox != NULL)
                {
                    tmp = tmp -> prox;
                }
                Vet_CRD[(int)tmp->val+1] = DBL_MAX;
                tmp = *Front;
                tmp = tmp -> prox;
                lista *tmp_after = tmp->prox;
                lista *tmp_before = *Front;
                int r;
                for(r = 2;r<number_solutions;r++)
                {
                    Vet_CRD[(int)tmp->val+1] = Vet_CRD[(int)tmp->val+1] + (ind[(int)tmp_after->val][0].power - ind[(int)tmp_before->val][0].power);
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

void inserir(lista** lst, double dado)
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

void remover(lista ** lst, double dado)
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

void insere_ord(lista** lst, double val)
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

void remover_R_EUC(lista_Euc ** lst)
{
    lista_Euc *ptr;
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

void insere_ord_Crowded_Comparison_Operator(lista_Euc** lst, double distance, int frente, int ind) ///jafoi
{
    lista_Euc *atual, *novo, *anterior;

    novo = (lista_Euc *) malloc(sizeof(lista_Euc));

    atual = *lst;
    anterior = NULL;

    novo->val = distance;
    novo->anterior = frente;
    novo->posterior = ind;

    if(atual == NULL)
    {
        novo->prox = NULL;
        *lst = novo;
    }
    else
    {

        while( (atual != NULL) && ( (atual->anterior < frente) || ( (atual->anterior == frente) && (atual->val >= distance)     )     )   )
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

