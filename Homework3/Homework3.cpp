#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <math.h>
#include <vector>

#define LMAX 1000 //lungime sir de biti
#define CMAX 290 //maxim de componente

#define ITERATII 30
#define DIM_POPULATIE 200
#define GENERATII_MAX 2000

using namespace std;
using namespace std::chrono;

ifstream fin("sursa.in");

float INCRUCISARE_PROB;
float MUTATIE_PROB;

default_random_engine generator;

float sel[400];
int vecin[290];
int solutiiSA[35], solutiiGA[35];
int ciclu[290];
float coordx[290], coordy[290];
int G[290][290];
int componente;
const int pre = 5; //precizia

struct membru {
    int n;
    int cost;
    int gene[290];
    float fitness;
    float prob_random;
}population[300], children1, children2, p1, p2;

membru vechia_populatie[300];
float valoare_minima, minim_vector[CMAX];
float interval_start, interval_end;
float valoarea_globala;

void Crossover(int&);


float simulateEvolution();
void Evaluate(int);
int MinScore(int);
float fitness(int, int);
void Select(int&);
float GenerareRandomFloat(float, float);
int GenerareRandomInt(int);
void generatePopulation(int&);
void SmallCrossOver(membru, membru);

void Mutate(int);
void MutateM(membru);
void MutatePop(membru, int, int);
void MutateGreedy(membru);
int RandomInt(int);
void StrongMutation(membru);

void init(int, int[300]);

long long CurrentTime();

int pivotFitness(int, int);
int pivotProbabilitate(int, int);
void SortarePopulatieDupaFitness(int, int);
void SortarePopulatieDupaProbabilitate(int, int);

int score(int[290], int, int[290][290]);
void decodificareGene(int[290], int);
void randomMutation(membru);
float saFitness(int);
float simulatedAnnealing();
void readGraph();

int minim_curent;

int main()
{
    readGraph();
    
    /*
    cout << "Aici " << componente << "\n";
    for (int i = 0; i < componente; i++) {
        for (int j = 0; j < componente; j++) {
            cout << G[i][j] << " ";
        }
        cout << "\n";
    }
    */
    

    /*
    int t1 = CurrentTime();
    float vmin = simulatedAnnealing();
    int t2 = CurrentTime();
    cout <<"Minimul este: "<< vmin<<"\n";
    float sum = 0;
    for (int i = 1; i <= 30; i++)
        sum += solutiiSA[i];
    float ma = sum / 30;
    cout << "Media este " << ma<<"\n";
    float sumDS = 0;
    for (int i = 1; i <= 30; i++)
        sumDS += (solutiiSA[i] - ma) * (solutiiSA[i] - ma);
    sumDS /= 29;
    sumDS = sqrt(sumDS);
    cout << "Deviatia Standard este " << sumDS << "\n";
    cout << "Timp(milisec) :" << t2 - t1 << "\n";
    int t1 = CurrentTime();
    float vming = simulateEvolution();
    float val;
    solutiiGA[1] = vming;
    printf("GA:  1  cost: ");
    cout << vming << "\n";
    for (int i = 1; i < 30; i++)
    {
        val = simulateEvolution();
        solutiiGA[i + 1] = val;
        printf("GA:  %d  cost: ", i+1);
        cout << val << "\n";
        if (val < vming) vming = val;
    }
    int t2 = CurrentTime();
    cout << "\nMinim: " << vming;
    float sum = 0;
    for (int i = 1; i <= 30; i++)
        sum += solutiiGA[i];
    float ma = sum / 30;
    cout << "Media este " << ma << "\n";
    float sumDS = 0;
    for (int i = 1; i <= 30; i++)
        sumDS += (solutiiGA[i] - ma) * (solutiiGA[i] - ma);
    sumDS /= 29;
    sumDS = sqrt(sumDS);
    cout << "Deviatia Standard este " << sumDS << "\n";
    cout << "Timp(milisec) :" << t2 - t1 << "\n";
    */
}

void readGraph()
{
    int j;
    fin >> componente;
    for (int i = 0; i < componente; i++)
    {
        fin >> j;
        fin >> coordx[j];
        fin >> coordy[j];
    }
    for (int i = 0; i < componente; i++) {
        for (int j = i+1; j < componente; j++) 
        {
            G[i][j] = sqrt(pow(coordx[j] - coordx[i], 2) + pow(coordy[j] - coordy[i], 2));
            G[j][i] = G[i][j];
        }
    }
}

int GenerareRandomInt(int b)
{
    uniform_int_distribution<int> distribution(0, b);
    int nr = distribution(generator);
    return nr;
}

void init(int n, int gene[290])
{
    for (int i = 0; i < n; i++) 
    {
        gene[i] = GenerareRandomInt(n - 1 - i);
    }
}

int score(int gene[290], int n, int G[290][290]) 
{
    int s = 0;
    gene[n - 1] = 0;
    decodificareGene(gene, n);
    for (int i = 0; i < n - 1; i++) {
        int c = G[ciclu[i]][ciclu[i + 1]];
        s += c;
    }
    int x = ciclu[n - 1];
    int y = ciclu[0];
    s += G[x][y];
    return s;
}

void decodificareGene(int gene[290], int n) 
{
    vector<int> vlist;
    for (int i = 0; i < n; i++) 
    {
        vlist.push_back(i);
    }
    for (int i = 0; i < n; i++) 
    {
        ciclu[i] = vlist[gene[i]];
        vlist.erase(vlist.begin() + gene[i]);
    }
}


void randomMutation(membru m) 
{
    int r = GenerareRandomInt(m.n - 1);

    m.gene[r]++;
    int mod = m.n - 1 - r;
    if (mod == 0) {
        m.gene[r] = 0;;
    }
    else {
        m.gene[r] %= mod;
    }
}

void Cautare_vecin(int vect[290], int n, int poz)
{
    int vmin = 2000000000, pozz=0;
    for (int i = 0; i < componente-1-poz; i++)
    {
        vect[poz] = i;
        int val = score(vect, componente, G);
        if (val < vmin) {
            vmin = val;
            pozz = i;
        }
    }
    vect[poz] = pozz;
}

float simulatedAnnealing() 
{
    int solutia_curenta[290];
    srand(time(0));
    int minCost = 2000000000;
    int valoarea_globala, valoarea_curenta, valoare_vecin;
    for (int i = 1; i <= ITERATII; i++) 
    {
        init(componente, solutia_curenta);
        valoarea_globala = score(solutia_curenta, componente, G);
        for (int k = 0; k < 1000; k++)
        {
            cout << k << " ";
            float t = 100;
            init(componente, solutia_curenta);
            valoarea_curenta=score(solutia_curenta, componente, G);
            while (t > 10e-8) 
            {
                for (int j = 0; j < componente; j++)
                    vecin[j] = solutia_curenta[j];
                int rand_poz = GenerareRandomInt(componente - 1);
                Cautare_vecin(vecin, componente, rand_poz);
                valoare_vecin = score(vecin, componente, G);
                float x = GenerareRandomFloat(0, 1 - 0.00001);
                if (valoare_vecin < valoarea_curenta) {
                    for (int j = 0; j < componente; j++)
                        solutia_curenta[j] = vecin[j]; //punem valuarea vecinului in valoarea curenta
                    valoarea_curenta = valoare_vecin;
                }
                else if (x < exp(-abs(valoarea_curenta - valoare_vecin) / t)) {
                    for (int j = 0; j < componente; j++)
                        solutia_curenta[j] = vecin[j];
                    valoarea_curenta = valoare_vecin;
                }
                t *= 0.99;
            }
            if (valoarea_curenta < valoarea_globala)
                valoarea_globala = valoarea_curenta;
        }
        solutiiSA[i] = valoarea_globala;
        minCost = min(minCost, valoarea_globala);
        printf("SA:  %d  cost: %d\n", i, valoarea_globala);
    }
    return minCost;
}

float simulateEvolution()
{
    MUTATIE_PROB = 1.0 / componente;
    INCRUCISARE_PROB = 0.3;
    int n = DIM_POPULATIE;
    srand(time(0));
    generatePopulation(n);
    Evaluate(n);
    for (int i = 1; i <= GENERATII_MAX; i++)
    {
        Select(n);
        Mutate(n);
        Crossover(n);
        Evaluate(n);
        MUTATIE_PROB *= 1.0005;
    }
    SortarePopulatieDupaFitness(0, n - 1);
    return population[0].cost;
}

void generatePopulation(int& n)
{
    n = DIM_POPULATIE;
    for (int i = 0; i < DIM_POPULATIE; i++)
    {
        population[i].n = componente;
        init(componente, population[i].gene);
    }
}


void Evaluate(int n)
{
    for (int i = 0; i < n; i++) {
        population[i].cost = score(population[i].gene, componente, G);
    }

    int minScore = MinScore(n);

    for (int i = 0; i < n; i++)
    {
        population[i].fitness = fitness(population[i].cost, minScore);
    }
}

int MinScore(int n) 
{
    int minScore = 2000000000;
    for (int i = 0; i < n; i++) 
    {
        minScore = min(minScore, population[i].cost);
    }
    return minScore;
}



float fitness(int sco, int minScore)
{
    float score = sco - minScore + 10;
    score = 1.0 / score + 2;
    float f = (pow(5.0, score) - 25) * 100;
    return f;
}

int pivotFitness(int st, int dr)
{
    int piv = population[dr].fitness;
    int i = st - 1;
    for (int j = st; j < dr; ++j)
        if (population[j].fitness < piv)
        {
            i++;
            membru aux = population[i];
            population[i] = population[j];
            population[j] = aux;
        }
    membru aux = population[i + 1];
    population[i + 1] = population[dr];
    population[dr] = aux;
    return i + 1;
}

void SortarePopulatieDupaFitness(int st, int dr)
{
    if (st < dr)
    {
        int p = pivotFitness(st, dr);
        SortarePopulatieDupaFitness(st, p - 1);
        SortarePopulatieDupaFitness(p + 1, dr);
    }
}

float GenerareRandomFloat(float a, float b)
{
    uniform_real_distribution<> distribution3(a, b); //generam o probabilitate
    float x = distribution3(generator);
    return x;
}

void MutateM(membru m) {

    for (int i = 0; i < m.n; i++) {
        float r = GenerareRandomFloat(0, 1);
        if (r < MUTATIE_PROB) {
            m.gene[i]++;
            int mod = m.n - 1 - i;
            if (mod == 0) {
                m.gene[i] = 0;;
            }
            else {
                m.gene[i] %= mod;
            }
        }
    }
}




void Select(int& n)
{
    float TOTAL_SCORE = 0;
    for (int i = 0; i < n; i++) 
    {
        TOTAL_SCORE += population[i].fitness;
    }
    sel[0] = 0;
    for (int i = 1; i <= n; i++) 
    {
        sel[i] = population[i - 1].fitness / TOTAL_SCORE;
    }
    for (int i = 1; i <= n; i++) 
    {
        sel[i] += sel[i - 1];
    }
    int oldn = n;
    for (int i = 0; i < oldn; i++)
    {
        vechia_populatie[i] = population[i];
    }
    for (int i = 0; i < DIM_POPULATIE; i++) {
        float r = GenerareRandomFloat(0.000001, 1);
        int s;
        for (s = 0; s < oldn - 1; s++) {
            if (sel[s] < r && r <= sel[s + 1]) {
                break;
            }
        }
        population[i] = vechia_populatie[s];
    }
    n = DIM_POPULATIE;
}

int RandomInt(int b)
{
    srand(time(0));
    int nr = rand() % (b + 1);
    while (nr == 0)
    {
        nr = rand() % (b + 1);
    }
    return nr;
}

void Mutate(int n)
{
    for (int j = 0; j < n; j++)
    {
        MutateM((population[j]));
    }
}

void SmallCrossOver(membru a, membru b) 
{
    int pivot = RandomInt(componente - 2);
    for (int i = 0; i <= pivot; i++)
    {
        children1.gene[i] = a.gene[i];
        children2.gene[i] = b.gene[i];
    }
    for (int i = pivot + 1; i < componente; i++)
    {
        children1.gene[i] = b.gene[i];
        children2.gene[i] = a.gene[i];
    }
}


void Crossover(int& n)
{
    for (int i = 0; i < n; i++) 
    {
        population[i].prob_random = GenerareRandomFloat(0, 1);
    }
    SortarePopulatieDupaProbabilitate(0, n - 1);
    int taietura = 0;
    while (population[taietura].prob_random < INCRUCISARE_PROB && taietura < n) 
    {
        taietura++;
    }
    taietura--;
    if (taietura % 2 == 1)
    {
        int nr = RandomInt(10);
        if (nr % 2 == 1)
            taietura++;
        else
            taietura--;
    }
    for (int i = 0; i < taietura; i += 2)
    {
        SmallCrossOver(population[i], population[i + 1]);
        population[n++] = children1;
        population[n++] = children2;
    }
}

long long CurrentTime() 
{
    milliseconds ms = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
        );
    return ms.count();
}




int pivotProbabilitate(int st, int dr)
{
    int piv = population[dr].prob_random;
    int i = st - 1;
    for (int j = st; j < dr; ++j)
        if (population[j].prob_random < piv)
        {
            i++;
            membru aux = population[i];
            population[i] = population[j];
            population[j] = aux;
        }
    membru aux = population[i + 1];
    population[i + 1] = population[dr];
    population[dr] = aux;
    return i + 1;
}

void SortarePopulatieDupaProbabilitate(int st, int dr)
{
    if (st < dr)
    {
        int p = pivotProbabilitate(st, dr);
        SortarePopulatieDupaProbabilitate(st, p - 1);
        SortarePopulatieDupaProbabilitate(p + 1, dr);
    }
}