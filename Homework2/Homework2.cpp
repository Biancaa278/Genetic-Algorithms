#include <iostream>
#include <random>
#include <chrono>
#include <math.h>

#define LMAX 1000 //lungime sir de biti
#define CMAX 50 //maxim de componente


#define ITERATII 30
#define DIM_POPULATIE 140
#define GENERATII_MAX 1400
#define INCRUCISARE_PROB 0.7
#define MUTATIE_PROB 0.001

using namespace std;
using namespace std::chrono;

default_random_engine generator;

int componente;
const int pre = 5; //precizia

struct membru {
    int gene[LMAX];
    float fitness;
    float prob_random;
}population[250], children1, children2, m;

float valoare_minima, minim_vector[CMAX];
float interval_start, interval_end;
float valoarea_globala;
int ld, L;

void init(float interval_start, float interval_end)
{
    ld = ceil(log2(pow(10, pre) * (interval_end - interval_start)));
    L = componente * ld;
}



float DeJong(float*);
float Schwefel(float*);
float Rastrigin(float*);
float Michalewicz(float*);

void DeJongsProcess();
void RastriginProcess();
void SchwefelProcess();
void MichalewiczProcess();




float simulateEvolution(char);
void Evaluate(int&, char);
void Select(int&);
float GenerareRandomFloat(float, float);
int GenerareRandomInt(int, int);
void GenerareSirDeBiti(int);
void generatePopulation(int&);
void SmallCrossOver(membru, membru, membru, membru);
void CrossOver(int&);
void Mutate(int&, char);
float getFitness(int*, float, char);


float f(int*, char);
void TransformareBitiInCoordonate(float*, int*);


long long CurrentTime();

int pivotFitness(int, int);
int pivotProbabilitate(int, int);
void SortarePopulatieDupaFitness(int, int);
void SortarePopulatieDupaProbabilitate(int, int);


int main()
{
    cin >> componente;


    //DeJongsProcess();
    //SchwefelProcess();
    RastriginProcess();
    //MichalewiczProcess();
}


float simulateEvolution(char ch)
{
    int n = 0;
    generatePopulation(n);
    Evaluate(n, ch);
    for (int g = 1; g <= GENERATII_MAX; g++)
    {
        Select(n);
        Mutate(n, ch);
        CrossOver(n);
        Evaluate(n, ch);
    }
    SortarePopulatieDupaFitness(0, n - 1);

    float highestFitnessValue = f(population[n - 1].gene, ch);

    return highestFitnessValue;

}




void Evaluate(int& n, char ch)
{

    float minim_curent = 1000000000;
    for (int i = 0; i < n; i++)
    {
        float val = f(population[i].gene, ch);
        for (int j = 0; j < LMAX; j++)

            if (val < minim_curent)
            {
                minim_curent = val;
            }
    }

    for (int i = 0; i < n; i++)
    {
        population[i].fitness = getFitness(population[i].gene, minim_curent, ch);
    }

}



void Mutate(int& n, char ch)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < L; i++) {
            float r = GenerareRandomFloat(0, 1);
            if (r < MUTATIE_PROB) {
                population[j].gene[i] = 1 - population[j].gene[i];
            }
        }
    }
}







void Select(int& n) //selectie turneu
{
    float T = 0;
    for (int i = 0; i < n; i++)
    {
        T += population[i].fitness;
    }
    float sel[DIM_POPULATIE * 10];

    sel[0] = 0;
    for (int i = 1; i <= n; i++)
    {
        sel[i] = population[i - 1].fitness / T;
    }

    for (int i = 1; i <= n; i++)
    {
        sel[i] += sel[i - 1];
    }

    membru vechia_populatie[250];
    int vechiul_n = n;

    for (int i = 0; i < vechiul_n; i++)
    {
        vechia_populatie[i] = population[i];
    }

    for (int i = 0; i < DIM_POPULATIE; i++)
    {
        float r = GenerareRandomFloat(0.000001, 1);
        int s;
        for (s = 0; s < n; s++)
        {
            if (sel[s] < r && r <= sel[s + 1])
            {
                break;
            }
        }
        population[i] = vechia_populatie[s];
    }
    n = DIM_POPULATIE;
}



float GenerareRandomFloat(float a, float b)
{
    uniform_real_distribution<> distribution3(a, b); //generam o probabilitate
    float x = distribution3(generator);
    return x;
}

int GenerareRandomInt(int a, int b)
{
    uniform_int_distribution<int> distribution1(a, b);
    int nr = distribution1(generator);
    return nr;
}


void GenerareSirDeBiti(int sir_biti[LMAX])
{
    for (int i = 0; i < L; i++)
    {
        int nr = GenerareRandomInt(1, 10);
        sir_biti[i] = nr % 2;

    }
}

void generatePopulation(int& n)
{
    for (int i = 0; i < DIM_POPULATIE; i++)
    {
        GenerareSirDeBiti(population[i].gene);
    }
    n = DIM_POPULATIE;
}


void SmallCrossOver(membru a, membru b)
{
    int pivot = GenerareRandomInt(1, L - 2);

    for (int i = 0; i <= pivot; i++)
    {
        children1.gene[i] = a.gene[i];
        children2.gene[i] = b.gene[i];
    }

    for (int i = pivot + 1; i < L; i++)
    {
        children1.gene[i] = b.gene[i];
        children2.gene[i] = a.gene[i];
    }
}

void CrossOver(int& n)
{
    for (int i = 0; i < n; i++)
    {
        population[i].prob_random = GenerareRandomFloat(0, 1); //probabilitate random
    }

    SortarePopulatieDupaProbabilitate(0, n - 1);

    int taietura = 0;

    while (population[taietura].prob_random < INCRUCISARE_PROB && taietura < n) {
        taietura++;
    }
    taietura--;

    if (taietura % 2 == 1)
    {
        int nr = GenerareRandomInt(1, 10);
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




float f(int sir_biti[LMAX], char ch)
{
    int ct = 0, j = 0;
    float X[CMAX];
    for (int i = 0; i < componente; i++)
        X[i] = 0;

    for (int i = L - 1; i >= 0; i--)
    {
        ct = ld;
        int p = 1;
        while (ct != 0)
        {
            X[j] = X[j] + sir_biti[i] * p;
            p = p * 2;
            ct--;
            i--;
        }
        j++;
        i++;
    }


    for (int i = 0; i < componente; i++)
    {
        float val = X[i] / (pow(2, ld) - 1);
        val *= (interval_end - interval_start);
        val += interval_start;
        X[i] = val;
    }

    if (ch == 'J') return DeJong(X);
    else if (ch == 'S') return Schwefel(X);
    else if (ch == 'R') return Rastrigin(X);
    else if (ch == 'M') return Michalewicz(X);
}

void TransformareBitiInCoordonate(float coordonate_minim[CMAX], int sir_biti[LMAX])
{
    int ct = 0, j = 0;
    for (int i = 0; i < componente; i++)
        coordonate_minim[i] = 0;

    for (int i = L - 1; i >= 0; i--)
    {
        ct = ld;
        int p = 1;
        while (ct != 0)
        {
            coordonate_minim[j] = coordonate_minim[j] + sir_biti[i] * p;
            p = p * 2;
            ct--;
            i--;
        }
        j++;
        i++;
    }

    for (int i = 0; i < componente; i++)
    {
        float val = coordonate_minim[i] / (pow(2, ld) - 1);
        val *= (interval_end - interval_start);
        val += interval_start;
        coordonate_minim[i] = val;
    }
}

long long CurrentTime() {
    milliseconds ms = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
        );
    return ms.count();
}


float DeJong(float X[LMAX])
{
    float val = 0;
    for (int i = 0; i < componente; i++)
        val += (X[i] * X[i]);
    return val;
}

float Schwefel(float X[CMAX])
{
    float val = 0;
    for (int i = 0; i < componente; i++)
        val += (( -1)* X[i] * sin(sqrt(abs(X[i]))));
    return val;
}

float Rastrigin(float X[CMAX])
{
    float val = 10 * componente;
    for (int i = 0; i < componente; i++)
        val += (X[i] * X[i] - 10 * cos(2.0 * 3.1415 * X[i]));
    return val;
}

float Michalewicz(float X[CMAX])
{
    float val = 0;
    for (int i = 0; i < componente; i++)
        val -= (sin(X[i]) * pow(sin((i * X[i] * X[i]) / 3.1415), 20));
    return val;
}

float getFitness(int sir_biti[LMAX], float minim_curent, char ch)
{
    float fitness = f(sir_biti, ch);
    if (minim_curent < 0)
        fitness += ((-1) * minim_curent);
    fitness += 100;
    fitness = 1 / fitness;
    fitness *= 2000;

    fitness = pow(10, fitness + 1);
    return fitness;
}

void DeJongsProcess()
{
    int t1 = CurrentTime();
    float s = 0;
    interval_start = -5.12;
    interval_end = 5.12;
    init(interval_start, interval_end);

    cout << "\nDe Jong\n";
    cout << "Valoarea minima este 0 in punctul x_i = 0  i = 1:n \n\n";

    float val = simulateEvolution('J');
    cout << val << "\n";
    s = val;
    float mmin = val;
    for (int i = 1; i < 30; i++)
    {
        val = simulateEvolution('J');
        cout << val << "\n";
        s += val;
        if (val < mmin) mmin = val;
    }
    cout << "\nValoarea minima: " << mmin << "\n\n";
    cout << "\nMedia: " << s / 30 << "\n\n";
    int t2 = CurrentTime();
    cout << "Time:" << " " << t2 - t1 << "\n\n";
}

void SchwefelProcess()
{
    int t1 = CurrentTime();
    float s = 0;
    interval_start = -500;
    interval_end = 500;
    init(interval_start, interval_end);

    cout << "\n\n\n\nSchwefel's\n";
    cout << "Valoarea minima este " << (-1) * componente * 418.9829 << " in punctul x_i = 420.9687  i = 1:n \n\n";

    float val = simulateEvolution('S');
    cout << val << "\n";
    s = val;
    float mmin = val;
    for (int i = 1; i < 30; i++)
    {
        val = simulateEvolution('S');
        cout << val << "\n";
        s += val;
        if (val < mmin) mmin = val;
    }
    cout << "\nValoarea minima: " << mmin << "\n\n";
    int t2 = CurrentTime();
    cout << "\nMedia: " << s / 30 << "\n\n";
    cout << "Time:" << " " << t2 - t1 << "\n\n";
}

void RastriginProcess()
{
    int t1 = CurrentTime();
    float s = 0;
    interval_start = -5.12;
    interval_end = 5.12;
    init(interval_start, interval_end);

    cout << "\n\n\n\nRastrigin\n";
    cout << "Valoarea minima este 0 in punctul x_i = 0  i = 1:n \n\n";

    float val = simulateEvolution('R');
    cout << val << "\n";
    s = val;
    float mmin = val;
    for (int i = 1; i < 30; i++)
    {
        val = simulateEvolution('R');
        cout << val << "\n";
        s += val;
        if (val < mmin) mmin = val;
    }
    cout << "\nValoarea minima: " << mmin << "\n\n";
    cout << "\nMedia: " << s / 30 << "\n\n";
    int t2 = CurrentTime();
    cout << "Time:" << " " << t2 - t1 << "\n\n";
}

void MichalewiczProcess()
{
    int t1 = CurrentTime();
    float s = 0;
    interval_start = 0;
    interval_end = 3.14;
    init(interval_start, interval_end);

    cout << "\n\n\n\nMichalewicz\n";
    cout << "valoarea minima este ?? in punctul x_i = ??  i = 1:n \n\n";

    float val = simulateEvolution('M');
    cout << val << "\n";
    s = val;
    float mmin = val;
    for (int i = 1; i < 30; i++)
    {
        val = simulateEvolution('M');
        cout << val << "\n";
        s += val;
        if (val < mmin) mmin = val;
    }
    cout << "\nValoarea minima: " << mmin << "\n\n";
    cout << "\nMedia: " << s / 30 << "\n\n";
    int t2 = CurrentTime();
    cout << "Time:" << " " << t2 - t1 << "\n\n";
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