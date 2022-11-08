#include <iostream>
#include <random>
#include <chrono>
#include <math.h>

#define LMAX 2000
#define NMAX 100

using namespace std;
using namespace std::chrono;

default_random_engine generator;

int N = 30;
const int pre = 5; //precizia

float valoare_minima, minim_vector[NMAX];
float interval_start, interval_end;
float valoarea_globala;
int ld, L;

void init() {
    ld = ceil(log2(pow(10, pre) * (interval_end - interval_start)));
    L = N * ld;
}

float DeJong(float*);
float Schwefel(float*);
float Rastrigin(float*);
float Michalewicz(float*);

void DeJongsProcess();
void RastriginProcess();
void SchwefelProcess();
void MichalewiczProcess();

void SimulatedAnnealing(char);
void SimulatedAnnealing_iterat(char);
void HillClimbingFI(char);
void HillClimbingBI(char);
void HillClimbingWI(char);
void HillClimbing(bool, bool, char);

float f(int*, char);
void TransformareBitiInCoordonate(float*, int*);
void Verificare_vecini(float&, int*, bool, bool t, char, bool&);

long long CurrentTime();

int main()
{
    DeJongsProcess();
    SchwefelProcess();
    RastriginProcess();
    MichalewiczProcess();
}

void SimulatedAnnealing_iterat(char ch)
{
    cout << "\nSimulated Anneling\n";
    float rulari[31];
    float minim;
    float vector_minim[NMAX];
    int t1 = CurrentTime();
    SimulatedAnnealing(ch);
    minim = valoare_minima;
    rulari[0] = valoare_minima;
    float s = valoare_minima;
    for (int i = 1; i < 30; i++)
    {
        SimulatedAnnealing(ch);
        rulari[i] = valoare_minima;
        s = s + valoare_minima;
        if (valoare_minima < minim)
        {
            minim = valoare_minima;
            for (int j = 0; j < N; j++)
                vector_minim[j] = minim_vector[j];
        }
    }

    cout << "Media: " << s / 30 << "\n";

    cout << "\n\n";
    cout << "******************************\n";
    cout << "Valoarea minima este " << minim << "\n\nIn punctul:(";
    for (int i = 0; i < N - 1; i++)
        cout << vector_minim[i] << ", ";
    cout << vector_minim[N - 1];
    cout << ")\n\n";
    cout << "Media: " << s / 30 << "\n\n";
    float d = 0;
    for (int i = 0; i < 30; i++)
        d = d + (rulari[i] * (s / 30)) * (rulari[i] * (s / 30));
    d = d / 29;
    d = sqrt(d);
    cout << "Deviatia standard: " << d << "\n\n";
    int t2 = CurrentTime();
    cout << "Time(milis): " << t2 - t1 << "\n";
    cout << "******************************\n\n\n\n";
}

void SimulatedAnnealing(char ch)
{
    cout << "Simulated Annealing\n";

    int t1 = CurrentTime(); //timpul la care incepe sa ruleze algoritmul
    int nr_iteratii = 10000;
    int solutia_curenta[LMAX]; //vectorul de biti pentru o solutie curenta, de la un anumit pas

    //selectam random un candidat
    for (int k = 0; k < L; k++)
    {
        uniform_int_distribution<int> distribution(1, 10);
        int nr = distribution(generator);
        if (nr % 2 == 0) solutia_curenta[k] = 0;
        else solutia_curenta[k] = 1;
    }

    float valoarea_globala = f(solutia_curenta, ch); //calculam valoarea pt sirul candidat

    for (int i = 0; i < nr_iteratii; i++)
    {
        float T = 100; //temperatura

        //generam un sir de biti
        for (int k = 0; k < L; k++)
        {
            uniform_int_distribution<int> distribution1(1, 10);
            int nr = distribution1(generator);
            if (nr % 2 == 0) solutia_curenta[k] = 0;
            else solutia_curenta[k] = 1;
        }

        float valoarea_curenta = f(solutia_curenta, ch);
        int vecin[LMAX]; //vecin

        while (T > 10e-8)
        {
            for (int k = 0; k < L; k++)
                vecin[k] = solutia_curenta[k]; // vecinul e egal cu valoarea curenta

            uniform_int_distribution<int> distribution2(0, L - 1); //generam o pozitie aleatoare din vecin
            int pozitie_rand = distribution2(generator);
            vecin[pozitie_rand] = 1 - vecin[pozitie_rand]; //modificam un bit al vecinului
            float valoare_vecin = f(vecin, ch); //calculam valoarea vecinului

            float a = 0.0, b = 1.0 - 0.00001;
            uniform_real_distribution<> distribution3(a, b); //generam o probabilitate
            float x = distribution3(generator);

            if (valoare_vecin < valoarea_curenta)
            { //daca valoarea vecinului e mai aproape de minim
                for (int j = 0; j < L; j++)
                    solutia_curenta[j] = vecin[j]; //punem valuarea vecinului in valoarea curenta
                valoarea_curenta = valoare_vecin;
            }
            else if (x < exp(-abs(valoarea_curenta - valoare_vecin) / T))
            { //dam sansa unei solutii mai slabe
                for (int j = 0; j < L; j++)
                    solutia_curenta[j] = vecin[j];
                valoarea_curenta = valoare_vecin;
            }
            // cout << T << "\n";
            T *= 0.99; //micsoram temperatura
        }

        if (valoarea_curenta < valoarea_globala)
            valoarea_globala = valoarea_curenta; //daca valoarea gasita e mai buna ca minimul global curent, modificam
    }


    valoare_minima = valoarea_globala; //valoarea minima
    for (int i = 0; i < N; i++)
    { //punctul unde este valoarea minima
        TransformareBitiInCoordonate(minim_vector, solutia_curenta);
    }

    cout << "\n\n";
    cout << "******************************\n";
    cout << "Valoarea minima este " << valoare_minima << "\n\nIn punctul:(";
    for (int i = 0; i < N - 1; i++)
        cout << minim_vector[i] << ", ";
    cout << minim_vector[N - 1];
    cout << ")\n\n";
    int t2 = CurrentTime(); //timpul la care se termina
    cout << "Time(milis): " << t2 - t1 << "\n";
    cout << "******************************\n\n\n\n";
}




void HillClimbingFI(char ch) 
{
    cout << "\nHill Climbing - FirstImprovement\n";
    float rulari[31];
    float minim;
    float vector_minim[NMAX];
    int t1 = CurrentTime();
    HillClimbing(true, false, ch);
    minim = valoare_minima;
    rulari[0] = valoare_minima;
    float s = valoare_minima;
    for (int i = 1; i < 30; i++)
    {
        HillClimbing(true, false, ch);
        rulari[i] = valoare_minima;
        s = s + valoare_minima;
        if (valoare_minima < minim)
        {
            minim = valoare_minima;
            for (int j = 0; j < N; j++)
                vector_minim[j] = minim_vector[j];
        }
    }

    cout << "Media: " << s / 30 << "\n";

    cout << "\n\n";
    cout << "******************************\n";
    cout << "Valoarea minima este " << minim << "\n\nIn punctul:(";
    for (int i = 0; i < N - 1; i++)
        cout << vector_minim[i] << ", ";
    cout << vector_minim[N - 1];
    cout << ")\n\n";
    cout << "Media: " << s / 30 << "\n\n";
    float d = 0;
    for (int i = 0; i < 30; i++)
        d = d + (rulari[i] * (s / 30)) * (rulari[i] * (s / 30));
    d = d / 29;
    d = sqrt(d);
    cout << "Deviatia standard: " << d << "\n\n";
    int t2 = CurrentTime();
    cout << "Time(milis): " << t2 - t1 << "\n";
    cout << "******************************\n\n\n\n";
}


void HillClimbingBI(char ch)
{
    cout << "Hill Climbing - BestImprovement\n";

    float rulari[31];
    float minim;
    float vector_minim[NMAX];
    int t1 = CurrentTime();
    HillClimbing(false, true, ch);
    minim = valoare_minima;
    rulari[0] = valoare_minima;
    float s = valoare_minima;
    for (int i = 1; i < 30; i++)
    {
        HillClimbing(false, true, ch);
        rulari[i] = valoare_minima;
        s = s + valoare_minima;
        if (valoare_minima < minim)
        {
            minim = valoare_minima;
            for (int j = 0; j < N; j++)
                vector_minim[j] = minim_vector[j];
        }
    }

    cout << "Media: " << s / 30 << "\n";

    cout << "\n\n";
    cout << "******************************\n";
    cout << "Valoarea minima este " << minim << "\n\nIn punctul:(";
    for (int i = 0; i < N - 1; i++)
        cout << vector_minim[i] << ", ";
    cout << vector_minim[N - 1];
    cout << ")\n\n";
    cout << "Media: " << s / 30 << "\n\n";
    float d = 0;
    for (int i = 0; i < 30; i++)
        d = d + (rulari[i] * (s / 30)) * (rulari[i] * (s / 30));
    d = d / 29;
    d = sqrt(d);
    cout << "Deviatia standard: " << d << "\n\n";
    int t2 = CurrentTime();
    cout << "Time(milis): " << t2 - t1 << "\n";
    cout << "******************************\n\n\n\n";
}


void HillClimbingWI(char ch) {
    cout << "Hill Climbing - WorstImprovement\n";

    float rulari[31];
    float minim;
    float vector_minim[NMAX];
    int t1 = CurrentTime();
    HillClimbing(false, false, ch);
    minim = valoare_minima;
    cout << valoare_minima << "\n";
    rulari[0] = valoare_minima;
    float s = valoare_minima;
    for (int i = 1; i < 30; i++)
    {
        HillClimbing(false, false, ch);
        rulari[i] = valoare_minima;
        s = s + valoare_minima;
        if (valoare_minima < minim)
        {
            minim = valoare_minima;
            for (int j = 0; j < N; j++)
                vector_minim[j] = minim_vector[j];
        }
    }

    cout << "Media: " << s / 30 << "\n";

    cout << "\n\n";
    cout << "******************************\n";
    cout << "Valoarea minima este " << minim << "\n\nIn punctul:(";
    for (int i = 0; i < N - 1; i++)
        cout << vector_minim[i] << ", ";
    cout << vector_minim[N - 1];
    cout << ")\n\n";
    cout << "Media: " << s / 30 << "\n\n";
    float d = 0;
    for (int i = 0; i < 30; i++)
        d = d + (rulari[i] * (s / 30)) * (rulari[i] * (s / 30));
    d = d / 29;
    d = sqrt(d);
    cout << "Deviatia standard: " << d << "\n\n";
    int t2 = CurrentTime();
    cout << "Time(milis): " << t2 - t1 << "\n";
    cout << "******************************\n\n\n\n";
}


void HillClimbing(bool firstImprovement, bool bestImprovement, char ch)
{
    const int iteratii = 10000;
    int solutia_curenta[LMAX];

    for (int k = 0; k < L; k++)
    {
        uniform_int_distribution<int> distribution(1, 10);
        int nr = distribution(generator);
        if (nr % 2 == 0) solutia_curenta[k] = 0;
        else solutia_curenta[k] = 1;
    }

    float valoarea_globala = f(solutia_curenta, ch);
    int vector_biti[LMAX];
    int vecin[LMAX];

    for (int i = 0; i < iteratii; i++)
    {
        for (int k = 0; k < L; k++)
        {
            uniform_int_distribution<int> distribution1(1, 10);
            int nr = distribution1(generator);
            if (nr % 2 == 0) solutia_curenta[k] = 0;
            else solutia_curenta[k] = 1;
        }
        float valoarea_curenta = f(solutia_curenta, ch);

        bool done = false;
        while (!done)
            Verificare_vecini(valoarea_curenta, solutia_curenta, firstImprovement, bestImprovement, ch, done);

        if (valoarea_curenta < valoarea_globala)
        {
            valoarea_globala = valoarea_curenta;
            for (int j = 0; j < L; j++)
                vector_biti[j] = solutia_curenta[j];
        }

        valoare_minima = valoarea_globala; //valoarea minima
        for (int i = 0; i < N; i++)
        { //punctul unde este valoarea minima
            TransformareBitiInCoordonate(minim_vector, vector_biti);
        }
    }
}




float f(int sir_biti[LMAX], char ch)
{
    int ct = 0, j = 0;
    float X[NMAX];
    for (int i = 0; i < N; i++)
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


    for (int i = 0; i < N; i++)
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

void TransformareBitiInCoordonate(float coordonate_minim[NMAX], int sir_biti[LMAX])
{
    int ct = 0, j = 0;
    for (int i = 0; i < N; i++)
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

    for (int i = 0; i < N; i++)
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

void Verificare_vecini(float& valoare_curenta, int solutia_curenta[LMAX], bool firstImprovement, bool bestImprovement, char ch, bool& done)
{
    int vecin[LMAX];
    float valoare_de_referinta = valoare_curenta;
    bool found = 0;
    for (int i = 0; i < L; i++)
        vecin[i] = solutia_curenta[i];
    for (int i = 0; i < L; i++)
    {
            vecin[i] = 1-vecin[i];
            
            float valoare_vecin = f(vecin, ch);
            if (firstImprovement == 1)
            {
                if (valoare_vecin < valoare_curenta)
                {
                    valoare_curenta = valoare_vecin;
                    for (int k = 0; k < L; k++)
                        solutia_curenta[k] = vecin[k];
                    found = 1;
                    break;
                }
                ;
            }
            else if (bestImprovement == 1)
            {
                if (valoare_vecin < valoare_curenta)
                {
                    valoare_curenta = valoare_vecin;
                    for (int k = 0; k < L; k++)
                        solutia_curenta[k] = vecin[k];
                    found = 1;
                }
                
            }
            else if (bestImprovement == 0)
            {
                if (valoare_vecin > valoare_curenta && valoare_vecin < valoare_de_referinta)
                {
                    valoare_curenta = valoare_vecin;
                    for (int k = 0; k < L; k++)
                        solutia_curenta[k] = vecin[k];
                    found = 1;
                }
                
            }
            if (firstImprovement && found == 1)
                break;
            vecin[i] = 1-vecin[i];
           
    }
    done = 1 - found;
}




float DeJong(float X[LMAX])
{
    float val = 0;
    for (int i = 0; i < N; i++)
        val += (X[i] * X[i]);
    return val;
}

float Schwefel(float X[NMAX])
{
    float val = 0;
    for (int i = 0; i < N; i++)
        val += (-X[i] * sin(sqrt(abs(X[i]))));
    return val;
}

float Rastrigin(float X[NMAX])
{
    float val = 10 * N;
    for (int i = 0; i < N; i++)
        val += (X[i] * X[i] - 10 * cos(2.0 * 3.1415 * X[i]));
    return val;
}

float Michalewicz(float X[NMAX])
{
    float val = 0;
    for (int i = 0; i < N; i++)
        val -= (sin(X[i]) * pow(sin((i * X[i] * X[i]) / 3.1415), 20));
    return val;
}




void DeJongsProcess()
{
    interval_start = -5.12;
    interval_end = 5.12;
    init();

    cout << "\nDe Jong\n";
    cout << "Valoarea minima este 0 in punctul x_i = 0  i = 1:n \n\n";

    HillClimbingFI('J');
    HillClimbingBI('J');
    HillClimbingWI('J');
    SimulatedAnnealing('J');
}

void SchwefelProcess()
{
    interval_start = -500;
    interval_end = 500;
    init();

    cout << "\n\n\n\nSchwefel's\n";
    cout << "Valoarea minima este " << -N * 418.9829 << " in punctul x_i = 420.9687  i = 1:n \n\n";

    HillClimbingFI('S');
    HillClimbingBI('S');
    HillClimbingWI('S');
    SimulatedAnnealing('S');
}

void RastriginProcess()
{
    interval_start = -5.12;
    interval_end = 5.12;
    init();

    cout << "\n\n\n\nRastrigin\n";
    cout << "Valoarea minima este 0 in punctul x_i = 0  i = 1:n \n\n";

    HillClimbingFI('R');
    HillClimbingBI('R');
    HillClimbingWI('R');
    SimulatedAnnealing('R');
}

void MichalewiczProcess()
{
    interval_start = 0;
    interval_end = 3.14;
    init();

    cout << "\n\n\n\nMichalewicz\n";
    cout << "valoarea minima este ?? in punctul x_i = ??  i = 1:n \n\n";

    HillClimbingFI('M');
    HillClimbingBI('M');
    HillClimbingWI('M');
    SimulatedAnnealing('M');
}