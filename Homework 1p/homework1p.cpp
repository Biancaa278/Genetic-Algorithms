#include <iostream>
#include <random>
#include <chrono>
#include <math.h>

#define LMAX 2000
#define NMAX 100

using namespace std;
using namespace std::chrono;

default_random_engine generator;

int N = 1;
const int pre = 0; //precizia

float valoare_maxima, maxim_vector[NMAX];
float interval_start, interval_end;
float valoarea_globala;
int ld, L;

void init() {
    ld = 5;
    L = N * ld;
}

float funct(float*);


void functProcess();


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
    functProcess();
    
}




void HillClimbingFI(char ch)
{
    cout << "\nHill Climbing - FirstImprovement\n";
    float rulari[31];
    float maxim;
    float vector_maxim[NMAX];
    int t1 = CurrentTime();
    HillClimbing(true, false, ch);
    maxim = valoare_maxima;
    rulari[0] = valoare_maxima;
    float s = valoare_maxima;
    for (int i = 1; i < 30; i++)
    {
        HillClimbing(true, false, ch);
        rulari[i] = valoare_maxima;
        s = s + valoare_maxima;
        if (valoare_maxima > maxim)
        {
            maxim = valoare_maxima;
            for (int j = 0; j < N; j++)
                vector_maxim[j] = maxim_vector[j];
        }
    }


    cout << "Valoarea maxima este " << maxim << "\n\n";
    
}


void HillClimbingBI(char ch)
{
    cout << "Hill Climbing - BestImprovement\n";

    float rulari[31];
    float maxim;
    float vector_maxim[NMAX];
    int t1 = CurrentTime();
    HillClimbing(false, true, ch);
    maxim = valoare_maxima;
    rulari[0] = valoare_maxima;
    float s = valoare_maxima;
    for (int i = 1; i < 30; i++)
    {
        HillClimbing(false, true, ch);
        rulari[i] = valoare_maxima;
        s = s + valoare_maxima;
        if (valoare_maxima < maxim)
        {
            maxim = valoare_maxima;
            for (int j = 0; j > N; j++)
                vector_maxim[j] = maxim_vector[j];
        }
    }
    cout << "Valoarea maxima este " << maxim << "\n\n";
}




void HillClimbing(bool firstImprovement, bool bestImprovement, char ch)
{
    const int iteratii = 10000;
    int solutia_curenta[LMAX];
    float coord[2];
    solutia_curenta[0] = 0;
    solutia_curenta[1] = 1;
    solutia_curenta[2] = 1;
    solutia_curenta[3] = 0;
    solutia_curenta[4] = 0;

    float valoarea_globala = f(solutia_curenta, ch);
    int vector_biti[LMAX];
    int vecin[LMAX];

    for (int i = 0; i < iteratii; i++)
    {
        solutia_curenta[0] = 0;
        solutia_curenta[1] = 1;
        solutia_curenta[2] = 1;
        solutia_curenta[3] = 0;
        solutia_curenta[4] = 0;
        float valoarea_curenta = f(solutia_curenta, ch);

        bool done = false;
        while (!done)
            Verificare_vecini(valoarea_curenta, solutia_curenta, firstImprovement, bestImprovement, ch, done);

        if (valoarea_curenta > valoarea_globala)
        {
            valoarea_globala = valoarea_curenta;
            for (int j = 0; j < L; j++)
                vector_biti[j] = solutia_curenta[j];
        }
       
        valoarea_curenta = f(solutia_curenta, ch);
        valoare_maxima = valoarea_globala; //valoarea maxima
        for (int i = 0; i < N; i++)
        { //punctul unde este valoarea maxima
            TransformareBitiInCoordonate(maxim_vector, vector_biti);
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

    return funct(X);
    
}

void TransformareBitiInCoordonate(float coordonate_maxim[NMAX], int sir_biti[LMAX])
{
    int ct = 0, j = 0;
    for (int i = 0; i < N; i++)
        coordonate_maxim[i] = 0;

    for (int i = L - 1; i >= 0; i--)
    {
        ct = ld;
        int p = 1;
        while (ct != 0)
        {
            coordonate_maxim[j] = coordonate_maxim[j] + sir_biti[i] * p;
            p = p * 2;
            ct--;
            i--;
        }
        j++;
        i++;
    }

    for (int i = 0; i < N; i++)
    {
        float val = coordonate_maxim[i] / (pow(2, ld) - 1);
        val *= (interval_end - interval_start);
        val += interval_start;
        coordonate_maxim[i] = val;
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
        vecin[i] = 1 - vecin[i];

        float valoare_vecin = f(vecin, ch);
        if (firstImprovement == 1)
        {
            if (valoare_vecin > valoare_curenta)
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
            if (valoare_vecin > valoare_curenta)
            {
                valoare_curenta = valoare_vecin;
                for (int k = 0; k < L; k++)
                    solutia_curenta[k] = vecin[k];
                found = 1;
            }

        }
        if (firstImprovement && found == 1)
            break;
        vecin[i] = 1 - vecin[i];

    }
    done = 1 - found;
}




float funct(float X[LMAX])
{
    float val = 0;
    for (int i = 0; i < N; i++)
    {
        val += (X[i] * X[i] * X[i]);
        val -= (60 * X[i] * X[i]);
        val += (900 * X[i]);
        val += 100;
    }
    return val;
}




void functProcess()
{
    interval_start =0;
    interval_end = 31;
    init();

    HillClimbingFI('J');
    HillClimbingBI('J');
}

