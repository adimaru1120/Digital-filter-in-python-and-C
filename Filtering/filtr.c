#include <math.h>
#include <stdio.h>


#define N 1024
#define Nmaska 1023
float x[N];

float output[N];
int outputI[N];

float state[N][2] = {0};
int stateI[N][2] = {0};

float e1[N] = {0};
float e0[N] = {0};

int e1I[N] = {0};
int e0I[N] = {0};

int n;

float Bln[500];
float Bldn[500];
float Bisec[500][6];
int BisecI[500][6];

int BlnI[500];
int BldnI[500];


__declspec(dllexport) int lenB = 0;
__declspec(dllexport) int lenA = 0;
__declspec(dllexport) int ftype = 0;
__declspec(dllexport) int fstructure = 0;
__declspec(dllexport) int q = 0;
__declspec(dllexport) int sectionsNo = 0;


__declspec(dllexport) void init_filtr();
__declspec(dllexport) int sampleProcessor(int input_adc);
__declspec(dllexport) void importNumerator(float *Bn, size_t size);
__declspec(dllexport) void importDenominator(float *Bdn, size_t sizeA);
__declspec(dllexport) void importNumeratorI(int *Bn, size_t size);
__declspec(dllexport) void importDenominatorI(int *Bdn, size_t sizeA);
__declspec(dllexport) void importBiquad(float *Bdn, int section);
__declspec(dllexport) void importBiquadI(int *Bdn, int section);

int main(){
    return 0;
}

__declspec(dllexport) void init_filtr(void)
{
    n=0;
    for(int i = 0; i < 500; i++)
    {
        Bln[i] = 0;
        Bldn[i] = 0;
        BlnI[i] = 0;
        BldnI[i] = 0;
        for(int j = 0; j < 6; j++)
        {
            Bisec[i][j] = 0;
            BisecI[i][j] = 0;
        }
    }
}

__declspec(dllexport) int sampleProcessor(int input_adc)
{
    int i, probka_wyj, sumaI;
    x[n]=input_adc;

    float suma = 0;
    sumaI = 0;
    output[n] = 0;
    outputI[n] = 0;

    if(ftype == 0)          //FIR    
    {
        switch (fstructure)
        {
            case 0:         //direct
                for (i=0; i<lenB; i++)
                {
                    suma = suma + Bln[i]*x[(n-i)&Nmaska]; 
                }
                probka_wyj=roundf(suma);
                break;

            case 1:         //direct fixed point 
                for (i=0; i<lenB; i++)
                {
                    sumaI = sumaI + BlnI[i]*x[(n-i)&Nmaska]; 
                }
                probka_wyj = sumaI>>q;
                break;
        
            default:
                break;
        } 
    }
    else                    //IIR
    {
        switch (fstructure)
        {
            case 0:         //direct fp
                for(i = 0; i < lenB; i++)
                {
                    suma = suma + Bln[i]*x[(n-i)&Nmaska];
                }
                
                for(i = 1; i < lenA; i++)
                {
                    suma = suma - Bldn[i]*output[(n-i)&Nmaska];
                }
                output[n] = suma;
                probka_wyj=roundf(suma);
                break;

                case 1:     //lattice fp
                    e0[lenB - 1] = x[n];

                    for(i = lenB - 1; i > 0; i--)
                    {
                        e0[i-1] = e0[i] + (Bldn[i-1]*e1[i-1]);
                        e1[i] = ((-Bldn[i-1])*e0[i-1]) + e1[i-1]; 
                    }

                    e1[0] = e0[0];

                    for(i = 0; i < lenB ; i++)
                    {
                        suma += Bln[i]*e1[i];                
                    }
                    probka_wyj=roundf(suma);
                break;
        
                case 2:         //direct fixed point
                    for(i = 0; i < lenB; i++)
                    {
                        sumaI = sumaI + BlnI[i]*x[(n-i)&Nmaska];
                    }
                
                    for(i = 1; i < lenA; i++)
                    {
                        sumaI = sumaI - BldnI[i]*outputI[(n-i)&Nmaska];
                    }
                    outputI[n] = sumaI>>q;
                    probka_wyj = sumaI>>q;
                    break;

                case 3:         //lattice fixed point
                    e0I[lenB - 1] = x[n];

                    for(i = lenB - 1; i > 0; i--)
                    {
                        e0I[i-1] = e0I[i] + (BldnI[i-1]*e1I[i-1])>>q;
                        e1I[i] = ((-BldnI[i-1])*e0I[i-1])>>q + e1I[i-1];
                    }

                    e1I[0] = e0I[0];

                    for(i = 0; i < lenB ; i++)
                    {
                    sumaI = sumaI + BlnI[i]*e1I[i];                              
                    }
                    probka_wyj = sumaI >> q;                    //q15
                    break;

                case 4:         //biquad floating point 
                    float we, wy;
                    we = x[n];

                    for(int i = 0; i <sectionsNo; i++)
                    {
                        wy = Bisec[i][0]*we + state[i][0];
                        state[i][0] = Bisec[i][1]*we - Bisec[i][4]*wy + state[i][1];
                        state[i][1] = Bisec[i][2]*we  - Bisec[i][5]*wy;
                        we = wy;
                    }
                        
                    suma = wy;
                    probka_wyj=roundf(suma);
                    break;  

                case 5:
                    int weI, wyI;
                    weI = x[n];

                    for(int i = 0; i <sectionsNo; i++)
                    {
                        wyI = (BisecI[i][0]*weI)>>q + stateI[i][0];
                        stateI[i][0] = (BisecI[i][1]*weI)>>q - (BisecI[i][4]*wyI)>>q + stateI[i][1];
                        stateI[i][1] = (BisecI[i][2]*weI)>>q  - (BisecI[i][5]*wyI)>>q;
                        weI = wyI;
                    }
                    sumaI = wyI;
                    probka_wyj= sumaI;
                    break;  


            default:
                break;
        }

    }

    n++;
    n&=Nmaska;
    return probka_wyj;
}

__declspec(dllexport) void importNumerator(float* Bn, size_t size)
{
    lenB = size;
    for(int i = 0; i < size; i++)
    {
        Bln[i] = Bn[i];
    }
}

__declspec(dllexport) void importDenominator(float* Bdn, size_t sizeA)
{
    lenA = sizeA;
    for(int i = 0; i < sizeA; i++)
    {
        Bldn[i] = Bdn[i];
    }
}

__declspec(dllexport) void importNumeratorI(int* Bn, size_t size)
{
    lenB = size;
    for(int i = 0; i < size; i++)
    {
        BlnI[i] = Bn[i];
    }
}

__declspec(dllexport) void importDenominatorI(int* Bdn, size_t sizeA)
{
    lenA = sizeA;
    for(int i = 0; i < sizeA; i++)
    {
        BldnI[i] = Bdn[i];
    }
}

__declspec(dllexport) void importBiquad(float *Bdn, int section)
{

    for(int i = 0; i < 6; i++)
    {
        Bisec[section][i] = Bdn[i];
    }

}

__declspec(dllexport) void importBiquadI(int *Bdn, int section)
{

    for(int i = 0; i < 6; i++)
    {
            BisecI[section][i] = Bdn[i];
    }

}

