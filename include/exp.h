// Glenn Ko giko@illinois.edu
// 2014.08.12
// Fixed-point Exponential Function

#ifndef __exp_h__
#define __exp_h__

#include <iostream>
#include <fixed_point.h>
#include <math.h>
#include <vector>

using namespace std;

typedef fpml::fixed_point<unsigned int, 16, 16> q16;
typedef fpml::fixed_point<unsigned int, 8, 24> q24;

void calculatelookuptable(vector<q16> &lut)
{
    // create a vector with all possible inputs

    FILE *file_data;
    char file_name[1000];
    sprintf(file_name, "lut.v");
    file_data = fopen(file_name,"w");

    q16 a;
    q16 b;
    float c;
    q16 d;
    unsigned int e;
    q16 f;
    unsigned int g;

    for (int i=0; i<256; i++)
    {
        a = q16(i);
        b = a >> q16(8);
        printf("a= %f\n",(float)a );
        printf("b= %f\n",(float)b );
        c = exp(float(b));
        printf("c= %f\n",(float)c );
        d = (q16)c;
        printf("d= %f\n",(float)d );

        //e = d <<= q16(16);
        e = unsigned((float)d * 65536);
        //printf("e= %f\n",(float)e );
        printf("e= %08x\n", e );
        fprintf(file_data, "%08x\n", e );

        f = d - (q16)1 - b;
        g = unsigned((float)f * 65536);
        printf("g= %08x\n", g );
        lut.push_back(f);
    }

    fclose(file_data);
}

int main()
{
    vector<q16> lut;
    calculatelookuptable(lut);

    q16 input   = 2.25125;
    q16 log2    = 0.6931471805599453994172321214;
    q16 divlog2 = 1.4426950408889634073599;

    unsigned int temp;
    temp = unsigned((float)log2 * 65536);
    printf("temp= %08x\n", temp);

    printf("input = %f\n", (float)input);
    printf("log2 = %f\n", (float)log2);
    printf("divlog2 = %f\n", (float)divlog2);

    q16 k;

    k = input * divlog2;
    printf("k = %f\n", (float)k);

    q16 rdK;
    rdK = (k >> q16(16)) << q16(16);
    printf("rdK = %f\n", (float)rdK);

    q16 X = rdK * log2;
    printf("X = %f\n", (float)X);

    q16 Y;
    Y = input - X;
    printf("Y = %f\n", (float)Y);
    // here it should actually be Q32 instead of Q16
    
    // this part should be replaced with LUT
    // e^Y = 1+Y+T(Y)
    q16 exY;
    exY = exp((float)Y);
    printf("exY = %f\n", (float)exY);

    q16 lut_val;
    q16 exYlut;
    unsigned int index;
    //index = Y <<= q16(16);
    index = (float)Y * 256;
    printf("index = %08x\n", (int)index);
    printf("index = %d\n", (int)index);
    lut_val = lut[index];
    printf("lut_val = %f\n", (float)lut_val);
    exYlut = (q16)1 + Y + lut_val;

    q16 twoK;
    q16 output;

    twoK = q16(1) << q16(3);
    printf("twoK = %f\n", (float)twoK);

    output = twoK * exY;
    printf("output = %f\n", (float)output);
    output = twoK * exYlut;
    printf("output lut = %f\n", (float)output);
    
    q16 real;
    real = exp((float)input);
    printf("read = %f\n", (float)real);

    return 0;
}

#endif // __exp_h__
