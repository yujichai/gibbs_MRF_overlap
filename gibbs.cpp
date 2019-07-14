/******************************************************************************/
//
//  sourcesep.cpp
//  Glenn Gihyun Ko
//  gko@seas.harvard.edu
//
//  source separation using mrf
//
/******************************************************************************/
// TODO:
// 1. randomize global labels cache

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <stdio.h>
#include "include/fixed_point.h"


using namespace fpml;
using namespace std;

typedef std::vector<std::vector<double> >    matrix_double;
typedef std::vector<std::vector<int> >      matrix_int;


template<typename OutputIt, typename Engine = std::mt19937>
void generate(OutputIt first, OutputIt last)
{
    static Engine eng;
    std::bernoulli_distribution dis;

    while (first != last)
    {
        *first++ = dis(eng);
    }
}

void Dump2DVector(std::vector<double> in_vec, int width, int height,
        char const kFileName[1000])
{
    ofstream vectorfile(kFileName);
    for(int h=0; h<height; h++)
    {
        for(int w=0; w<width; w++)
        {
            //vectorfile << in_vec[h*width+w] << "\t";
            vectorfile << in_vec[h*width+w] << " ";
        }
        vectorfile << endl;
    }
}

void Dump2DVectorTranspose(std::vector<double> in_vec, int width, int height,
        char const kFileName[1000])
{
    ofstream vectorfile(kFileName);
    for(int w=0; w<width; w++)
    {
        for(int h=0; h<height; h++)
        {
            //vectorfile << in_vec[h*width+w] << "\t";
            vectorfile << in_vec[h*width+w] << " ";
        }
        vectorfile << endl;
    }
}

void Dump3DVector(std::vector<double> in_vec, int width, int height, int num_lbl,
        char const kFileName[1000])
{
    ofstream vectorfile(kFileName);
    for(int h=0; h<height; h++)
    {
        for(int w=0; w<width; w++)
        {
            vectorfile << "[";
            for(int n=0; n<num_lbl; n++)
            {
                vectorfile << in_vec[(h*width+w)*num_lbl+n] << " ";
            }
            vectorfile << "]\t";
        }
        vectorfile << endl;
    }
}

void Copy2DVector(std::vector<double> in_vec, int f_w, int f_h, int s_w, int s_h,
        int e_w, int e_h, std::vector<double> &out_vec)
{
    int width;
    int height;
    width  = e_w;
    height = e_h;
    for(int h=0; h<height; h++)
    {
        for(int w=0; w<width; w++)
        {
            out_vec.at(h*width+w) = in_vec[((h+s_h)*f_w)+(w+s_w)];
            //cout << out_vec.at(h*width+w) << "\t" << in_vec[((h+s_h)*f_w)+(w+s_w)] << endl;
        }
    }
}

void Copy2DVector_(std::vector<double> in_vec, int f_w, int f_h, int s_w, int s_h,
        int e_w, int e_h, std::vector<double> &out_vec)
{
    int width;
    int height;
    width  = f_w;
    height = f_h;
    for(int h=0; h<height; h++)
    {
        for(int w=0; w<width; w++)
        {
            out_vec.at(((h+s_h)*e_w)+(w+s_w)) = in_vec[h*width+w];
            //cout << out_vec.at(h*width+w) << "\t" << in_vec[((h+s_h)*f_w)+(w+s_w)] << endl;
        }
    }
}

void Copy3DVector(std::vector<double> in_vec, int f_w, int f_h, int s_w, int s_h,
        int e_w, int e_h, int num_lbl, std::vector<double> &out_vec)
{
    int width;
    int height;
    width  = e_w - s_w;
    height = e_h - s_h;
    for(int h=0; h<height; h++)
    {
        for(int w=0; w<width; w++)
        {
            for(int n=0; n<num_lbl; n++)
            {
                out_vec.at(h*width*num_lbl+(w*num_lbl)+n) = in_vec[((h+s_h)*f_w*num_lbl)+(w+s_w)*num_lbl+n];
                //cout << n << "\t" << h*width*num_lbl+w+n << "\t" << in_vec[((h+s_h)*f_w)+(w+s_w)+n] << endl;
            }
        }
    }
}

void GibbsSampling(std::vector<double> datacost, std::vector<double> smoothcost,
        std::vector<double> &totalcost, std::vector<double> &label,
        int gibbs_iter, int width, int height, int num_lbl,
        double beta, double lambda, std::mt19937 gen)
{
    std::vector<double> p(num_lbl, 0);
    std::vector<double> p_sum(num_lbl, 0);
    double sum_neighbors;
    int north_label;
    int west_label;
    int east_label;
    int south_label;
    int new_label;
    double tcost;
    //double tcost_sum;
    double psum;
    double curr_uni_rand;
    double adj_uni_rand;

    north_label = 0;
    west_label = 0;
    east_label = 0;
    south_label = 0;
    new_label = 0;

    int label0 =0;
    int label1 =0;

    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for(int g=0; g<gibbs_iter; g++)
    {
        for(int h=0; h<height; h++)
        {
            for(int w=0; w<width; w++)
            {

                //cout << "Node: [" << h << "," << w << "]" << endl;
                // fetch nieghboring labels
                // need to check boundary conditions
                //tcost_sum = 0;
                psum = 0;
                curr_uni_rand = dis(gen);

                for(int n=0; n<num_lbl; n++)
                {
                    sum_neighbors = 0;
                    tcost = 0;

                    if (h>0)
                    {
                        north_label = label.at((h-1)*width+w);
                        sum_neighbors += smoothcost.at(n*num_lbl+north_label);
                    }
                    if (w>0)
                    {
                        west_label = label.at(h*width+(w-1));
                        sum_neighbors += smoothcost.at(n*num_lbl+west_label);
                    }
                    if (w<width-1)
                    {
                        east_label = label.at(h*width+(w+1));
                        sum_neighbors += smoothcost.at(n*num_lbl+east_label);
                    }
                    if (h<height-1)
                    {
                        south_label = label.at(((h+1)*width)+w);
                        sum_neighbors += smoothcost.at(n*num_lbl+south_label);
                    }

                    // there may be way to skip some computations
                    tcost =
                        datacost.at(h*width*num_lbl+(w*num_lbl)+n)
                        + (lambda * sum_neighbors);

                    totalcost.at((h*width*num_lbl)+(w*num_lbl)+n) = tcost;

                    p.at(n) = exp(-beta*tcost);
                    //cout << "D(" << n << ")= " << datacost.at(h*width*num_lbl+(w*num_lbl)+n) << endl;
                    //cout << "T(" << n << ")= " << tcost << endl;
                    //cout << "P(" << n << ")= " << p.at(n) << endl;
                    psum = psum  + p.at(n);
                    p_sum.at(n) = psum;
                }

                adj_uni_rand = curr_uni_rand * psum;
                //adj_uni_rand = 0.5 * psum;
                //adj_uni_rand = 0.000000001;

                //cout << "P Sum: [" << p_sum.at(0) << "," << p_sum.at(1) << "]" << endl;
                //cout << "Random = " << adj_uni_rand << endl;

                for(int n=0; n<num_lbl; n++)
                {
                    new_label = n;
                    if (p_sum.at(n) > adj_uni_rand)
                    {
                        break;
                    }
                }

                //cout << "New Label: " << new_label << endl;
                //if (new_label == 0)
                //    label0 ++;
                //else if (new_label ==1)
                //    label1 ++;

                label.at(h*width+w) = new_label;
            }
        }
        //cout << "Label0: " << label0 << endl;
        //cout << "Label1: " << label1 << endl;
    }
}

void GibbsSampling_c(std::vector<double> datacost, std::vector<double> smoothcost,
        std::vector<double> &totalcost, std::vector<double> &label,
        int gibbs_iter, int width, int height, int num_lbl,
        double beta, double lambda, std::mt19937 &gen)
{
    std::vector<double> p(num_lbl, 0);
    std::vector<double> p_sum(num_lbl, 0);
    double sum_neighbors;
    int north_label;
    int west_label;
    int east_label;
    int south_label;
    int new_label;
    double tcost;
    //double tcost_sum;
    double psum;
    double curr_uni_rand;
    double adj_uni_rand;

    fixed_point<unsigned int, 16, 16> beta_fp = beta;
    fixed_point<unsigned int, 16, 16> lambda_fp = lambda;
    std::vector<fixed_point<unsigned int, 8, 24>> p_fp(num_lbl, 0);
    std::vector<fixed_point<unsigned int, 8, 24>> p_sum_fp(num_lbl, 0);
    fixed_point<unsigned int, 16, 16> sum_neighbors_fp;
    fixed_point<unsigned int, 16, 16> tcost_fp;
    fixed_point<unsigned int, 8, 24> psum_fp;
    unsigned int curr_uni_rand_fp;
    fixed_point<unsigned int, 8, 24> adj_uni_rand_fp;
	fixed_point<unsigned int, 16, 16> dcost_fp;
	int new_label_fp = 0;

    double datacost_fp_temp = 0;
	double beta_fp_temp = 0;
	double tcost_fp_temp =0;
    double rand_fp_temp = 0;
    double psum_fp_temp = 0;

    // Testing variables
    int zero_cnt = 0;
    int diff_cnt = 0;

    north_label = 0;
    west_label = 0;
    east_label = 0;
    south_label = 0;
    new_label = 0;

    int label0 =0;
    int label1 =0;

    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for(int g=0; g<gibbs_iter; g++)
    {
        for(int color=0; color<2; color++)
        {
            for(int w=0; w<width; w++)
            {
                for(int h=(color==w%2)?0:1; h<height; h=h+2)
                {

                    //cout << "Node: [" << h << "," << w << "]" << endl;
                    // fetch nieghboring labels
                    // need to check boundary conditions
                    //tcost_sum = 0;
                    psum = 0;
                    curr_uni_rand = dis(gen);

                    psum_fp = 0;
                    curr_uni_rand_fp = curr_uni_rand * pow(2,32);

                    for(int n=0; n<num_lbl; n++)
                    {
                        sum_neighbors = 0;
                        tcost = 0;

                        sum_neighbors_fp = 0;
                        tcost_fp = 0;

                        if (h>0)
                        {
                            north_label = label.at((h-1)*width+w);
                            sum_neighbors += smoothcost.at(n*num_lbl+north_label);
                        }
                        if (w>0)
                        {
                            west_label = label.at(h*width+(w-1));
                            sum_neighbors += smoothcost.at(n*num_lbl+west_label);
                        }
                        if (w<width-1)
                        {
                            east_label = label.at(h*width+(w+1));
                            sum_neighbors += smoothcost.at(n*num_lbl+east_label);
                        }
                        if (h<height-1)
                        {
                            south_label = label.at(((h+1)*width)+w);
                            sum_neighbors += smoothcost.at(n*num_lbl+south_label);
                        }

                        // there may be way to skip some computations
                        tcost =
                            datacost.at(h*width*num_lbl+(w*num_lbl)+n)
                            + (lambda * sum_neighbors);

                        // totalcost.at((h*width*num_lbl)+(w*num_lbl)+n) = tcost;

                        p.at(n) = exp(-beta*tcost);
                        //cout << "D(" << n << ")= " << datacost.at(h*width*num_lbl+(w*num_lbl)+n) << endl;
                        //cout << "T(" << n << ")= " << tcost << endl;
                        //cout << "P(" << n << ")= " << p.at(n) << endl;
                        psum = psum  + p.at(n);
                        p_sum.at(n) = psum;

                        // Fixed point verion of finding p
                        datacost_fp_temp = datacost.at(h*width*num_lbl+(w*num_lbl)+n);

                        dcost_fp = datacost_fp_temp;
                        sum_neighbors_fp = sum_neighbors;

                        tcost_fp =
                            dcost_fp
                            + (lambda_fp * sum_neighbors_fp);

                        totalcost.at((h*width*num_lbl)+(w*num_lbl)+n) = tcost_fp;

                        beta_fp_temp = beta_fp;
                        tcost_fp_temp = tcost_fp;

                        p_fp.at(n) = exp(-beta_fp_temp*tcost_fp_temp);
                        psum_fp = psum_fp  + p_fp.at(n);
                        p_sum_fp.at(n) = psum_fp;

                        // Test the fix point p
                        //cout << "Label: " << n << endl;
                        //cout << "dcost_fp: " << dcost_fp << endl;
                        //cout << "tcost: " << tcost << endl;
                        //cout << "tcost_fp: " << tcost_fp << endl;
                        //cout << "psum: " << psum << endl;
                        //cout << "psum_fp: " << psum_fp << endl;

                        /*
                        if(psum_fp == (fixed_point<unsigned int, 8, 24>)0) {
                            cout << "Label: " << n << endl;
                            cout << "dcost_fp: " << dcost_fp << endl;
                            cout << "tcost: " << tcost << endl;
                            cout << "tcost_fp: " << tcost_fp << endl;
                            cout << "psum: " << psum << endl;
                            cout << "psum_fp: " << psum_fp << endl;
                        }
                        //*/

                    }

                    adj_uni_rand = curr_uni_rand * psum;
                    //adj_uni_rand = 0.000000001;
                    //adj_uni_rand = 0.5 * psum;
                    //adj_uni_rand = 0.0078125;

                    //cout << "P Sum: [" << p_sum.at(0) << "," << p_sum.at(1) << "]" << endl;
                    //cout << "Random = " << adj_uni_rand << endl;

                    for(int n=0; n<num_lbl; n++)
                    {
                        new_label = n;
                        if (p_sum.at(n) > adj_uni_rand)
                        {
                            break;
                        }
                    }

                    //cout << "New Label: " << new_label << endl;
                    //if (new_label == 0)
                    //    label0 ++;
                    //else if (new_label ==1)
                    //    label1 ++;

                    // label.at(h*width+w) = new_label;


                    // Fixed point verion of finding label
                    rand_fp_temp = curr_uni_rand_fp;
                    psum_fp_temp = psum_fp;
                    adj_uni_rand_fp = rand_fp_temp * psum_fp_temp / pow(2,32);

                    for(int n=0; n<num_lbl; n++)
                    {
                        new_label_fp = n;
                        if(p_sum_fp.at(n) > adj_uni_rand_fp)
                        {
                            break;
                        }
                    }

                    label.at(h*width+w) = new_label_fp;

                    // Test the fix point label
                    //cout << "curr_uni_rand: " << curr_uni_rand << endl;
                    //cout << "rand_fp_temp/2^32: " << rand_fp_temp / pow(2,32) << endl;
                    //cout << "psum: " << psum << endl;
                    //cout << "psum_fp_temp: " << psum_fp_temp << endl;
                    //cout << "adj_uni_rand: " << adj_uni_rand << endl;
                    //cout << "adj_uni_rand_fp: " << adj_uni_rand_fp << endl;
                    //cout << endl;

                    if(double(adj_uni_rand_fp) == 0) {
                    	zero_cnt ++;
                    }

                    if(new_label_fp != new_label) {
                    	diff_cnt ++;
                    }
                }
            }
        }
    }

    // Fixed point verion testing
    cout << "Zero Count: " << zero_cnt << endl;
    cout << "Difference Count: " << diff_cnt << endl;

}

void GibbsSampling_cn(std::vector<double> datacost, std::vector<double> smoothcost,
        std::vector<double> &totalcost, std::vector<double> &label,
        int gibbs_iter, int width, int height, int num_lbl,
        double beta, double lambda, std::mt19937 &gen)
{
    std::vector<double> p(num_lbl, 0);
    std::vector<double> p_sum(num_lbl, 0);
    double sum_neighbors;
    int north_label;
    int west_label;
    int east_label;
    int south_label;
    int new_label;
    double tcost;
    //double tcost_sum;
    double psum;
    double curr_uni_rand;
    double adj_uni_rand;

    fixed_point<unsigned int, 16, 16> beta_fp = beta;
    fixed_point<unsigned int, 16, 16> lambda_fp = lambda;
    std::vector<fixed_point<unsigned int, 8, 24>> p_fp(num_lbl, 0);
    std::vector<fixed_point<unsigned int, 8, 24>> p_sum_fp(num_lbl, 0);
    std::vector<fixed_point<unsigned int, 16, 16>> tcost_fp_v(num_lbl, 0);
    fixed_point<unsigned int, 16, 16> sum_neighbors_fp;
    fixed_point<unsigned int, 16, 16> tcost_fp;
    fixed_point<unsigned int, 8, 24> psum_fp;
    fixed_point<unsigned int, 16, 16> tcost_min_fp;
    unsigned int curr_uni_rand_fp;
    fixed_point<unsigned int, 8, 24> adj_uni_rand_fp;
    fixed_point<unsigned int, 16, 16> dcost_fp;
    int new_label_fp = 0;

    double beta_fp_temp = beta_fp;
    double tcost_fp_temp =0;
    double rand_fp_temp = 0;
    double psum_fp_temp = 0;

    // Testing variables
    int zero_cnt = 0;
    int diff_cnt = 0;

    north_label = 0;
    west_label = 0;
    east_label = 0;
    south_label = 0;
    new_label = 0;

    int label0 =0;
    int label1 =0;

    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for(int g=0; g<gibbs_iter; g++)
    {
        for(int color=0; color<2; color++)
        {
            for(int w=0; w<width; w++)
            {
                for(int h=(color==w%2)?0:1; h<height; h=h+2)
                {

                    //cout << "Node: [" << h << "," << w << "]" << endl;
                    // fetch nieghboring labels
                    // need to check boundary conditions
                    //tcost_sum = 0;
                    psum = 0;
                    curr_uni_rand = dis(gen);

                    psum_fp = 0;
                    curr_uni_rand_fp = curr_uni_rand * pow(2,32);

                    for(int n=0; n<num_lbl; n++)
                    {
                        sum_neighbors = 0;
                        tcost = 0;

                        sum_neighbors_fp = 0;
                        tcost_fp = 0;

                        if (h>0)
                        {
                            north_label = label.at((h-1)*width+w);
                            sum_neighbors += smoothcost.at(n*num_lbl+north_label);
                        }
                        if (w>0)
                        {
                            west_label = label.at(h*width+(w-1));
                            sum_neighbors += smoothcost.at(n*num_lbl+west_label);
                        }
                        if (w<width-1)
                        {
                            east_label = label.at(h*width+(w+1));
                            sum_neighbors += smoothcost.at(n*num_lbl+east_label);
                        }
                        if (h<height-1)
                        {
                            south_label = label.at(((h+1)*width)+w);
                            sum_neighbors += smoothcost.at(n*num_lbl+south_label);
                        }

                        // there may be way to skip some computations
                        tcost =
                            datacost.at(h*width*num_lbl+(w*num_lbl)+n)
                            + (lambda * sum_neighbors);

                        // totalcost.at((h*width*num_lbl)+(w*num_lbl)+n) = tcost;

                        p.at(n) = exp(-beta*tcost);
                        //cout << "D(" << n << ")= " << datacost.at(h*width*num_lbl+(w*num_lbl)+n) << endl;
                        //cout << "T(" << n << ")= " << tcost << endl;
                        //cout << "P(" << n << ")= " << p.at(n) << endl;
                        psum = psum  + p.at(n);
                        p_sum.at(n) = psum;

                        // Fixed point verion of finding p
                        dcost_fp = datacost.at(h*width*num_lbl+(w*num_lbl)+n);
                        sum_neighbors_fp = sum_neighbors;

                        tcost_fp =
                            dcost_fp
                            + (lambda_fp * sum_neighbors_fp);

                        totalcost.at((h*width*num_lbl)+(w*num_lbl)+n) = tcost_fp;

                        if(n == 0) {
                            tcost_min_fp = tcost_fp;
                        } 
                        else if(tcost_min_fp > tcost_fp) {
                            tcost_min_fp = tcost_fp;
                        }

                        tcost_fp_v.at(n) = tcost_fp;


                        // Test the fix point p
                        //cout << "Label: " << n << endl;
                        //cout << "tcost: " << tcost << endl;
                        //cout << "tcost_fp: " << tcost_fp << endl;
                        //cout << "psum: " << psum << endl;
                        //cout << "psum_fp: " << psum_fp << endl;

                    }

                    for(int n=0; n<num_lbl; n++) {
                        tcost_fp_temp = tcost_fp_v.at(n) - tcost_min_fp;

                        p_fp.at(n) = exp(-beta_fp_temp*tcost_fp_temp);
                        psum_fp = psum_fp  + p_fp.at(n);
                        p_sum_fp.at(n) = psum_fp;
                    }

                    adj_uni_rand = curr_uni_rand * psum;
                    //adj_uni_rand = 0.000000001;
                    //adj_uni_rand = 0.5 * psum;
                    //adj_uni_rand = 0.0078125;

                    //cout << "P Sum: [" << p_sum.at(0) << "," << p_sum.at(1) << "]" << endl;
                    //cout << "Random = " << adj_uni_rand << endl;

                    for(int n=0; n<num_lbl; n++)
                    {
                        new_label = n;
                        if (p_sum.at(n) > adj_uni_rand)
                        {
                            break;
                        }
                    }

                    //cout << "New Label: " << new_label << endl;
                    //if (new_label == 0)
                    //    label0 ++;
                    //else if (new_label ==1)
                    //    label1 ++;

                    // label.at(h*width+w) = new_label;


                    // Fixed point verion of finding label
                    rand_fp_temp = curr_uni_rand_fp;
                    psum_fp_temp = psum_fp;
                    adj_uni_rand_fp = rand_fp_temp * psum_fp_temp / pow(2,32);

                    for(int n=0; n<num_lbl; n++)
                    {
                        new_label_fp = n;
                        if(p_sum_fp.at(n) > adj_uni_rand_fp)
                        {
                            break;
                        }
                    }

                    label.at(h*width+w) = new_label_fp;

                    // Test the fix point label
                    //cout << "curr_uni_rand: " << curr_uni_rand << endl;
                    //cout << "rand_fp_temp/2^32: " << rand_fp_temp / pow(2,32) << endl;
                    //cout << "psum: " << psum << endl;
                    //cout << "psum_fp_temp: " << psum_fp_temp << endl;
                    //cout << "adj_uni_rand: " << adj_uni_rand << endl;
                    //cout << "adj_uni_rand_fp: " << adj_uni_rand_fp << endl;
                    //cout << endl;

                    if(double(adj_uni_rand_fp) == 0) {
                        zero_cnt ++;
                    }

                    if(new_label_fp != new_label) {
                        diff_cnt ++;
                    }
                }
            }
        }
    }

    // Fixed point verion testing
    cout << "Zero Count: " << zero_cnt << endl;
    cout << "Difference Count: " << diff_cnt << endl;

}

int main(int argc, char *argv[])
{

	// Fix point lib test
	//fixed_point<int, 16> a = 250;
    //fixed_point<int, 16> b = sqrt(a);

    //cout << "a:" << a << endl;
    //cout << "b:" << b << endl;

    int width  = 513;
    int height = 125;
    int num_lbl = 2;

    if (argc != 9)
    {
        cerr << "Not enough arguments!" << endl;
        return 0;
    }

    int em_iter = atoi(argv[1]);
    int gibbs_iter_tot = atoi(argv[2]);
    double beta = atof(argv[3]);
    double lambda = atof(argv[4]);
    int step = atof(argv[5]);
    int HW = atof(argv[6]);
    int run_num = atof(argv[7]);
    int loop = atof(argv[8]);


    cout << "========================================" << endl;
    cout << "EM iter: " << em_iter << endl;
    cout << "Gibbs iter: " << gibbs_iter_tot << endl;
    cout << "Beta: " << beta << endl;
    cout << "Lambda: " << lambda << endl;
    cout << "Step: " << step << endl;
    cout << "CPU/FPGA: " << HW << endl;
    cout << "Number of runs: " << run_num << endl;
    cout << "Loop: " << loop << endl;

    //double beta = 0.5;
    //double lambda = 0.5;

    // read in ild file
    char const *kILDfileName = "input/ild-noise.txt";
    ifstream ild_stream(kILDfileName);

    // store ild values in
    //matrix_double input_ild(height, vector<double>(width,0));
    std::vector<double> input_ild (height * width, 0);

    double buffer;
    for(int w=0; w<width; w++)
    {
        for(int h=0; h<height; h++)
        {
            ild_stream >> buffer;
            //input_ild[h][w] = buffer;
            // transposed input is stored into vector
            input_ild[h*width+w] = buffer;
        }
    }

    // for checking only
    //DumpMatrix(input_ild, height, width, "output/ild.txt");
    Dump2DVector(input_ild, width, height, "output/ild.txt");

    for(int gibbs_iter=(loop==1)?1:gibbs_iter_tot; gibbs_iter<=gibbs_iter_tot; gibbs_iter++) {
    for(int run=0; run<run_num; run++){

    //--------------------------------------------------------------------------
    // Initialize cost vectors and Gaussian parameters
    //--------------------------------------------------------------------------
//    double mean[] = {10, -30}; // binary
//    double vari[] = {30, 30}; // binary
    //double mean[] = {8, -8}; // binary
    //double vari[] = {8, 8}; // binary
    double mean_[] = {8, -8}; // binary
    double vari_[] = {8, 8}; // binary
    cout << "Initial mean : [" << mean_[0] << "," << mean_[1] << "]" << endl;

    std::vector<double> datacost(height * width * num_lbl, 0);
    std::vector<double> smoothcost(num_lbl * num_lbl, 0);

    //--------------------------------------------------------------------------
    // Finding data cost
    //--------------------------------------------------------------------------
//    for(int h=0; h<height; h++)
//    {
//        for(int w=0; w<width; w++)
//        {
//            for(int n=0; n<num_lbl; n++)
//            {
//                datacost.at((h*width+w)*num_lbl+n) =
//                    (input_ild.at(h*width+w) - mean[n]) *
//                    (input_ild.at(h*width+w) - mean[n]) /
//                    (2 * vari[n]);
//            }
//        }
//    }
//
//    Dump3DVector(datacost, width, height, num_lbl, "output/datacost.txt");

    //--------------------------------------------------------------------------
    // Finding smoothness cost function
    //--------------------------------------------------------------------------
    for(int n=0; n<num_lbl; n++)
    {
        for(int m=0; m<num_lbl; m++)
        {
            smoothcost.at(n*num_lbl+m) = (n-m)*(n-m);
        }
    }

    //Dump2DVector(smoothcost, num_lbl, num_lbl, "output/smoothcost.txt");


    //--------------------------------------------------------------------------
    // Create local copies
    //--------------------------------------------------------------------------

    // it is best if these values are randomized, global labels can be flexible
    // random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    // bernoulli
    std::bernoulli_distribution b_dis(0.5);

    //std::vector<double> global_labels(height * width);
    //std::generate(global_labels.begin(), global_labels.end(), int(d(gen));
    //std::generate(global_labels.begin(), global_labels.end(),
    //        [&gen, &b_dis] { return b_dis(gen);});
    //std::generate(global_labels.begin(), global_labels.end(),
    //         b_dis(gen));
    //std::vector<double> global_labels(height * width, 0);

    // random vs zero initialization
    //std::vector<double> global_labels(height * width);
    //generate(global_labels.begin(), global_labels.end());
    std::vector<double> global_labels(height * width, 0);

    //char const *LBLfileName = "input/init_labels.txt";
    //ifstream lbl_stream(LBLfileName);

    /*
    for(int h=0; h<height; h++)
    {
        for(int w=0; w<width; w++)
        {
            lbl_stream >> buffer;
            //input_ild[h][w] = buffer;
            // transposed input is stored into vector
            global_labels[h*width+w] = buffer;
        }
    }
    */

    //Dump2DVector(global_labels, 513, 125, "output/global_labels.txt");

    int local_width;
    int local_height;
    local_width  = 513;     // smallest power of 4 larger than 513
    local_height = step;
    //int step = 1;
    int s_w = 0;
    int s_h = 0;
    int diff = height - local_height;
    int shift = (diff / step) + 1;
    if(diff % step != 0) {
        shift++;
    }

    std::vector<double> local_labels(local_height * local_width,0);
    std::vector<double> local_dcosts(local_height * local_width * num_lbl,0);
    std::vector<double> totalcost(local_height * local_width * num_lbl, 0);

    double mean[num_lbl];
    double vari[num_lbl];

    for(int s = 0; s < shift; s++) {

        for(int k = 0; k < num_lbl; k++) {
            mean[k] = mean_[k];
            vari[k] = vari_[k];
        }

        if(s == shift - 1 ) {
            s_h = diff;
        }
        //s_h = 24;
        Copy2DVector(global_labels, width, height, s_w, s_h, local_width, local_height, local_labels);
        //Dump2DVector(local_labels, local_width, local_height, "output/local_labels.txt");

        //Copy3DVector(datacost, 513, 125, 0, 0,
        //        local_width, local_height, num_lbl, local_dcosts);
        //Dump3DVector(local_dcosts, local_width, local_height, num_lbl,
        //        "output/local_dcosts.txt");

        //--------------------------------------------------------------------------
        // EM
        //--------------------------------------------------------------------------

        double energy;
        double data_e;
        double smooth_e;
        int   node_label;
        int label_count[2] = {0, 0};
        double new_mean[2] = {0, 0};
        double new_vari[2] = {0, 0};

        double dcost_cap = 15;
        double dcost_diff = 0;
        double dcost_diff_temp = 0;

        for (int e=0; e<em_iter; e++)
        {
            // update datacostgit clone https://glennko@bitbucket.org/glennko/ssgibbs-matlab.git
            for(int h=0; h<local_height; h++)
            {
                for(int w=0; w<local_width; w++)
                {
                    dcost_diff = 0;

                    for(int n=0; n<num_lbl; n++)
                    {
                        local_dcosts.at((h*local_width+w)*num_lbl+n) =
                            (input_ild.at((h+s_h)*width+w) - mean[n]) *
                            (input_ild.at((h+s_h)*width+w) - mean[n]) /
                            (vari[n]);
                        dcost_diff_temp = local_dcosts.at((h*local_width+w)*num_lbl+n) - dcost_cap;
                        if(dcost_diff_temp > dcost_diff) {
                            dcost_diff = dcost_diff_temp;
                        }
                        //cout << local_dcosts.at((h*local_width+w)*num_lbl+n) << " ";
                    }

                    if(HW == 1 && dcost_diff > 0) {
                        for(int n=0; n<num_lbl; n++){
                            local_dcosts.at((h*local_width+w)*num_lbl+n) -= dcost_diff;
                            if(local_dcosts.at((h*local_width+w)*num_lbl+n) < 0) {
                                local_dcosts.at((h*local_width+w)*num_lbl+n) = 0;
                            }
                            //cout << local_dcosts.at((h*local_width+w)*num_lbl+n) << " ";
                        }
                    }

                    //cout << endl;
                }
            }
            //Dump3DVector(local_dcosts, local_width, local_height, num_lbl,"output/local_dcosts.txt");

            // Gibbs sampling
            // p for each label must be found
            // p = exp(-B*E)
            // E = DC*lambda*SC

        	if(HW == 0) {
				GibbsSampling(local_dcosts, smoothcost, totalcost, local_labels,
                    gibbs_iter, local_width, local_height, num_lbl, beta, lambda, gen);
			} 
			else if (HW == 1) {
				GibbsSampling_c(local_dcosts, smoothcost, totalcost, local_labels,
                    gibbs_iter, local_width, local_height, num_lbl, beta, lambda, gen);
			}
            else if (HW == 2) {
                GibbsSampling_cn(local_dcosts, smoothcost, totalcost, local_labels,
                    gibbs_iter, local_width, local_height, num_lbl, beta, lambda, gen);
            }

            // Update Gaussian parameters (source separation only)
            label_count[0] = 0;
            label_count[1] = 0;
            new_mean[0] = 0;
            new_mean[1] = 0;
            new_vari[0] = 0;
            new_vari[1] = 0;

            for(int h=0; h<local_height; h++)
            {
                for(int w=0; w<local_width; w++)
                {
                    if (local_labels.at(h*local_width+w) == 0)
                    {
                        new_mean[0] += input_ild[(h+s_h)*width+w];
                        new_vari[0] += input_ild[(h+s_h)*width+w]*input_ild[(h+s_h)*width+w];
                        label_count[0]++;
                    }
                    else if (local_labels.at(h*local_width+w) == 1)
                    {
                        new_mean[1] += input_ild[(h+s_h)*width+w];
                        new_vari[1] += input_ild[(h+s_h)*width+w]*input_ild[(h+s_h)*width+w];
                        label_count[1]++;
                    }
                }
            }

            for (int n=0; n<num_lbl; n++)
            {
                new_mean[n] = new_mean[n] / label_count[n];
                new_vari[n] = (new_vari[n] / label_count[n])
                    - (new_mean[n] * new_mean[n]);
            }
            //cout << "Label count = [" << label_count[0] << "," << label_count[1] << "]" << endl;
            //cout << "New mean = [" << new_mean[0] << "," << new_mean[1] << "]" << endl;
            //cout << "New variance = [" << new_vari[0] << "," << new_vari[1] << "]" << endl;

            for (int n=0; n<num_lbl; n++)
            {
                mean[n] = new_mean[n];
                vari[n] = new_vari[n];
            }

            //cout << fixed;
            //cout << "data_e\t\t\t" << "smooth_e\t\t\t" << "total_e\t\t\t" << endl;

            data_e = 0;
            for(int h=0; h<local_height; h++)
            {
                for(int w=0; w<local_width; w++)
                {
                    node_label = local_labels.at(h*local_width+w);
                    data_e += local_dcosts.at((h*local_width*num_lbl)+(w*num_lbl)+node_label);
                }
            }
            //cout << "data_e = " << data_e << endl;

            int other_label;
            smooth_e = 0;
            for(int h=0; h<local_height-1; h++)
            {
                for(int w=0; w<local_width-1; w++)
                {
                    node_label  = local_labels.at(h*local_width+w);
                    other_label = local_labels.at(h*local_width+(w+1));
                    smooth_e += smoothcost.at(node_label*num_lbl+other_label);
                    other_label = local_labels.at(((h+1)*local_width)+w);
                    smooth_e += smoothcost.at(node_label*num_lbl+other_label);
                }
            }
            node_label  = local_labels.at((local_height-1)*local_width+(local_width-1));
            other_label = local_labels.at((local_height-1)*local_width+(local_width-2));
            smooth_e += smoothcost.at(node_label*num_lbl+other_label);
            node_label  = local_labels.at((local_height-1)*local_width+(local_width-1));
            other_label = local_labels.at((local_height-2)*local_width+(local_width-1));
            smooth_e += smoothcost.at(node_label*num_lbl+other_label);
            //cout << "smooth_e = " << smooth_e << endl;
            smooth_e = smooth_e * lambda;
            //cout << "smooth_e = " << smooth_e << endl;

            //cout << "smooth_e = " << smooth_e << endl;
            //cout << "total_e = " << data_e + smooth_e << endl;

        }
        cout << endl;
        cout << s_h << endl;
        cout << "Shift: " << (s + 1) << "/" << shift << endl;

        // Report Gaussian parameters
        cout << "Label count = [" << label_count[0] << "," << label_count[1] << "]" << endl;
        cout << "New mean = [" << new_mean[0] << "," << new_mean[1] << "]" << endl;
        cout << "New variance = [" << new_vari[0] << "," << new_vari[1] << "]" << endl;

        // Report energies
        cout << fixed;
        cout << "data_e\t\t" << "smooth_e\t" << "total_e\t\t" << endl;
        cout << data_e << "\t";
        cout << smooth_e << "\t";
        cout << data_e + smooth_e << "\t" << endl;

        // Updating the global_labels using new local_labels
        Copy2DVector_(local_labels, local_width, local_height, s_w, s_h, width, height, global_labels);

        // Updating the shift amount in height
        s_h = s_h + step;

    }

    //Dump3DVector(totalcost, local_width, local_height, num_lbl, "output/totalcosts.txt");
    string s;
    char buffer [100];
	if(HW == 1) {
		// s = "output/graph_overlap_" + to_string((int)em_iter) + "_" + to_string((int)gibbs_iter) + "_" + to_string((int)step) + "_" + to_string((int)(run+1)) + ".txt";
		sprintf (buffer, "output/labels_FPGA_FP/graph_overlap_%g_%g_%d_%d_%d_%d.txt", beta, lambda, em_iter, gibbs_iter, step, (run+1));
		s = buffer;
		//s = "output/SDR_FPGA_FP/graph_overlap_" + to_string((float)beta) + "_" + to_string((float)lambda) + "_" + to_string((int)em_iter) + "_" + to_string((int)gibbs_iter) + "_" + to_string((int)step) + "_" + to_string((int)(run+1)) + ".txt";
	} 
	else if (HW == 0){
		// s = "output/graph_overlap_" + to_string((int)em_iter) + "_" + to_string((int)gibbs_iter) + "_" + to_string((int)step) + "_" + to_string((int)(run+1)) + ".txt";
        sprintf (buffer, "output/labels_CPU/graph_overlap_%g_%g_%d_%d_%d_%d.txt", beta, lambda, em_iter, gibbs_iter, step, (run+1));
        s = buffer;
	}
    else if(HW == 2) {
        // s = "output/graph_overlap_" + to_string((int)em_iter) + "_" + to_string((int)gibbs_iter) + "_" + to_string((int)step) + "_" + to_string((int)(run+1)) + ".txt";
        sprintf (buffer, "output/labels_FPGA_FP_CN/graph_overlap_%g_%g_%d_%d_%d_%d.txt", beta, lambda, em_iter, gibbs_iter, step, (run+1));
        s = buffer;
        //s = "output/SDR_FPGA_FP/graph_overlap_" + to_string((float)beta) + "_" + to_string((float)lambda) + "_" + to_string((int)em_iter) + "_" + to_string((int)gibbs_iter) + "_" + to_string((int)step) + "_" + to_string((int)(run+1)) + ".txt";
    } 
   
    //s = "output/new_labels.txt";

    Dump2DVector(global_labels, width, height, s.c_str());
    //Dump2DVectorTranspose(global_labels, width, height, "output/new_labels_t.txt");
    }

    // Gibbs Iter Loop
    }
}
