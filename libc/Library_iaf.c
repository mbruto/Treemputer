#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

struct mat
{
 
    double value;
    int index;
    
};

struct tuple
{
 
    double value_1;
    double value_2;
    
};

int cmp(const void *a, const void *b)
{
    struct mat *a1 = (struct mat *)a;
    struct mat *a2 = (struct mat *)b;
    if ((*a1).value < (*a2).value)
        return -1;
    else if ((*a1).value > (*a2).value)
        return 1;
    else
        return 0;
    
}

double min_dist_imputation(int i, int j, int number, double sub_add_dist[number][number])
{
    
    double MinDist= 1000.0;
    double test_dist = 0.0;
    double opti_dist = -99.0;
    double dist_1 = 0.0;
    double dist_2 = 0.0;
    int min_k = 0;
    int min_l = 0;
            
    for (int k=0; k<number; k++)
    {
        if(sub_add_dist[i][k] == -99.0 || sub_add_dist[j][k] == -99.0) continue;
                
        for (int l=k+1; l<number; l++)
        {
            if(k==l || k==i || k==j || l==i || l==j || sub_add_dist[k][l]==-99.0 || sub_add_dist[i][l]==-99.0 || sub_add_dist[j][l]==-99.0) continue;
                    
            test_dist = sub_add_dist[i][k] + sub_add_dist[j][k] + sub_add_dist[i][l] + sub_add_dist[j][l] + sub_add_dist[k][l];
                    
            if (test_dist < MinDist)
            {
                        
                dist_1 = sub_add_dist[i][k] + sub_add_dist[j][l];
                dist_2 = sub_add_dist[j][k] + sub_add_dist[i][l];
                        
                if (fabs(dist_1 - dist_2) < 0.00001)
                {
                }
                else if (dist_1 > dist_2)
                {
                    MinDist = test_dist;
                    opti_dist = sub_add_dist[i][k] + sub_add_dist[j][l] - sub_add_dist[k][l];
                    min_k = k;
                    min_l = l;
                    
                }
                else if (dist_1 < dist_2)
                {
                    MinDist = test_dist;
                    opti_dist = sub_add_dist[j][k] + sub_add_dist[i][l] - sub_add_dist[k][l];
                    min_k = k;
                    min_l = l;
                    
                }
            }
        }
    }
    return opti_dist;
}

/*
double distance_estimation_quartet(int m, int n, int number, double sub_matrix[number][number])
{
    
    FILE* log_file = NULL ;
    log_file = fopen("log_imput.txt", "a");
    
    double temp_dist = 0.0;
    double opti_dist = -99.0;
    int num_kl = 0;
    
    struct mat ik_min_list[number];
    struct mat jl_min_list[number];
    
    for (int x=0; x<number; x++)
    {
        if(sub_matrix[m][x] != -99.0 && sub_matrix[n][x] != -99.0)
        {
                
            ik_min_list[num_kl].value = sub_matrix[m][x];
            ik_min_list[num_kl].index = x;
            
            jl_min_list[num_kl].value = sub_matrix[n][x];
            jl_min_list[num_kl].index = x;
            
            num_kl += 1;
            
        }
    }

    qsort(ik_min_list, num_kl, sizeof(ik_min_list[0]), cmp);
    qsort(jl_min_list, num_kl, sizeof(jl_min_list[0]), cmp);
    
    int k_index = 0;
    int l_index = 0;
    
    
    for (int k=1; k<num_kl; k++)
    {
        k_index = ik_min_list[k].index;
        
        for (int l=1; l<num_kl; l++)
        {
            
            l_index = jl_min_list[l].index;
            
            if(k_index == l_index) continue;
            if(sub_matrix[k_index][l_index] < 0 || sub_matrix[m][k_index] < 0 || sub_matrix[m][l_index] < 0 || sub_matrix[n][k_index] < 0 || sub_matrix[n][l_index] < 0) continue;
            if(sub_matrix[k_index][l_index] < sub_matrix[m][k_index] && sub_matrix[k_index][l_index] < sub_matrix[m][l_index] && sub_matrix[k_index][l_index] < sub_matrix[n][k_index] && sub_matrix[k_index][l_index] < sub_matrix[n][l_index]) continue;
            
            double d_1 = 0.0;
            double d_2 = 0.0;
            double d_3 = 0.0;
            double d_4 = 0.0;
            double d_5 = 0.0;
            
            fprintf(log_file, "-----------\n");
            fprintf(log_file, "m = %d - n = %d\n", m, n);
            fprintf(log_file, "%d, %d\n", k, l);
            fprintf(log_file, "%d %d\n", k_index, l_index);
            fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", sub_matrix[m][k_index], sub_matrix[m][l_index], sub_matrix[n][k_index], sub_matrix[n][l_index], sub_matrix[k_index][l_index]);
            fprintf(log_file, "-----------\n");
            
            if(sub_matrix[n][k_index] < sub_matrix[n][l_index] && sub_matrix[m][l_index] < sub_matrix[m][k_index] && sub_matrix[n][k_index] < sub_matrix[k_index][l_index] && sub_matrix[m][l_index] < sub_matrix[k_index][l_index])
            {
                
                if(sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index] > 0)
                {
                    //DIVISER PAR DEUX BORDEL !!!!
                    d_1 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index]);
                    d_2 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[k_index][l_index] - sub_matrix[m][k_index]);
                    d_3 = 0.5 * (sub_matrix[n][l_index] + sub_matrix[n][k_index] - sub_matrix[k_index][l_index]);
                    d_4 = 0.5 * (sub_matrix[n][k_index] + sub_matrix[k_index][l_index] - sub_matrix[n][l_index]);
                    d_5 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[n][l_index] - sub_matrix[m][l_index] - sub_matrix[n][k_index]);
                    
                    if (d_5 < d_1 && d_5 < d_2 && d_5 < d_3 && d_5 < d_4 ) continue;
                    
                    if (d_1 > 0 && d_2 > 0 && d_3 > 0 && d_4 > 0 && d_5 > 0)
                    {
                        
                        opti_dist = sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index];
                        
                        fprintf(log_file, "-----------\n");
                        fprintf(log_file, "!!!UN %d\n", n);
                        fprintf(log_file, "m = %d - n = %d\n", m, n);
                        fprintf(log_file, "%d, %d\n", k, l);
                        fprintf(log_file, "%d %d\n", k_index, l_index);
                        fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", d_1, d_2, d_3, d_4, d_5);
                        fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", sub_matrix[m][k_index], sub_matrix[m][l_index], sub_matrix[n][k_index], sub_matrix[n][l_index], sub_matrix[k_index][l_index]);
                        fprintf(log_file, "%lf\n", opti_dist);
                        fprintf(log_file, "-----------\n");
                        
                    }
                }
            }
            
            else if(sub_matrix[n][k_index] > sub_matrix[n][l_index] && sub_matrix[m][l_index] > sub_matrix[m][k_index] && sub_matrix[n][l_index] < sub_matrix[k_index][l_index] && sub_matrix[m][k_index] < sub_matrix[k_index][l_index])
            {
                
                if(sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index] > 0)
                {
                    //DIVISER PAR DEUX BORDEL !!!!
                    d_1 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index]);
                    d_2 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[k_index][l_index] - sub_matrix[m][l_index]);
                    d_3 = 0.5 * (sub_matrix[n][k_index] + sub_matrix[n][l_index] - sub_matrix[k_index][l_index]);
                    d_4 = 0.5 * (sub_matrix[n][l_index] + sub_matrix[k_index][l_index] - sub_matrix[n][k_index]);
                    d_5 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[n][k_index] - sub_matrix[m][k_index] - sub_matrix[n][l_index]);
                    
                    if (d_5 < d_1 && d_5 < d_2 && d_5 < d_3 && d_5 < d_4 ) continue;
                    
                    if(d_1 > 0 && d_2 > 0 && d_3 > 0 && d_4 > 0 && d_5 > 0)
                    {
                            
                        opti_dist = sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index];
                        
                        fprintf(log_file, "-----------\n");
                        fprintf(log_file, "!!!DEUX %d\n", n);
                        fprintf(log_file, "m = %d - n = %d\n", m, n);
                        fprintf(log_file, "%d, %d\n", k, l);
                        fprintf(log_file, "%d %d\n", k_index, l_index);
                        fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", d_1, d_2, d_3, d_4, d_5);
                        fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", sub_matrix[m][k_index], sub_matrix[m][l_index], sub_matrix[n][k_index], sub_matrix[n][l_index], sub_matrix[k_index][l_index]);
                        fprintf(log_file, "%lf\n", opti_dist);
                        fprintf(log_file, "-----------\n");
                        
                    }
                }
            }
            
            if(opti_dist != -99.0) break;
            
        }
        
        if(opti_dist != -99.0) break;
        
    }
    fclose(log_file);
    return opti_dist;
    
}
*/

double distance_estimation_quartet(int m, int n, int number, double sub_matrix[number][number])
{
    
    FILE* log_file = NULL ;
    log_file = fopen("log_imput.txt", "a");
    
    double temp_dist = 0.0;
    double opti_dist = -99.0;
    int num_kl = 0;
    
    struct mat ik_min_list[number];
    struct mat jl_min_list[number];
    
    for (int x=0; x<number; x++)
    {
        if(sub_matrix[m][x] != -99.0 && sub_matrix[n][x] != -99.0)
        {
                
            ik_min_list[num_kl].value = sub_matrix[m][x];
            ik_min_list[num_kl].index = x;
            
            jl_min_list[num_kl].value = sub_matrix[n][x];
            jl_min_list[num_kl].index = x;
            
            num_kl += 1;
            
        }
    }

    qsort(ik_min_list, num_kl, sizeof(ik_min_list[0]), cmp);
    qsort(jl_min_list, num_kl, sizeof(jl_min_list[0]), cmp);
    
    int k_index = 0;
    int l_index = 0;
    
    
    for (int k=1; k<num_kl; k++)
    {
        k_index = ik_min_list[k].index;
        
        for (int l=1; l<num_kl; l++)
        {
            
            l_index = jl_min_list[l].index;
            
            if(k_index == l_index) continue;
            if(sub_matrix[k_index][l_index] < 0 || sub_matrix[m][k_index] < 0 || sub_matrix[m][l_index] < 0 || sub_matrix[n][k_index] < 0 || sub_matrix[n][l_index] < 0) continue;
            if(sub_matrix[k_index][l_index] < sub_matrix[m][k_index] && sub_matrix[k_index][l_index] < sub_matrix[m][l_index] && sub_matrix[k_index][l_index] < sub_matrix[n][k_index] && sub_matrix[k_index][l_index] < sub_matrix[n][l_index]) continue;
            
            double d_1 = 0.0;
            double d_2 = 0.0;
            double d_3 = 0.0;
            double d_4 = 0.0;
            double d_5 = 0.0;
            
            if(sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index] > 0 && sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index] < 0)
            {
                    //DIVISER PAR DEUX BORDEL !!!!
                    d_1 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index]);
                    d_2 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[k_index][l_index] - sub_matrix[m][k_index]);
                    d_3 = 0.5 * (sub_matrix[n][l_index] + sub_matrix[n][k_index] - sub_matrix[k_index][l_index]);
                    d_4 = 0.5 * (sub_matrix[n][k_index] + sub_matrix[k_index][l_index] - sub_matrix[n][l_index]);
                    d_5 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[n][l_index] - sub_matrix[m][l_index] - sub_matrix[n][k_index]);
                    
                    if (d_5 < d_1 && d_5 < d_2 && d_5 < d_3 && d_5 < d_4 ) continue;
                    
                    if (d_1 > 0 && d_2 > 0 && d_3 > 0 && d_4 > 0 && d_5 > 0)
                    {
                        
                        opti_dist = sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index];
                        
                        fprintf(log_file, "-----------\n");
                        fprintf(log_file, "!!!UN %d\n", n);
                        fprintf(log_file, "m = %d - n = %d\n", m, n);
                        fprintf(log_file, "%d, %d\n", k, l);
                        fprintf(log_file, "%d %d\n", k_index, l_index);
                        fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", d_1, d_2, d_3, d_4, d_5);
                        fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", sub_matrix[m][k_index], sub_matrix[m][l_index], sub_matrix[n][k_index], sub_matrix[n][l_index], sub_matrix[k_index][l_index]);
                        fprintf(log_file, "%lf\n", opti_dist);
                        fprintf(log_file, "-----------\n");
                        
                    }
            }
            
            if(sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index] > 0 && sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index] > 0)
            {
                    //DIVISER PAR DEUX BORDEL !!!!
                    d_1 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index]);
                    d_2 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[k_index][l_index] - sub_matrix[m][l_index]);
                    d_3 = 0.5 * (sub_matrix[n][k_index] + sub_matrix[n][l_index] - sub_matrix[k_index][l_index]);
                    d_4 = 0.5 * (sub_matrix[n][l_index] + sub_matrix[k_index][l_index] - sub_matrix[n][k_index]);
                    d_5 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[n][k_index] - sub_matrix[m][k_index] - sub_matrix[n][l_index]);
                    
                    if (d_5 < d_1 && d_5 < d_2 && d_5 < d_3 && d_5 < d_4 ) continue;
                    
                    if(d_1 > 0 && d_2 > 0 && d_3 > 0 && d_4 > 0 && d_5 > 0)
                    {
                            
                        opti_dist = sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index];
                        
                        fprintf(log_file, "-----------\n");
                        fprintf(log_file, "!!!DEUX %d\n", n);
                        fprintf(log_file, "m = %d - n = %d\n", m, n);
                        fprintf(log_file, "%d, %d\n", k, l);
                        fprintf(log_file, "%d %d\n", k_index, l_index);
                        fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", d_1, d_2, d_3, d_4, d_5);
                        fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", sub_matrix[m][k_index], sub_matrix[m][l_index], sub_matrix[n][k_index], sub_matrix[n][l_index], sub_matrix[k_index][l_index]);
                        fprintf(log_file, "%lf\n", opti_dist);
                        fprintf(log_file, "-----------\n");
                        
                    }
            }
            
            if(opti_dist != -99.0) break;
            
        }
        
        if(opti_dist != -99.0) break;
        
    }
    fclose(log_file);
    return opti_dist;
    
}

struct tuple distance_estimation_quartet_mean(int m, int n, int number, double sub_matrix[number][number], int threshold)
{

    fflush(stdout);
    
    FILE* log_file = NULL ;
    log_file = fopen("log_imput.txt", "a");
    
    double temp_dist = 0.0;
    double opti_dist = -99.0;
    double total_dist = 0.0;
    int count = 0.0;
    double sum = 0.0;
    double var_dist = -99.0;
    int num_k = 0;
    int num_l = 0;
    int k_index = 0;
    int l_index = 0;
    
    struct mat ik_min_list[number];
    struct mat jl_min_list[number];
    double imputed_values[threshold];
    
    for (int x=0; x<number; x++)
    {
        if(sub_matrix[m][x] != -99.0 && sub_matrix[n][x] != -99.0)
        {
                
            ik_min_list[num_k].value = sub_matrix[m][x];
            ik_min_list[num_k].index = x;
            num_k += 1;
            
        }
    }
    
    for (int x=0; x<number; x++)
    {
            
        if(sub_matrix[n][x] != -99.0 && sub_matrix[m][x] != -99.0)
        {
                
            jl_min_list[num_l].value = sub_matrix[n][x];
            jl_min_list[num_l].index = x;
            num_l +=1;
                
        }
    }
    
    qsort(ik_min_list, num_k, sizeof(ik_min_list[0]), cmp);
    qsort(jl_min_list, num_l, sizeof(jl_min_list[0]), cmp);
    
    for (int k=1; k<num_k; k++)
    {
        k_index = ik_min_list[k].index;
        
        for (int l=1; l<num_l; l++)
        {
            
            fprintf(log_file, "-----------\n");
            fprintf(log_file, "m = %d - n = %d\n", m, n);
            fprintf(log_file, "%d, %d\n", k, l);
            fprintf(log_file, "%d %d\n", k_index, l_index);
            fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", sub_matrix[m][k_index], sub_matrix[m][l_index], sub_matrix[n][k_index], sub_matrix[n][l_index], sub_matrix[k_index][l_index]);
            fprintf(log_file, "-----------\n");
            
            l_index = jl_min_list[l].index;
            
            if(k_index == l_index) continue;
            if(sub_matrix[k_index][l_index] < 0 || sub_matrix[m][k_index] < 0 || sub_matrix[m][l_index] < 0 || sub_matrix[n][k_index] < 0 || sub_matrix[n][l_index] < 0) continue;
            if(sub_matrix[k_index][l_index] < sub_matrix[m][k_index] && sub_matrix[k_index][l_index] < sub_matrix[m][l_index] && sub_matrix[k_index][l_index] < sub_matrix[n][k_index] && sub_matrix[k_index][l_index] < sub_matrix[n][l_index]) continue;
            
            double d_1 = 0.0;
            double d_2 = 0.0;
            double d_3 = 0.0;
            double d_4 = 0.0;
            double d_5 = 0.0;
            
            if(sub_matrix[n][k_index] < sub_matrix[n][l_index] && sub_matrix[m][l_index] < sub_matrix[m][k_index] && sub_matrix[n][k_index] < sub_matrix[k_index][l_index] && sub_matrix[m][l_index] < sub_matrix[k_index][l_index])
            {
                
                if(sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index] > 0)
                {
                    //DIVISER PAR DEUX BORDEL !!!!
                    d_1 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index]);
                    d_2 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[k_index][l_index] - sub_matrix[m][k_index]);
                    d_3 = 0.5 * (sub_matrix[n][l_index] + sub_matrix[n][k_index] - sub_matrix[k_index][l_index]);
                    d_4 = 0.5 * (sub_matrix[n][k_index] + sub_matrix[k_index][l_index] - sub_matrix[n][l_index]);
                    d_5 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[n][l_index] - sub_matrix[m][l_index] - sub_matrix[n][k_index]);
                    
                    if (d_5 < d_1 && d_5 < d_2 && d_5 < d_3 && d_5 < d_4 ) continue;
                    
                    if (d_1 > 0 && d_2 > 0 && d_3 > 0 && d_4 > 0 && d_5 > 0)
                    {
                        
                        total_dist += sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index];
                        imputed_values[count] = sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index];
                        count += 1.0;
                        
                    }
                }
            }
            
            else if(sub_matrix[n][k_index] > sub_matrix[n][l_index] && sub_matrix[m][l_index] > sub_matrix[m][k_index] && sub_matrix[n][l_index] < sub_matrix[k_index][l_index] && sub_matrix[m][k_index] < sub_matrix[k_index][l_index])
            {
                
                if(sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index] > 0)
                {
                    //DIVISER PAR DEUX BORDEL !!!!
                    d_1 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index]);
                    d_2 = 0.5 * (sub_matrix[m][k_index] + sub_matrix[k_index][l_index] - sub_matrix[m][l_index]);
                    d_3 = 0.5 * (sub_matrix[n][k_index] + sub_matrix[n][l_index] - sub_matrix[k_index][l_index]);
                    d_4 = 0.5 * (sub_matrix[n][l_index] + sub_matrix[k_index][l_index] - sub_matrix[n][k_index]);
                    d_5 = 0.5 * (sub_matrix[m][l_index] + sub_matrix[n][k_index] - sub_matrix[m][k_index] - sub_matrix[n][l_index]);
                    
                    if (d_5 < d_1 && d_5 < d_2 && d_5 < d_3 && d_5 < d_4 ) continue;
                    
                    if(d_1 > 0 && d_2 > 0 && d_3 > 0 && d_4 > 0 && d_5 > 0)
                    {
                            
                        total_dist += sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index];
                        imputed_values[count] = sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index];
                        count += 1.0;
                        
                    }
                }
            }
            
            if (count == threshold) break;
            
        }
        
        if (count == threshold) break;
        
    }
    if (total_dist != 0 && count > 1)
    {
        //Mean computation
        opti_dist = total_dist / count;
        
        for (int a=0; a<count; a++)
        {
            sum = sum + pow((imputed_values[a] - opti_dist), 2);
            
        }
        
        var_dist = sum / count;
        
        fprintf(log_file, "-1---------\n");
        fprintf(log_file, "m = %d - n = %d\n", m, n);
        fprintf(log_file, "%lf\n", total_dist);
        fprintf(log_file, "%d\n", count);
        fprintf(log_file, "%lf\n", opti_dist);
        fprintf(log_file, "%lf\n", var_dist);
        
    }
    
    if (total_dist != 0 && count == 1)
    {
        //Mean computation
        opti_dist = total_dist / count;
        var_dist = opti_dist;
        
        fprintf(log_file, "-1---------\n");
        fprintf(log_file, "m = %d - n = %d\n", m, n);
        fprintf(log_file, "%lf\n", total_dist);
        fprintf(log_file, "%d\n", count);
        fprintf(log_file, "%lf\n", opti_dist);
        fprintf(log_file, "%lf\n", var_dist);
        
    }
    
    fclose(log_file);
    
    struct tuple ret = { opti_dist, var_dist };
    return ret;
    
}
