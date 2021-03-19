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
            
            if(sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index] > sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index])
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
            
            else if(sub_matrix[n][k_index] + sub_matrix[m][l_index] - sub_matrix[k_index][l_index] > sub_matrix[n][l_index] + sub_matrix[m][k_index] - sub_matrix[k_index][l_index])
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
            else
            {
                fprintf(log_file, "-----------\n");
                fprintf(log_file, "!!!TROIS %d\n", n);
                fprintf(log_file, "m = %d - n = %d\n", m, n);
                fprintf(log_file, "%d, %d\n", k, l);
                fprintf(log_file, "%d %d\n", k_index, l_index);
                fprintf(log_file, "%lf, %lf, %lf, %lf, %lf\n", sub_matrix[m][k_index], sub_matrix[m][l_index], sub_matrix[n][k_index], sub_matrix[n][l_index], sub_matrix[k_index][l_index]);
                fprintf(log_file, "-----------\n");
            }
            
            if(opti_dist != -99.0) break;
            
        }
        
        if(opti_dist != -99.0) break;
        
    }
    fclose(log_file);
    return opti_dist;
    
}
