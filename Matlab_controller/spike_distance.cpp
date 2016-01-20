#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/****************************************************************************/
/* Input arguments                                                          */
/****************************************************************************/
#define IN_spike_matrix     prhs[0]	//ST1
#define IN_spike_vector     prhs[1]	//ST2
#define IN_Time_Length		prhs[2]	//T
#define IN_Max_Num_Spikes	prhs[3]	//MaxSpike                       //Don't need this??!!!

/****************************************************************************/
/* Output arguments                                                         */
/****************************************************************************/
#define OUT					plhs[0]

	
	//int M, N, m, n;

	/* Since this is a custome C code, I will not checks for data type, assume user knows the valid type */

	double TimeLength, MaxSpike;
	int i,j, M, N, Total_ROW, Num_rows_single_spike;
	int col, row, left_i, right_i;
	double *output, *spike_matrix, *spike_vector;
	double *Sorted_SpikeTimings, *spiketime_label;  
    double sum, squaredCount;

	TimeLength = mxGetScalar(IN_Time_Length);
	MaxSpike = mxGetScalar(IN_Max_Num_Spikes);
	
	spike_matrix = mxGetPr(IN_spike_matrix);
	spike_vector = mxGetPr(IN_spike_vector);
	
	M = mxGetM(IN_spike_matrix);
	N = mxGetN(IN_spike_matrix);
	Num_rows_single_spike = mxGetM(IN_spike_vector);

	Total_ROW = M + Num_rows_single_spike;

	OUT = mxCreateDoubleMatrix(1, N, mxREAL); /* Create the output matrix */
	output = mxGetPr(OUT);
	//=================================SORTING=================================Merge=================================
	
    //Sorted_SpikeTimings = (double**) malloc(Total_ROW * sizeof(double*));
    //for (i = 0; i < Total_ROW; i++)  
    //Sorted_SpikeTimings[i] = (double*) malloc(N * sizeof(double));
    Sorted_SpikeTimings = (double*) malloc(Total_ROW * sizeof(double));
    

	// Can be cut down to 1D array
    spiketime_label = (double*) malloc(Total_ROW * sizeof(double)); 

	for(col = 0; col < N; col++)
	{
		left_i = 0; right_i = 0; row = 0;

        if(spike_matrix[col*M+left_i] == -1 && spike_vector[right_i] == -1)
            output[col] = 0.0;
        else{
			while(left_i < M || right_i < Num_rows_single_spike){
                //Check for empty spike trains
				if(spike_matrix[col*M+left_i] == -1)
					left_i = M;
                if(spike_vector[right_i] == -1)
					right_i = Num_rows_single_spike;

				if (left_i < M && right_i < Num_rows_single_spike){
					if(spike_matrix[col*M+left_i] <= spike_vector[right_i]){
						Sorted_SpikeTimings[row] = spike_matrix[col*M+left_i];
						spiketime_label[row] = 1;
						left_i++;
						row++;
					}
					else{
						Sorted_SpikeTimings[row] = spike_vector[right_i];
						spiketime_label[row] = -1;
						right_i++;
						row++;
					}
				}

				else if (left_i < M){
                    while(left_i < M && spike_matrix[col*M+left_i] != -1){
						Sorted_SpikeTimings[row] = spike_matrix[col*M+left_i];
						spiketime_label[row] = 1;                        
                        left_i++;
						row++;
                    }
                    break;
				}
				
				else if (right_i < Num_rows_single_spike){
                    while(right_i < Num_rows_single_spike && spike_vector[right_i] != -1){
						Sorted_SpikeTimings[row] = spike_vector[right_i];
						spiketime_label[row] = -1;                        
                        right_i++;
						row++;
                    }
                    break;
				}
			}
            
            //what happens when ROW is < 2 (i.e., = 1)?
            sum = 0.0; squaredCount = spiketime_label[0];
            for(i = 0; i<row-1; i++)
            {
                sum += (Sorted_SpikeTimings[i+1]-Sorted_SpikeTimings[i])*squaredCount*squaredCount;
                squaredCount += spiketime_label[i+1];
            }
            sum += (TimeLength-Sorted_SpikeTimings[row-1])*squaredCount*squaredCount;
            //mexPrintf("Value We Need: %f\n\n",sum);         
            output[col] = sum;
        }
	}
    
    //for (i = 0; i < Total_ROW; i++)      
    //    free(Sorted_SpikeTimings[i]);
    free(Sorted_SpikeTimings);
    
    free(spiketime_label);
    return;
}