///////////
// Changes made:
// 1. Added another parameter: number of populations
// 2. Added a hierarchical merging step based on density change between centroids of two hyper-regions
// 3. Picked the longest dimensions between the population centroids to judge whether the two parts should be merged
// 4. Removed checking time parameter
// 5. Output error to stderr
// 6. Fixed the bug of density threshold always = 3
// 7. Added another error (select_num_bin<min_grid) || (select_num_bin>max_grid) to STDERR
// 8. Fixed a bug for 2D data by using K=e*K
// 9. Added some header files, may not be necessary
// 10. Added a lower bound (at least two) for number of populations
/***************************************************************************************************************************************
	
	FLOCK: FLOw cytometry Clustering without K (Named by: Jamie A. Lee and Richard H. Scheuermann) 
	
	Author: (Max) Yu Qian, Ph.D.
	
	Copyright: Scheuermann Lab, Dept. of Pathology, UTSW
	
	Development: November 2005 ~ Forever

	Status: July 2010: Release 2.0

	Usage: flock data_file
		    Note: the input file format must be channel values and the delimiter between two values must be a tab.


    	
****************************************************************************************************************************************/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>



#define DEBUG 0
#define LINE_LEN 1024
#define FILE_NAME_LEN 128
#define PARA_NAME_LEN 64
#define MAX_VALUE 1000000000
#define MIN_GRID 6
#define MAX_GRID 50
#define E_T 1.0

#define NORM_METHOD 2 //2 if z-score; 0 if no normalization; 1 if min-max 
#define KMEANS_TERM 10 
#define MAX_POP_NUM 128

#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

static long **Gr=0;
static long *gcID = 0;          /* grid cluster IDs */
static long *cluster_count=0;   /* count of nodes per cluster */
static long ndense=0;
static long ndim=0;
/* cid changes between depth-first searches, but is constant within a
   single search, so it goes here. */
static long cid=0;
/* Do a depth-first search for a single connected component in graph
 * G.  Start from node, tag the nodes found with cid, and record
 * the tags in grid_clusterID.  Also, record the node count in
 * cluster_count.  If we find a node that has already been assigned to
 * a cluster, that means we're merging two clusters, so zero out the
 * old cid's node count.
 *
 * Note that our graph is constructed as a DAG, so we can be
 * guaranteed to terminate without checking for cycles.  
 *
 * Note2: this function can potentially recurse to depth = ndense.
 * Plan stack size accordingly.  
 *
 * Output:
 *
 * grid_clusterID[] -- array where we tag the nodes.
 * cluster_count[]  -- count of the number of nodes per cluster.
 */
static void merge_cluster(long from, long into)
{
  int i;

  for(i=0; i<ndense; ++i)
    if(gcID[i] == from)
      gcID[i] = into;
}
   
void depth_first(long node)
{
  long i;

  if(gcID[node] == cid)         // we're guaranteed no cycles, but it is possible to reach a node 
    return;                     // through two different paths in the same cluster.  This early
                                // return saves us some unneeded work.

  /* Check to see if we are merging a cluster */
  if(gcID[node] >= 0) {
    /* We are, so zero the count for the old cluster. */
    cluster_count[ gcID[node] ] = 0;  
    merge_cluster(gcID[node], cid);
    return;
  }

  /* Update for this node */
  gcID[node] = cid;
  cluster_count[cid]++;

  /* Recursively search the child nodes */
  for(i=0; i<ndim; ++i)
    if(Gr[node][i] >= 0)      /* This is a child node */
      depth_first(Gr[node][i]);
}

void bail(const char *msg)
{
  fprintf(stderr,"%s",msg);
  exit(0);
}


static void check_clusters(long *gcID, int ndense)
{
  int i;

  for(i=0; i<ndense; ++i)
    if(gcID[i] < 0) {
      fprintf(stderr,"faulty cluster id at i= %d\n",i);
      exit(0);
    }
}



long find_connected(long **G, long n_dense_grids, long num_dm, long *grid_clusterID)
{
  long nclust=0;                  /* number of clusters found */
  long i;
  long *subfac;
  long clustid=0;
  int subval=0,nempty=0;
  
  int sz = n_dense_grids*sizeof(long); 
  cluster_count = malloc(sz);
  if(!cluster_count)
    bail("find_connected: Unable to allocate %zd bytes.\n");
  memset(cluster_count,0,sz);

  /* set up the statics that will be used in the DFS */
  Gr=G;
  gcID = grid_clusterID;
  ndense = n_dense_grids;
  ndim = num_dm;
  
  for(i=0;i<ndense;++i)
    grid_clusterID[i] = -1;

  for(i=0;i<ndense;++i) {
    if(grid_clusterID[i] < 0) {  /* grid hasn't been assigned yet */
      cid = nclust++;
      depth_first(i);
    }
  }


  
  /* At this point we probably have some clusters that are empty due to merging.
     We want to compact the cluster numbering to eliminate the empty clusters. */

  subfac = malloc(sz);
  if(!subfac)
    bail("find_connected: Unable to allocate %zd bytes.\n");
  subval=0;
  nempty=0;

  /* cluster #i needs to have its ID decremented by 1 for each empty cluster with
     ID < i.  Precaclulate the decrements in this loop: */
  for(i=0;i<nclust;++i) {
    //clustid = grid_clusterID[i];
    if(cluster_count[i] == 0) {
      subval++;
      nempty++;
    }
    subfac[i] = subval;
  }

  //printf("nempty is %d\n",nempty);
  
  /* Now apply the decrements to all of the dense grids */
  for(i=0;i<ndense;++i) {
    clustid = grid_clusterID[i];
    grid_clusterID[i] -= subfac[clustid];
  }


  
  /* correct the number of clusters found */
  nclust -= nempty;

  //printf("nclust is %d\n",nclust);

  return nclust;
}

/************************************* Read basic info of the source file **************************************/
/************************************* Read basic info of the source file **************************************/
void getfileinfo(FILE *f_src, long *file_Len, long *num_dm, char *name_string, long *time_ID)
{
  char src[LINE_LEN];
  char current_name[64];
  char prv;

  long num_rows=0;
  long num_columns=0;
  long ch='\n';
  long prev='\n';
  long time_pos=0;
  long i=0;
  long j=0;
  int sw=0;

  src[0]='\0';
  fgets(src, LINE_LEN, f_src);

  name_string[0]='\0';
  current_name[0]='\0';
  prv='\n';

  while ((src[i]==' ') || (src[i]=='\t')) //skip space and tab characters
    i++;

  while ((src[i]!='\0') && (src[i]!='\n')) //repeat until the end of the line
    {
      current_name[j]=src[i];
		
      if ((src[i]=='\t') && (prv!='\t')) //a complete word
        {
          current_name[j]='\0';

          if (0!=strcmp(current_name,"Time"))
            {
              num_columns++; //num_columns does not inlcude the column of Time
              time_pos++;
              if (sw) {
                  strcat(name_string,"\t");
              }
              strcat(name_string,current_name); 
              sw = 1;
            }
          else
            {
              *time_ID=time_pos;
            }

          current_name[0]='\0';
          j=0;			
        }		
		
      if ((src[i]=='\t') && (prv=='\t')) //a duplicate tab or space
        {
          current_name[0]='\0';
          j=0;
        }
		
      if (src[i]!='\t')
        j++;
		
      prv=src[i];
      i++;
    }
	
  //name_string[j]='\0';
  if (prv!='\t') //the last one hasn't been retrieved
    {
      current_name[j]='\0';
      if (0!=strcmp(current_name,"Time"))
        {
          num_columns++;
          strcat(name_string,"\t");
          strcat(name_string,current_name);
          time_pos++;
        }
      else
        {
          *time_ID=time_pos;
        }
    }
  if (DEBUG==1)
    {
      printf("time_ID is %ld\n",*time_ID);
      printf("name_string is %s\n",name_string);
    }

  //start computing # of rows

  while ((ch = fgetc(f_src))!= EOF )
    {
      if (ch == '\n')
        {
          ++num_rows;
        }
      prev = ch;
    }
  if (prev!='\n')
    ++num_rows;

   if (num_rows<50)
  {
    fprintf(stderr,"Number of events in the input file is too few and should not be processed!\n"); //modified on July 23, 2010
	exit(0);
  }
	
  *file_Len=num_rows;
  *num_dm=num_columns; 

  printf("original file size is %ld; number of dimensions is %ld\n", *file_Len, *num_dm);
}



/************************************* Read the source file into uncomp_data **************************************/
void readsource(FILE *f_src, long file_Len, long num_dm, double **uncomp_data, long time_ID)
{
  //long time_pass=0; //to mark whether the time_ID has been passed
  long index=0;

  long i=0;
  long j=0;
  long t=0;

  char src[LINE_LEN];
  char xc[LINE_LEN/10];

  src[0]='\0';
  fgets(src,LINE_LEN, f_src); //skip the first line about parameter names

  while (!feof(f_src) && (index<file_Len)) //index = 0, 1, ..., file_Len-1
    {
      src[0]='\0';	    
      fgets(src,LINE_LEN,f_src);
      i=0;
      //time_pass=0;
						
      //if (time_ID==-1)  //there is no time_ID
      //  {
          for (t=0;t<num_dm;t++) 
            {
              xc[0]='\0';
              j=0;
              while ((src[i]!='\0') && (src[i]!='\n') && (src[i]!=' ') && (src[i]!='\t'))
                {
                  xc[j]=src[i];
                  i++;
                  j++;
                }
		
              xc[j]='\0';	    
              i++;

              uncomp_data[index][t]=atof(xc);
            }	
       // }
      /*else
        {
          for (t=0;t<=num_dm;t++) //the time column needs to be skipped, so there are num_dm+1 columns
            {
              xc[0]='\0';
              j=0;
              while ((src[i]!='\0') && (src[i]!='\n') && (src[i]!=' ') && (src[i]!='\t'))
                {
                  xc[j]=src[i];
                  i++;
                  j++;
                }
		
              xc[j]='\0';	    
              i++;

              if (t==time_ID)
                {
                  time_pass=1;
                  continue;
                }
				
              if (time_pass)
                uncomp_data[index][t-1]=atof(xc);
              else
                uncomp_data[index][t]=atof(xc);
            }
        }*/ //commented by Yu Qian on Aug 31, 2010 as we do not want to check time_ID anymore        	
      index++;     	
      //fprintf(fout_ID,"%s",src);
    } //end of while
	
  if (DEBUG == 1)
    {
      printf("the last line of the source data is:\n");
      for (j=0;j<num_dm;j++)
        printf("%f ",uncomp_data[index-1][j]);
      printf("\n");
    }
}


/**************************************** Normalization ******************************************/
void tran(double **orig_data, long file_Len, long num_dm, long norm_used, double **matrix_to_cluster)
{
  long i=0;
  long j=0;

  double biggest=0;
  double smallest=MAX_VALUE;

  double *aver; //average of each column
  double *std; //standard deviation of each column

  aver=(double*)malloc(sizeof(double)*file_Len);
  memset(aver,0,sizeof(double)*file_Len);

  std=(double*)malloc(sizeof(double)*file_Len);
  memset(std,0,sizeof(double)*file_Len);	
		
  if (norm_used==2) //z-score normalization
    {
      for (j=0;j<num_dm;j++)
        {
          aver[j]=0;
          for (i=0;i<file_Len;i++)
            aver[j]=aver[j]+orig_data[i][j];
          aver[j]=aver[j]/(double)file_Len;

          std[j]=0;
          for (i=0;i<file_Len;i++)
            std[j]=std[j]+((orig_data[i][j]-aver[j])*(orig_data[i][j]-aver[j]));
          std[j]=sqrt(std[j]/(double)file_Len);
			
          for (i=0;i<file_Len;i++)
            matrix_to_cluster[i][j]=(orig_data[i][j]-aver[j])/std[j];  //z-score normalization
        }
    }

  if (norm_used==1) //0-1 min-max normalization
    {
      for (j=0;j<num_dm;j++)
        {
          biggest=0;
          smallest=MAX_VALUE;
          for (i=0;i<file_Len;i++)
            {
              if (orig_data[i][j]>biggest)
                biggest=orig_data[i][j];
              if (orig_data[i][j]<smallest)
                smallest=orig_data[i][j];
            }
			
          for (i=0;i<file_Len;i++)
            {
              if (biggest==smallest)
                matrix_to_cluster[i][j]=biggest;
              else
                matrix_to_cluster[i][j]=(orig_data[i][j]-smallest)/(biggest-smallest);
            }
        }
    }

  if (norm_used==0) //no normalization
    {
      for (i=0;i<file_Len;i++)
        for (j=0;j<num_dm;j++)
          matrix_to_cluster[i][j]=orig_data[i][j];
    }

}



/********************************************** RadixSort *******************************************/
/* Perform a radix sort using each dimension from the original data as a radix.
 * Outputs:
 * sorted_seq   -- a permutation vector mapping the ordered list onto the original data.
 *                  (sorted_seq[i] -> index in the original data of the ith element of the ordered list)
 * grid_ID      -- mapping between the original data and the "grids" (see below) found as a byproduct
 *                  of the sorting procedure.
 * num_nonempty -- the number of grids that occur in the data (= the number of distinct values assumed
 *                  by grid_ID)
 */

void radixsort_flock(long **position,long file_Len,long num_dm,long num_bin,long *sorted_seq,long *num_nonempty,long *grid_ID)
{
  long i=0;
  long length=0; 
  long start=0;
  long prev_ID=0;
  long curr_ID=0;
	
  long j=0;
  long t=0;
  long p=0;
  long loc=0;
  long temp=0;
  long equal=0;
	
  long *count; //count[i]=j means there are j numbers having value i at the processing digit
  long *index; //index[i]=j means the starting position of grid i is j
  long *cp; //current position
  long *mark; //mark[i]=0 means it is not an ending point of a part, 1 means it is (a "part" is a group of items with identical bins for all dimensions)
  long *seq; //temporary sequence

  count=(long*)malloc(sizeof(long)*num_bin);
  memset(count,0,sizeof(long)*num_bin);

  cp=(long*)malloc(sizeof(long)*num_bin);
  memset(cp,0,sizeof(long)*num_bin);

  index=(long*)malloc(sizeof(long)*num_bin); // initialized below

  seq=(long*)malloc(sizeof(long)*file_Len);
  memset(seq,0,sizeof(long)*file_Len);
	
  mark=(long*)malloc(sizeof(long)*file_Len);
  memset(mark,0,sizeof(long)*file_Len);
	
  for (i=0;i<file_Len;i++)
    {
      sorted_seq[i]=i;
      mark[i]=0;
      seq[i]=0;
    }
  for (i=0;i<num_bin;i++)
    {
      index[i]=0;
      cp[i]=0;
      count[i]=0;
    }

  for (j=0;j<num_dm;j++)
    {
      if (j==0) //compute the initial values of mark
        {
          for (i=0;i<file_Len;i++)
            count[position[i][j]]++; // initialize the count to the number of items in each bin of the 0th dimension

          index[0] = 0;
          for (i=0;i<num_bin-1;i++)
            {
              index[i+1]=index[i]+count[i];  //index[k]=x means k segment starts at x (in the ordered list)
              if ((index[i+1]>0) && (index[i+1]<=file_Len))
                {
                  mark[index[i+1]-1]=1; // Mark the end of the segment in the ordered list
                }
              else
                {
                  printf("out of myboundary for mark at index[i+1]-1.\n");
                }
            }
          mark[file_Len-1]=1;
			
          for (i=0;i<file_Len;i++)
            {
              /* Build a permutation vector for the partially ordered data.  Store the PV in sorted_seq */
              loc=position[i][j];
              temp=index[loc]+cp[loc]; //cp[i]=j means the offset from the starting position of grid i is j 
              sorted_seq[temp]=i;  //sorted_seq[i]=temp is also another way to sort
              cp[loc]++;
            }
        }
      else
        {
          //reset count, index, loc, temp, cp, start, and length
          length=0;
          loc=0;
          temp=0;
          start=0;
          for (p=0;p<num_bin;p++)
            {
              cp[p]=0;
              count[p]=0;
              index[p]=0;
            }

          for (i=0;i<file_Len;i++)
            {
              long iperm = sorted_seq[i]; // iperm allows us to traverse the data in sorted order.
              if (mark[i]!=1)
                {
                  /* Count the number of items in each bin of
                     dimension j, BUT we are going to reset at the end
                     of each "part".  Thus, the total effect is to do
                     a sort by bin on the jth dimension for each group
                     of data that has been identical for the
                     dimensions processed up to this point.  This is
                     the standard radix sort procedure, but doing it
                     this way saves us having to allocate buckets to
                     hold the data in each group of "identical-so-far"
                     elements. */
                  count[position[iperm][j]]++;  //count[position[i][j]]++;
                  length++;                     // This is the total length of the part, irrespective of the value of the jth component
                                                // (actually, less one, since we don't increment for the final element below)
                }
              if (mark[i]==1)
                {
                  //length++;
                  count[position[iperm][j]]++;//count[position[i][j]]++;  //the current point must be counted in
                  start=i-length; //this part starts from start to i: [start,i]
                  /* Now we sort on the jth radix, just like we did for the 0th above, but we restrict it to just the current part.
                     This would be a lot more clear if we broke this bit of code out into a separate function and processed recursively,
                     plus we could multi-thread over the parts.  (Hmmm...)
                  */
                  index[0] = start; // Let index give the offset within the whole ordered list.
                  for (t=0;t<num_bin-1;t++)
                    {
                      index[t+1]=index[t]+count[t];
						
                      if ((index[t+1]<=file_Len) && (index[t+1]>0))
                        {
                          mark[index[t+1]-1]=1; // update the part boundaries to include the differences in the current radix.
                        }
						
                    }
                  mark[i]=1;

                  /* Update the permutation vector for the current part (i.e., from start to i).  By the time we finish the loop over i
                     the PV will be completely updated for the partial ordering up to the current radix. */
                  for (t=start;t<=i;t++)
                    {
                      loc=position[sorted_seq[t]][j];//loc=position[t][j];
                      temp=index[loc]+cp[loc];
                      if ((temp<file_Len) && (temp>=0)) 
                        {
                          // seq is a temporary because we have to maintain the old PV until we have finished this step.
                          seq[temp]=sorted_seq[t];  //sorted_seq[i]=temp is also another way to sort
                          cp[loc]++;
                        }
                      else
                        {
                          printf("out of myboundary for seq at temp.\n");
                        }
                    }

                  for (t=start;t<=i;t++)
                    {
                      // copy the temporary back into sorted_seq.  sorted_seq is now updated for radix j up through
                      // entry i in the ordered list.
                      if ((t>=0) && (t<file_Len))
                        sorted_seq[t]=seq[t];
                      else
                        printf("out of myboundary for seq and sorted_seq at t.\n");
                    }
                  //reset count, index, seq, length, and cp
                  length=0;
                  loc=0;
                  temp=0;
                  for (p=0;p<num_bin;p++)
                    {
                      cp[p]=0;
                      count[p]=0;
                      index[p]=0;
                    }
                }
            }//end for i
        }//end else
    }//end for j

  /* sorted_seq[] now contains the ordered list for all radices.  mark[] gives the boundaries between groups of elements that are
     identical over all radices (= dimensions in the original data) (although it appears we aren't going to make use of this fact) */
  
  for (i=0;i<file_Len;i++)
    grid_ID[i]=0; //in case the initial value hasn't been assigned
  *num_nonempty=1; //starting from 1!	

  /* assign the "grid" identifiers for all of the data.  A grid will be what we were calling a "part" above.  We will number them
     serially and tag the *unordered* data with the grid IDs.  We will also count the number of populated grids (in general there will
     be many possible combinations of bin values that simply never occur) */
  
  for (i=1;i<file_Len;i++)
    {
      equal=1;
      prev_ID=sorted_seq[i-1];
      curr_ID=sorted_seq[i];
      for (j=0;j<num_dm;j++)
        {
          if (position[prev_ID][j]!=position[curr_ID][j])
            {	
              equal=0;  //not equal
              break;
            }
        }
		
      if (equal)
        {
          grid_ID[curr_ID]=grid_ID[prev_ID];
        }
      else
        {
          *num_nonempty=*num_nonempty+1;
          grid_ID[curr_ID]=grid_ID[prev_ID]+1;
        }
      //all_grid_vol[grid_ID[curr_ID]]++;
    }

  //free memory
  free(count);
  free(index);	
  free(cp);	
  free(seq);
  free(mark); 
  
}

/********************************************** Compute Position of Events ************************************************/
void compute_position(double **data_in, long file_Len, long num_dm, long num_bin, long **position, double *interval)
{
  /* What we are really doing here is binning the data, with the bins
     spanning the range of the data and number of bins = num_bin */
  long i=0;
  long j=0;

  double *small; //small[j] is the smallest value within dimension j
  double *big; //big[j] is the biggest value within dimension j
		
  small=(double*)malloc(sizeof(double)*num_dm);
  memset(small,0,sizeof(double)*num_dm);

  big=(double*)malloc(sizeof(double)*num_dm);
  memset(big,0,sizeof(double)*num_dm);
	
	
  for (j=0;j<num_dm;j++)
    {
      big[j]=MAX_VALUE*(-1);
      small[j]=MAX_VALUE;
      for (i=0;i<file_Len;i++)
        {
          if (data_in[i][j]>big[j])
            big[j]=data_in[i][j];

          if (data_in[i][j]<small[j])
            small[j]=data_in[i][j];
        }
		
      interval[j]=(big[j]-small[j])/(double)num_bin;	//interval is computed using the biggest value and smallest value instead of the channel limit
      /* XXX: I'm pretty sure the denominator of the fraction above should be num_bin-1. */
    }
    
  for (j=0;j<num_dm;j++)
  {
     for (i=0;i<file_Len;i++)
     {	
        if (data_in[i][j]>=big[j])
           position[i][j]=num_bin-1;
        else
        {
           position[i][j]=(long)((data_in[i][j]-small[j])/interval[j]); //position[i][j]=t means point i is at the t grid of dimensional j
           if ((position[i][j]>=num_bin) || (position[i][j]<0))
           {
               //printf("position mis-computed in density analysis!\n");
               //exit(0);
			   fprintf(stderr,"Incorrect input file format or input parameters (number of bins overflows)!\n"); //modified on July 23, 2010
				exit(0);
           }
        }
     }
  }


  free(small);
  free(big);
}

/********************************************** select_bin to select the number of bins **********************************/
//num_bin=select_bin(normalized_data, file_Len, num_dm, MIN_GRID, MAX_GRID, position, sorted_seq, all_grid_ID, &num_nonempty);
/* Determine the number of bins to use in each dimension.  Additionally sort the data elements according to the binned
 * values, and partition the data into "grids" with identical (binned) values.  We try progressively more bins until we
 * maximize a merit function, then return the results obtained using the optimal number of bins. 
 *
 * Outputs:
 * position     -- binned data values
 * sorted_seq   -- permutation vector mapping the ordered list to the original data
 * all_grid_ID  -- grid to which each data element was assigned.
 * num_nonempty -- number of distinct values assumed by all_grid_ID
 * interval     -- bin width for each data dimension
 * return value -- the number of bins selected.
 */

long select_bin(double **normalized_data, long file_Len, long num_dm, long min_grid, long max_grid, long **position, long *sorted_seq, 
                long *all_grid_ID, long *num_nonempty, double *interval, long user_num_bin)
{
 
  long num_bin=0;
  long select_num_bin=0;
  long m=0;
  long n=0;
	
  long i=0;
  long bin_scope=0;
  long temp_num_nonempty=0;

  long *temp_grid_ID;
  long *temp_sorted_seq;
  long **temp_position;

  //sorted_seq[i]=j means the event j ranks i

  double temp_index=0;
  double *bin_index;	
  double *temp_interval;
	
	
  temp_grid_ID=(long *)malloc(sizeof(long)*file_Len);
  memset(temp_grid_ID,0,sizeof(long)*file_Len);
	
  temp_sorted_seq=(long *)malloc(sizeof(long)*file_Len);
  memset(temp_sorted_seq,0,sizeof(long)*file_Len);

  temp_position=(long **)malloc(sizeof(long*)*file_Len);
  memset(temp_position,0,sizeof(long*)*file_Len);
  for (m=0;m<file_Len;m++)
    {
      temp_position[m]=(long*)malloc(sizeof(long)*num_dm);
      memset(temp_position[m],0,sizeof(long)*num_dm);
    }

  temp_interval=(double*)malloc(sizeof(double)*num_dm);
  memset(temp_interval,0,sizeof(double)*num_dm);

  bin_scope=max_grid-min_grid+1;
  bin_index=(double *)malloc(sizeof(double)*bin_scope);
  memset(bin_index,0,sizeof(double)*bin_scope);

  i=0;

  for (num_bin=min_grid;num_bin<=max_grid;num_bin++)
    {
      
	  temp_num_nonempty=0;
	   /* compute_position bins the data into num_bin bins.  Each
         dimension is binned independently.

         Outputs:
         temp_position[i][j] -- bin for the jth component of data element i.
         temp_interval[j]    -- bin-width for the jth component
      */
      compute_position(normalized_data, file_Len, num_dm, num_bin, temp_position, temp_interval);
      radixsort_flock(temp_position,file_Len,num_dm,num_bin,temp_sorted_seq,&temp_num_nonempty,temp_grid_ID);

      /* our figure of merit is the number of non-empty grids divided by number of bins per dimension.
         We declare victory when we have found a local maximum */
      bin_index[i]=((double)temp_num_nonempty)/((double)num_bin);
	  if ((double)(temp_num_nonempty)>=(double)(file_Len)*0.95)
		  break;
      if ((bin_index[i]<temp_index) && (user_num_bin==0))
         break;
	  if ((user_num_bin==num_bin-1) && (user_num_bin!=0))
		 break;
      
      /* Since we have accepted this trial bin, copy all the temporary results into
         the output buffers */
      memcpy(all_grid_ID,temp_grid_ID,sizeof(long)*file_Len);
      memcpy(sorted_seq,temp_sorted_seq,sizeof(long)*file_Len);
      memcpy(interval,temp_interval,sizeof(double)*num_dm);
		
      for (m=0;m<file_Len;m++)
        for (n=0;n<num_dm;n++)
          position[m][n]=temp_position[m][n];

      temp_index=bin_index[i];
      select_num_bin=num_bin;
      num_nonempty[0]=temp_num_nonempty;
      i++;
    }

   if ((select_num_bin<min_grid) || (select_num_bin>max_grid))
  {
    fprintf(stderr,"Number of events collected is too few in terms of number of markers used. The file should not be processed!\n"); //modified on Nov 04, 2010
	exit(0);
  }
	
  if (temp_index==0)
  {
	 fprintf(stderr,"Too many dimensions with too few events in the input file, or a too large number of bins used.\n"); //modified on July 23, 2010
	 exit(0);
  }

  free(temp_grid_ID);
  free(temp_sorted_seq);
  free(bin_index);
  free(temp_interval);
	
  for (m=0;m<file_Len;m++)
    free(temp_position[m]);
  free(temp_position);

  return select_num_bin; 
}

/************************************* Select dense grids **********************************/
// compute num_dense_grids, num_dense_events, dense_grid_reverse, and all_grid_vol
// den_cutoff=select_dense(file_Len, all_grid_ID, num_nonempty, &num_dense_grids, &num_dense_events, dense_grid_reverse);
/*
 * Prune away grids that are insufficiently "dense" (i.e., contain too few data items)
 *
 * Outputs:
 * num_dense_grids    -- number of dense grids
 * num_dense_events   -- total number of data items in all dense grids
 * dense_grid_reverse -- mapping from list of all grids to list of dense grids.
 * return value       -- density cutoff for separating dense from non-dense grids.
 */

int select_dense(long file_Len, long *all_grid_ID, long num_nonempty, long *num_dense_grids, long *num_dense_events, long *dense_grid_reverse, int den_t_event)
{
  

  long i=0;
  long vol_ID=0;
  long biggest_size=0; //biggest grid_size, to define grid_size_index
  long biggest_index=0;
  //long actual_threshold=0; //the actual threshold on grid_size, e.g., 1 implies 1 event in the grid
  //long num_remain=0; //number of remaining grids with different density thresholds
  long temp_num_dense_grids=0;
  long temp_num_dense_events=0;
	
  long *grid_size_index;
  long *all_grid_vol;
  long *grid_density_index;

  //double den_average=0;
 // double avr_index=0;
 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute all_grid_vol
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  all_grid_vol=(long *)malloc(sizeof(long)*num_nonempty);
  memset(all_grid_vol,0,sizeof(long)*num_nonempty);

  /* Grid "volume" is just the number of data contained in the grid. */
  for (i=0;i<file_Len;i++)
    {
      vol_ID=all_grid_ID[i]; //vol_ID=all_grid_ID[sorted_seq[i]];
      all_grid_vol[vol_ID]++;  //all_grid_vol[i]=j means grid i has j points
    }

 
	
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute grid_size_index (histogram of grid sizes)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (i=0;i<num_nonempty;i++)
    {
      if (biggest_size<all_grid_vol[i])
        {
          biggest_size=all_grid_vol[i];
        }
    }

   if (biggest_size<3)
  {
	 fprintf(stderr,"Too many dimensions with too few events in the input file, or a too large number of bins used.\n"); //modified on July 23, 2010
	 exit(0);
  }


  grid_size_index=(long*)malloc(sizeof(long)*biggest_size);
  memset(grid_size_index,0,sizeof(long)*biggest_size);
	
  for (i=0;i<num_nonempty;i++)
    {
      grid_size_index[all_grid_vol[i]-1]++; //grid_size_index[0]=5 means there are 5 grids having size 1
    }



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute den_cutoff
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  grid_density_index=(long *)malloc(sizeof(long)*(biggest_size-2));//from event 2 to biggest_size-1, i.e., from 0 to biggest_size-3
  memset(grid_density_index,0,sizeof(long)*(biggest_size-2));

  if (den_t_event==0)
  {
	  biggest_index=0;

	  for (i=2;i<biggest_size-1;i++) //the grid with 1 event will be skipped, i.e., grid_density_index[0] won't be defined
	  {
		  grid_density_index[i-1]=(grid_size_index[i-1]+grid_size_index[i+1]-2*grid_size_index[i]);
		  if (biggest_index<grid_density_index[i-1])
		  {
			biggest_index=grid_density_index[i-1];
			den_t_event=i+1;
		  }
	  }
  }

  if (den_t_event==0) //if biggest_size==3
  {
	 fprintf(stderr,"the densest hyperregion has only 3 events, smaller than the user-specified value. Therefore the density threshold is automatically changed to 3.\n");
	 den_t_event=3;
  }

  for (i=0;i<num_nonempty;i++)
	  if (all_grid_vol[i]>=den_t_event)
		temp_num_dense_grids++;

  if (temp_num_dense_grids<=1)
  {
	  //exit(0);
	  //printf("a too high density threshold is set! Please decrease your density threshold.\n");
	  fprintf(stderr,"a too high density threshold! Please decrease your density threshold.\n"); //modified on July 23, 2010
	  exit(0);
  }

   if (temp_num_dense_grids>=100000)
  {
	  //modified on July 23, 2010
	  //printf("a too low density threshold is set! Please increase your density threshold.\n");
	  //exit(0);
	  fprintf(stderr,"a too low density threshold! Please increase your density threshold.\n"); //modified on July 23, 2010
	  exit(0);
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute dense_grid_reverse
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  temp_num_dense_grids=0;
  temp_num_dense_events=0;
  for (i=0;i<num_nonempty;i++)
    {
      dense_grid_reverse[i]=-1;
		
      if (all_grid_vol[i]>=den_t_event) 
        {			
          dense_grid_reverse[i]=temp_num_dense_grids;  //dense_grid_reverse provides a mapping from all nonempty grids to dense grids.
          temp_num_dense_grids++;
          temp_num_dense_events+=all_grid_vol[i];						
        }
    }

  num_dense_grids[0]=temp_num_dense_grids;
  num_dense_events[0]=temp_num_dense_events;	

  free(grid_size_index);
  free(all_grid_vol);

  return den_t_event;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute densegridID_To_gridevent and eventID_To_denseventID
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	grid_To_event(file_Len, dense_grid_reverse, all_grid_ID, eventID_To_denseventID, densegridID_To_gridevent);
/*
 * Filter the data so that only the data belonging to dense grids is left
 *
 * Output:
 * eventID_To_denseeventID   -- mapping from original event ID to ID in list containing only events contained in dense grids.
 * densegridID_To_gridevent  -- mapping from dense grids to prototype members of the grids.
 *
 */

void grid_To_event(long file_Len, long *dense_grid_reverse, long *all_grid_ID, long *eventID_To_denseventID, long *densegridID_To_gridevent)
{
  long i=0;
  long dense_grid_ID=0;
  long dense_event_ID=0;

  for (i=0;i<file_Len;i++)
    {
      dense_grid_ID=dense_grid_reverse[all_grid_ID[i]];
      eventID_To_denseventID[i]=-1;
      if (dense_grid_ID!=-1) //for point (i) belonging to dense grids 
        {			
          eventID_To_denseventID[i]=dense_event_ID;
          dense_event_ID++;
		
          if (densegridID_To_gridevent[dense_grid_ID]==-1) //for point i that hasn't been selected
            densegridID_To_gridevent[dense_grid_ID]=i; //densegridID_To_gridevent maps dense_grid_ID to its representative point			
        }		
    }
	
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute dense_grid_seq
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	generate_grid_seq(file_Len, num_dm, sorted_seq, num_dense_grids, densegridID_To_gridevent, position, dense_grid_rank, dense_grid_seq);
/* Construct a table of binned data values for each dense grid.
 *
 * Output:
 *
 * dense_grid_seq  -- table of binned data values for each dense grid (recall that all members of a grid have identical binned data values)
 */

void generate_grid_seq(long num_dm, long num_dense_grids, long *densegridID_To_gridevent, long **position, long **dense_grid_seq)
{  

  long i=0;
  long j=0;
  long ReventID=0; //representative event ID of the dense grid

  for (i=0;i<num_dense_grids;i++)
    {
      ReventID = densegridID_To_gridevent[i];

      for (j=0;j<num_dm;j++)
        dense_grid_seq[i][j]=position[ReventID][j];
    }
}

//compare two vectors
long compare_value(long num_dm, long *search_value, long *seq_value)
{
  long i=0;

  for (i=0;i<num_dm;i++)
    {
      if (search_value[i]<seq_value[i])
        return 1;
      if (search_value[i]>seq_value[i])
        return -1;
      if (search_value[i]==seq_value[i])
        continue;
    }
  return 0;
}

//binary search the dense_grid_seq to return the dense grid ID if it exists
long binary_search(long num_dense_grids, long num_dm, long *search_value, long **dense_grid_seq)
{

  long low=0;
  long high=0;
  long mid=0;
  long comp_result=0;
  long match=0;
  //long found=0;
	
  low = 0;
  high = num_dense_grids-1;
    
  while (low <= high) 
    {
      mid = (long)((low + high)/2);
	
      comp_result=compare_value(num_dm, search_value,dense_grid_seq[mid]);
	
		
      switch(comp_result)
        {
        case 1:
          high=mid-1;
          break;
        case -1:
          low=mid+1;
          break;
        case 0:
          match=1;
          break;
        }
      if (match==1)
        break;
    }
	


  if (match==1)
    {
      return mid;
    }
  else
    return -1;   
}


/********************************************** Computing Centers Using IDs **********************************************/

void ID2Center(double **data_in, long file_Len, long num_dm, long *eventID_To_denseventID, long num_clust, long *cluster_ID, double **population_center)
{
  long i=0;
  long j=0;
  long ID=0;
  long eventID=0;
  long *size_c;



  size_c=(long *)malloc(sizeof(long)*num_clust);
  memset(size_c,0,sizeof(long)*num_clust);

  for (i=0;i<num_clust;i++)
    for (j=0;j<num_dm;j++)
      population_center[i][j]=0;

  for (i=0;i<file_Len;i++)
    {
      eventID=eventID_To_denseventID[i];

      if (eventID!=-1) //only events in dense grids count
        {
          ID=cluster_ID[eventID];
			
          if (ID==-1)
            {
              //printf("ID==-1! in ID2Center\n");
              //exit(0);
			  fprintf(stderr,"Incorrect file format or input parameters (no dense regions found!)\n"); //modified on July 23, 2010
			  exit(0);
            }

          for (j=0;j<num_dm;j++)
            population_center[ID][j]=population_center[ID][j]+data_in[i][j];
			
          size_c[ID]++;				
        }
    }
	

  for (i=0;i<num_clust;i++)
    {
      for (j=0;j<num_dm;j++)
        if (size_c[i]!=0)
          population_center[i][j]=(population_center[i][j]/(double)(size_c[i]));
        else
		{
          population_center[i][j]=0;
		  printf("size_c[%ld]=0 in ID2center\n",i);
		}
    }

  free(size_c);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Compute Population Center with all events
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ID2Center_all(double **data_in, long file_Len, long num_dm, long num_clust, long *cluster_ID, double **population_center)
{
  long i=0;
  long j=0;
  long ID=0;
  long *size_c;



  size_c=(long *)malloc(sizeof(long)*num_clust);
  memset(size_c,0,sizeof(long)*num_clust);

  for (i=0;i<num_clust;i++)
    for (j=0;j<num_dm;j++)
      population_center[i][j]=0;

	
  for (i=0;i<file_Len;i++)
    {
         ID=cluster_ID[i];
			
         if (ID==-1)
         {
            //printf("ID==-1! in ID2Center_all\n");
            //exit(0);
			fprintf(stderr,"Incorrect file format or input parameters (resulting in incorrect population IDs)\n"); //modified on July 23, 2010
			exit(0);
         }

         for (j=0;j<num_dm;j++)
           population_center[ID][j]=population_center[ID][j]+data_in[i][j];
			
         size_c[ID]++;        
    }
	
 
  for (i=0;i<num_clust;i++)
    {
      for (j=0;j<num_dm;j++)
        if (size_c[i]!=0)
          population_center[i][j]=(population_center[i][j]/(double)(size_c[i]));
        else
		{
          population_center[i][j]=0;
		  printf("size_c[%ld]=0 in ID2center_all\n",i);
		}
    }


  free(size_c);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Merge neighboring grids to clusters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

long merge_grids(double **normalized_data, double *interval, long file_Len, long num_dm, long num_bin, long **position, long num_dense_grids, 
                 long *dense_grid_reverse, long **dense_grid_seq, long *eventID_To_denseventID, long *densegridID_To_gridevent, long *all_grid_ID,
                 long *cluster_ID, long *grid_ID, long *grid_clusterID)
{
  

  long i=0;
  long j=0;
  long t=0;
  long p=0;
  long num_clust=0;
  long ReventID=0;
  long denseID=0;
  long neighbor_ID=0;
  //long temp_grid=0;

  long *grid_value;
  long *search_value;
	
  long **graph_of_grid; //the graph constructed for dense grids: each dense grid is a graph node

  double real_dist=0;
  double **norm_grid_center;
	
  norm_grid_center=(double **)malloc(sizeof(double*)*num_dense_grids);
  memset(norm_grid_center,0,sizeof(double*)*num_dense_grids);
	
  for (i=0;i<num_dense_grids;i++)
    {
      norm_grid_center[i]=(double *)malloc(sizeof(double)*num_dm);
      memset(norm_grid_center[i],0,sizeof(double)*num_dm);
    }

  for (i=0;i<file_Len;i++)
    {
      denseID=eventID_To_denseventID[i];
      if (denseID!=-1) //only dense events can enter
        {
          grid_ID[denseID]=dense_grid_reverse[all_grid_ID[i]];
			
          if (grid_ID[denseID]==-1)
            {
              fprintf(stderr,"Incorrect input file format or input parameters (no dense region found)!\n");
              exit(0);
            }			
        }
    }

	
  /* Find centroid (in the normalized data) for each dense grid */
  ID2Center(normalized_data,file_Len,num_dm,eventID_To_denseventID,num_dense_grids,grid_ID,norm_grid_center);	
 
  //printf("pass the grid ID2 center\n");

	
  graph_of_grid=(long **)malloc(sizeof(long*)*num_dense_grids);
  memset(graph_of_grid,0,sizeof(long*)*num_dense_grids);
  for (i=0;i<num_dense_grids;i++)
    {
      graph_of_grid[i]=(long *)malloc(sizeof(long)*num_dm);
      memset(graph_of_grid[i],0,sizeof(long)*num_dm);		
		

      for (j=0;j<num_dm;j++)
        graph_of_grid[i][j]=-1;
    }	
	
  grid_value=(long *)malloc(sizeof(long)*num_dm);
  memset(grid_value,0,sizeof(long)*num_dm);

  search_value=(long *)malloc(sizeof(long)*num_dm);
  memset(search_value,0,sizeof(long)*num_dm);

  
  for (i=0;i<num_dense_grids;i++)
    {
      ReventID=densegridID_To_gridevent[i];

      for (j=0;j<num_dm;j++)
        {
          grid_value[j]=position[ReventID][j];
          
        }
     

      /* For each dimension, find the neighbor, if any, that is equal in all other dimensions and 1 greater in
         the chosen dimension.  If the neighbor's centroid is not too far away, add it to this grid's neighbor
         list. */
      for (t=0;t<num_dm;t++)
        {
          for (p=0;p<num_dm;p++)
            search_value[p]=grid_value[p];

          if (grid_value[t]==num_bin-1)
            continue;

          search_value[t]=grid_value[t]+1; //we only consider the neighbor at the bigger side

          neighbor_ID=binary_search(num_dense_grids, num_dm, search_value, dense_grid_seq);
			
          if (neighbor_ID!=-1) 
            {
              real_dist=norm_grid_center[i][t]-norm_grid_center[neighbor_ID][t];
	
              if (real_dist<0)
                real_dist=-real_dist;
				
              if (real_dist<2*interval[t]) //by using 2*interval, this switch is void
                graph_of_grid[i][t]=neighbor_ID;			
            }
        }
      grid_clusterID[i]=i; //initialize grid_clusterID
    } 
	
	
  //graph constructed
  //DFS-based search begins

  /* Use a depth-first search to construct a list of connected subgraphs (= "clusters").
     Note our graph as we have constructed it is a DAG, so we can use that to our advantage
     in our search. */
  //  num_clust=dfs(graph_of_grid,num_dense_grids,num_dm,grid_clusterID);
  num_clust=find_connected(graph_of_grid, num_dense_grids, num_dm, grid_clusterID);

	
  //computes grid_ID and cluster_ID
  for (i=0;i<file_Len;i++)
    {
      denseID=eventID_To_denseventID[i];
      if (denseID!=-1) //only dense events can enter
	  {
        cluster_ID[denseID]=grid_clusterID[grid_ID[denseID]];
		//if (cluster_ID[denseID]==-1)
		//	printf("catch you!\n");
	  }
    }
	
  free(search_value);
  free(grid_value);

  for (i=0;i<num_dense_grids;i++)
    {
      free(graph_of_grid[i]);
      free(norm_grid_center[i]);
    }
  free(graph_of_grid);
  free(norm_grid_center);

  return num_clust;
}

/********************************************* Merge Clusters to Populations *******************************************/

long merge_clusters(long num_clust, long num_dm, double *interval, double **cluster_center, long *cluster_populationID, long max_num_pop)
{
  long num_population=0;
  long temp_num_population=0;

  long i=0;
  long j=0;
  long t=0;
  long merge=0;
  long smid=0;
  long bgid=0;
  double merge_dist=1.1;

  long *map_ID;

  double diff=0;  

  map_ID=(long*)malloc(sizeof(long)*num_clust);
  memset(map_ID,0,sizeof(long)*num_clust);

  for (i=0;i<num_clust;i++)
  {
    	cluster_populationID[i]=i;
		map_ID[i]=i;
  }
    

  while (((num_population>max_num_pop) && (merge_dist<5.0)) || ((num_population<=1) && (merge_dist>0.1)))
  {
  
	  if (num_population<=1)
	  	  merge_dist=merge_dist-0.1;

	  if (num_population>max_num_pop)
          merge_dist=merge_dist+0.1;
	  
	    
	 for (i=0;i<num_clust;i++)
		cluster_populationID[i]=i;

    for (i=0;i<num_clust-1;i++)
    {
      for (j=i+1;j<num_clust;j++)
        {
          merge=1;
			
          for (t=0;t<num_dm;t++)
            {
              diff=cluster_center[i][t]-cluster_center[j][t];
				
              if (diff<0)
                diff=-diff;
              if (diff>(merge_dist*interval[t]))
                merge=0;
            }

          if ((merge) && (cluster_populationID[i]!=cluster_populationID[j]))
            {
              if (cluster_populationID[i]<cluster_populationID[j])  //they could not be equal
                {
                  smid = cluster_populationID[i];
                  bgid = cluster_populationID[j];
                }
              else
                {
                  smid = cluster_populationID[j];
                  bgid = cluster_populationID[i];
                }
              for (t=0;t<num_clust;t++)
                {
                  if (cluster_populationID[t]==bgid)
                    cluster_populationID[t]=smid;
                }
            }
        }
    }

  

  for (i=0;i<num_clust;i++)
    map_ID[i]=-1;

  for (i=0;i<num_clust;i++)
    map_ID[cluster_populationID[i]]=1;

  num_population=0;
  for (i=0;i<num_clust;i++)
    {
      if (map_ID[i]==1)
        {
          map_ID[i]=num_population;
          num_population++;
        }
    }

  if ((temp_num_population>max_num_pop) && (num_population==1))
	  break;
  else
	  temp_num_population=num_population;

  if (num_clust<=1)
	break;
  } //end while

  for (i=0;i<num_clust;i++)
    cluster_populationID[i]=map_ID[cluster_populationID[i]];
	
  free(map_ID);

  return num_population; 
}

///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double kmeans(double **Matrix, long k, double kmean_term, long file_Len, long num_dm, long *shortest_id, double **center)
{
	
	long i=0;
	long j=0;
	long t=0;
	long random=0;
	long random1=0;
	long random2=0;
	long times=0;
	long dist_used=0; //0 is Euclidean and 1 is Pearson
	long random_init=0; //0: not use random seeds; 
	long real_Len=0;
	long skipped=0;
	
 	long *num;  //num[i]=t means the ith cluster has t points
	
	double vvv=1.0; // the biggest variation;
	double distance=0.0;
	double xv=0.0;
	double variation=0.0;

	double mean_dx=0;
	double mean_dy=0;
	double sum_var=0;
	double dx=0;
	double dy=0;
	double sd_x=0;
	double sd_y=0;	
	double diff=0;
	double distortion=0;
	double sd=0; //standard deviation
	double shortest_distance;

	double *temp_center;	
			
	double **sum;	
 
	temp_center = (double *)malloc(sizeof(double)*num_dm);
	memset(temp_center,0,sizeof(double)*num_dm);

	if (random_init)
	{
		for (i=0;i<k;i++)
		{	
			random1=rand()*rand();
			random2=abs((random1%5)+1);
			for (t=0;t<random2;t++)
				random2=random2*rand()+rand();
	
			random=abs(random2%file_Len);
			for (j=0;j<num_dm;j++)
				center[i][j]=Matrix[random][j];				
		}
	}


	num = (long *)malloc(sizeof(long)*k);
	memset(num,0,sizeof(long)*k);

	sum = (double **)malloc(sizeof(double*)*k);
	memset(sum,0,sizeof(double*)*k);
	for (i=0;i<k;i++)
	{
		sum[i] = (double *)malloc(sizeof(double)*num_dm);
		memset(sum[i],0,sizeof(double)*num_dm);
	}


	times=0;
	real_Len=0;

	while (((vvv>kmean_term) && (kmean_term<1)) || ((times<kmean_term) && (kmean_term>=1)))
	{
		for (i=0;i<k;i++)
		{
			num[i]=0;
			for (j=0;j<num_dm;j++)
				sum[i][j]=0.0;  
		}
		
		for (i=0;i<file_Len;i++)  //for each data point i, we compute the distance between Matrix[i] and center[j]
		{		
			skipped = 0;
			shortest_distance=MAX_VALUE;
			for (j=0;j<k;j++)  //for each center j
			{
				distance=0.0;
								
				if (dist_used==0)  //Euclidean distance
				{
					for (t=0;t<num_dm;t++) //for each dimension here num_dm is always 1 as we consider individual dimensions
					{
					
						diff=center[j][t]-Matrix[i][t];
					
						diff=diff*diff;  
						//this K-means is only used for dimension selection, and therefore for 1-dimensional data only, no need to use cube distance

						distance = distance+diff; //here we have a weight for each dimension
					}
				}
				else  //pearson correlation
				{
					mean_dx=0.0;
					mean_dy=0.0;
					sum_var=0.0;
					dx=0.0;
					dy=0.0;
					sd_x=0.0;
					sd_y=0.0;
					for (t=0;t<num_dm;t++)
					{
						mean_dx+=center[j][t];
						mean_dy+=Matrix[i][t];
					}
					mean_dx=mean_dx/(double)num_dm;
					mean_dy=mean_dy/(double)num_dm;
			
					for (t=0;t<num_dm;t++)
					{
						dx=center[j][t]-mean_dx;
						dy=Matrix[i][t]-mean_dy;
						sum_var+=dx*dy;
						sd_x+=dx*dx;
						sd_y+=dy*dy;
					}
					if (sqrt(sd_x*sd_y)==0)
						distance = 1.0;
					else
						distance = 1.0 - (sum_var/(sqrt(sd_x*sd_y))); // distance ranges from 0 to 2;
					//printf("distance=%f\n",distance);			
				}	//pearson correlation ends 				

				if ((distance<shortest_distance) && (skipped == 0))
				{
					shortest_distance=distance;						
					shortest_id[i]=j;  
				}			
			}//end for j
				real_Len++;
				num[shortest_id[i]]=num[shortest_id[i]]+1;
				for (t=0;t<num_dm;t++)
					sum[shortest_id[i]][t]=sum[shortest_id[i]][t]+Matrix[i][t];		
		}//end for i
	/* recompute the centers */
	//compute_mean(group);		
		vvv=0.0;
		for (j=0;j<k;j++)
		{
			memcpy(temp_center,center[j],sizeof(double)*num_dm);
			variation=0.0;
			if (num[j]!=0)
			{
				for (t=0;t<num_dm;t++)
				{
					center[j][t]=sum[j][t]/(double)num[j];
					xv=(temp_center[t]-center[j][t]);
					variation=variation+xv*xv;
				}
			}
			
			if (variation>vvv)
				vvv=variation;  //vvv is the biggest variation among the k clusters;			
		}
	//compute_variation;
		times++;
	} //end for while




	free(num);

	for (i=0;i<k;i++)
		free(sum[i]);
	free(sum);
	free(temp_center);


	return distortion;
		
}

//////////////////////////////////////////////////////
/*************************** Show *****************************/
void show(double **Matrix, long *cluster_id, long file_Len, long k, long num_dm, char *name_string)
{
	int situ1=0;
	int situ2=0;

	long i=0;
	long id=0;
	long j=0;
	long info_id=0;
	long nearest_id=0;
	long insert=0;
	long temp=0;
	long m=0;
	long n=0;
	long t=0;
	
	long *size_c;
	


	long **size_mybound_1;
	long **size_mybound_2;
	long **size_mybound_3;
	long **size_mybound_0;

	double interval=0.0;

	double *big;
	double *small;


	double **center;
	double **mybound;
	
	long **prof; //prof[i][j]=1 means population i is + at parameter j
	
	FILE *fpcnt_id; //proportion id
	//FILE *fcent_id; //center_id, i.e., centers of clusters within the original data
	FILE *fprof_id; //profile_id

	big=(double *)malloc(sizeof(double)*num_dm);
	memset(big,0,sizeof(double)*num_dm);

	small=(double *)malloc(sizeof(double)*num_dm);
	memset(small,0,sizeof(double)*num_dm);

	for (i=0;i<num_dm;i++)
	{
		big[i]=0.0;
		small[i]=(double)MAX_VALUE;
	}
	
	
	size_c=(long *)malloc(sizeof(long)*k);
	memset(size_c,0,sizeof(long)*k);

	center=(double**)malloc(sizeof(double*)*k);
	memset(center,0,sizeof(double*)*k);
	for (i=0;i<k;i++)
	{
		center[i]=(double*)malloc(sizeof(double)*num_dm);
		memset(center[i],0,sizeof(double)*num_dm);
	}

	mybound=(double**)malloc(sizeof(double*)*num_dm);
	memset(mybound,0,sizeof(double*)*num_dm);
	for (i=0;i<num_dm;i++) //there are 3 mybounds for 4 categories
	{
		mybound[i]=(double*)malloc(sizeof(double)*3);
		memset(mybound[i],0,sizeof(double)*3);
	}

	prof=(long **)malloc(sizeof(long*)*k);
	memset(prof,0,sizeof(long*)*k);
	for (i=0;i<k;i++)
	{
		prof[i]=(long *)malloc(sizeof(long)*num_dm);
		memset(prof[i],0,sizeof(long)*num_dm);
	}


	for (i=0;i<file_Len;i++)
	{
		id=cluster_id[i];
		for (j=0;j<num_dm;j++)
		{
			center[id][j]=center[id][j]+Matrix[i][j];
			if (big[j]<Matrix[i][j])
				big[j]=Matrix[i][j];
			if (small[j]>Matrix[i][j])
				small[j]=Matrix[i][j];
		}
		
		size_c[id]++;		
	}

	for (i=0;i<k;i++)
		for (j=0;j<num_dm;j++)
		{			
			if (size_c[i]!=0)
				center[i][j]=(center[i][j]/(double)(size_c[i]));
			else
				center[i][j]=0;	
		}

	for (j=0;j<num_dm;j++)
	{
		interval=((big[j]-small[j])/4.0);
		//printf("interval[%d] is %f\n",j,interval);
		for (i=0;i<3;i++)
			mybound[j][i]=small[j]+((double)(i+1)*interval);
	}
	

	size_mybound_0=(long **)malloc(sizeof(long*)*k);
	memset(size_mybound_0,0,sizeof(long*)*k);
	
	for (i=0;i<k;i++)
	{
		size_mybound_0[i]=(long*)malloc(sizeof(long)*num_dm);
		memset(size_mybound_0[i],0,sizeof(long)*num_dm);		
	}

	size_mybound_1=(long **)malloc(sizeof(long*)*k);
	memset(size_mybound_1,0,sizeof(long*)*k);
	
	for (i=0;i<k;i++)
	{
		size_mybound_1[i]=(long*)malloc(sizeof(long)*num_dm);
		memset(size_mybound_1[i],0,sizeof(long)*num_dm);		
	}

	size_mybound_2=(long **)malloc(sizeof(long*)*k);
	memset(size_mybound_2,0,sizeof(long*)*k);
	
	for (i=0;i<k;i++)
	{
		size_mybound_2[i]=(long*)malloc(sizeof(long)*num_dm);
		memset(size_mybound_2[i],0,sizeof(long)*num_dm);	
	}

	size_mybound_3=(long **)malloc(sizeof(long*)*k);
	memset(size_mybound_3,0,sizeof(long*)*k);

	for (i=0;i<k;i++)
	{
		size_mybound_3[i]=(long*)malloc(sizeof(long)*num_dm);
		memset(size_mybound_3[i],0,sizeof(long)*num_dm);
	}
	
	for (i=0;i<file_Len;i++)
		for (j=0;j<num_dm;j++)
			{
				if (Matrix[i][j]<mybound[j][0])// && ((Matrix[i][j]-small[j])>0)) //the smallest values excluded
					size_mybound_0[cluster_id[i]][j]++;
				else
				{
					if (Matrix[i][j]<mybound[j][1])
						size_mybound_1[cluster_id[i]][j]++;
					else
					{
						if (Matrix[i][j]<mybound[j][2])
							size_mybound_2[cluster_id[i]][j]++;
						else
							//if (Matrix[i][j]!=big[j]) //the biggest values excluded
								size_mybound_3[cluster_id[i]][j]++;
					}

				}
			}

	fprof_id=fopen("profile.txt","w");
	fprintf(fprof_id,"Population_ID\t");
	fprintf(fprof_id,"%s\n",name_string);
	
	for (i=0;i<k;i++)
	{
		fprintf(fprof_id,"%ld\t",i+1); //i changed to i+1 to start from 1 instead of 0: April 16, 2009
		for (j=0;j<num_dm;j++)
		{
			
			if (size_mybound_0[i][j]>size_mybound_1[i][j])
				situ1=0;
			else
				situ1=1;
			if (size_mybound_2[i][j]>size_mybound_3[i][j])
				situ2=2;
			else
				situ2=3;

			if ((situ1==0) && (situ2==2))
			{
				if (size_mybound_0[i][j]>size_mybound_2[i][j])
					prof[i][j]=0;
				else
					prof[i][j]=2;
			}
			if ((situ1==0) && (situ2==3))
			{
				if (size_mybound_0[i][j]>size_mybound_3[i][j])
					prof[i][j]=0;
				else
					prof[i][j]=3;
			}
			if ((situ1==1) && (situ2==2))
			{
				if (size_mybound_1[i][j]>size_mybound_2[i][j])
					prof[i][j]=1;
				else
					prof[i][j]=2;
			}
			if ((situ1==1) && (situ2==3))
			{
				if (size_mybound_1[i][j]>size_mybound_3[i][j])
					prof[i][j]=1;
				else
					prof[i][j]=3;
			}
			
			//begin to output profile
			if (j==num_dm-1)
			{
				if (prof[i][j]==0)
					fprintf(fprof_id,"1\n");
				if (prof[i][j]==1)
					fprintf(fprof_id,"2\n");
				if (prof[i][j]==2)
					fprintf(fprof_id,"3\n");
				if (prof[i][j]==3)
					fprintf(fprof_id,"4\n");
			}
			else
			{
				if (prof[i][j]==0)
					fprintf(fprof_id,"1\t");
				if (prof[i][j]==1)
					fprintf(fprof_id,"2\t");
				if (prof[i][j]==2)
					fprintf(fprof_id,"3\t");
				if (prof[i][j]==3)
					fprintf(fprof_id,"4\t");
			}
		}
	}
	fclose(fprof_id);

	///////////////////////////////////////////////////////////
	

	fpcnt_id=fopen("percentage.txt","w");
	fprintf(fpcnt_id,"Population_ID\tPercentage\n");

	for (t=0;t<k;t++)
	{
		fprintf(fpcnt_id,"%ld\t%.2f\n",t+1,(double)size_c[t]*100.0/(double)file_Len);	//t changed to t+1 to start from 1 instead of 0: April 16, 2009									
	}
	fclose(fpcnt_id);

	free(big);
	free(small);
	free(size_c);

	for (i=0;i<k;i++)
	{
		free(center[i]);
		free(prof[i]);
		free(size_mybound_0[i]);
		free(size_mybound_1[i]);
		free(size_mybound_2[i]);
		free(size_mybound_3[i]);
	}
	free(center);
	free(prof);
	free(size_mybound_0);
	free(size_mybound_1);
	free(size_mybound_2);
	free(size_mybound_3);

	for (i=0;i<num_dm;i++)
		free(mybound[i]);
	free(mybound);
	
}
///////////////////////////////////////////////////////////////////////////

double get_avg_dist(double *population_center, long smaller_pop_ID, long larger_pop_ID, long *population_ID, long num_population, long file_Len, 
				   long num_dm, double **normalized_data, long d1, long d2, long d3, long *size_c)
{
	long k=0;
	long temp=0;
	long i=0;
	long j=0;
	long t=0;
	long current_biggest=0;


	double total_dist=0.0;
	double distance=0.0;
	double dist=0.0;
	double biggest_distance=0.0;

	long *nearest_neighbors;
	double *nearest_distance;

	k=min(size_c[smaller_pop_ID],size_c[larger_pop_ID]);
	
	if (k>0)
	{
		k=(long)sqrt((double)k);
		if (num_dm<=3)
			k=(long)(2.7183*(double)k);  //e*k for low-dimensional space, added Nov. 4, 2010
		//printf("the k-nearest k is %d\n",k);
	}
	else
	{
		printf("error in k\n");
		exit(0);
	}

	nearest_neighbors=(long *)malloc(sizeof(long)*k);
	memset(nearest_neighbors,0,sizeof(long)*k);

	nearest_distance=(double *)malloc(sizeof(double)*k);
	memset(nearest_distance,0,sizeof(double)*k);

	temp=0;

	for (i=0;i<file_Len;i++)
	{
		if ((population_ID[i]==smaller_pop_ID) || (population_ID[i]==larger_pop_ID)) //the event belongs to the pair of populations
		{
			distance=0.0;
			for (t=0;t<num_dm;t++)
			{
				if ((t!=d1) && (t!=d2) && (t!=d3))
					continue;
				else
				{
					dist=population_center[t]-normalized_data[i][t];
					distance=distance+dist*dist;
				}
			}

			if (temp<k)
			{
				nearest_neighbors[temp]=i;
				nearest_distance[temp]=distance;
				temp++;
			}
			else
			{
			    biggest_distance=0.0;
				for (j=0;j<k;j++)
				{
					if (nearest_distance[j]>biggest_distance)
					{
						biggest_distance=nearest_distance[j];
						current_biggest=j;
					}
				}

				if (biggest_distance>distance)
				{
					nearest_neighbors[current_biggest]=i;
					nearest_distance[current_biggest]=distance;
				}
			}
		}
	}

    for (i=0;i<k;i++)
		total_dist=total_dist+nearest_distance[i];

	total_dist=total_dist/(double)k;

	free(nearest_distance);
	free(nearest_neighbors);
	
	return total_dist;
}


/******************************************************** Main Function **************************************************/
//for normalized data, there are five variables:
//cluster_ID
//population_center
//grid_clusterID
//grid_ID
//grid_Center

//the same five variables exist for the original data
//however, the three IDs (cluster_ID, grid_ID, grid_clusterID) don't change in either normalized or original data
//also, data + cluster_ID -> population_center
//data + grid_ID -> grid_Center

/* what is the final output */
//the final things we want are grid_Center in the original data and grid_clusterID //or population_center in the original data
//Sparse grids will not be considered when computing the two centroids (centroids of grids and centroids of clusters)

/*  what information should select_bin output */
//the size of all IDs are unknown to function main because we only consider the events in dense grids, and also the number of dense grids
//is unknown, therefore I must use a prescreening to output 
//how many bins I should use
//the number of dense grids
//total number of events in the dense grids

/* basic procedure of main function */ 
//1. read raw file and normalize the raw file
//2. select_bin
//3. allocate memory for eventID_To_denseventID, grid_clusterID, grid_ID, cluster_ID. 
//4. call select_dense and merge_grids with grid_clusterID, grid_ID, cluster_ID.
//5. release normalized data; allocate memory for grid_Center and population_center
//6. output grid_Center and population_center using ID2Center together with grid_clusterID //from IDs to centers

int main (int argc, char **argv)
{
  //inputs
  FILE *f_src; //source file pointer
 
  FILE *f_out; //coordinates
  FILE *f_cid; //population-ID of events
  FILE *f_ctr; //centroids of populations
  FILE *f_results; //coordinates file event and population column
  FILE *f_mfi; //added April 16, 2009 for mean fluorescence intensity
  FILE *f_parameters; //number of bins and density calculated by
                      //the algorithm. Used to update the database
  FILE *f_properties; //Properties file used by Image generation software

  //char tmpbuf[128];

  char para_name_string[LINE_LEN];

  long time_ID=-1;
  long num_bin=0; //the bins I will use on each dimension
  long num_pop=0;	
  long file_Len=0; //number of events		
  long num_dm=0;
  long num_clust=0;
  long num_dense_events=0;
  long num_dense_grids=0;
  long num_nonempty=0;
  long num_population=0;
  long num_real_pop=0;
  long keep_merge=0;
  long num_checked_range=0;

  //below are read from configuration file
  long i=0;
  long j=0;
  long p=0;
  long t=0;
  long q=0;
  long temp_i=0;
  long temp_j=0;
  long first_i=0;
  long first_j=0;
  long d1=0;
  long d2=0;
  long d3=0;
  long d_d=0;
  long max_num_pop=0; //upperbound of the number of populations that is equal to 2^d
  long index_id=0;
  long num_computed_population=0;
  long tot_size=0;

  int den_t_event=0;

  double distance=0.0;
  double dist=0.0;
  double current_smallest_dist=0;
  double max_d_dist=0.0;
  double Ei=0.0;
  double Ej=0.0;
  double E1=0.0;
  double E2=0.0;
  double E3=0.0;
  double Ep=0.0;
  double temp=0.0;
  double tmp=0.0;

  long *grid_clusterID; //(dense)gridID to clusterID
  long *grid_ID; //(dense)eventID to gridID
  long *cluster_ID; //(dense)eventID to clusterID
  long *eventID_To_denseventID; //eventID_To_denseventID[i]=k means event i is in a dense grid and its ID within dense events is k
  long *all_grid_ID; //gridID for all events
  long *densegridID_To_gridevent;
  long *sorted_seq;
  long *dense_grid_reverse;
  long *population_ID; //denseeventID to populationID
  long *cluster_populationID; //clusterID to populationID
  long *grid_populationID; //gridID to populationID
  long *all_population_ID; //populationID of event
  long *size_c;
  long *size_p_2;
  long *partit;
  long *size_p;
  long *temp_size_j;
  long *all_computed_pop_ID;


  long **position;
  long **dense_grid_seq;
  long **temp_population_ID;

  double *interval;
  double *center_1;
  double *center_2;
  double *center_3;
  double *aver;
  double *std;
  double *center_p_1;
  double *center_p_2;
  double *center_p_3;
	
  double **population_center; //population centroids in the raw/original data
  double **cluster_center; //cluster centroids in the raw/original data
  double **new_population_center;
  double **temp_population_center;
  double **min_pop;
  double **max_pop;

  double **input_data;
  double **normalized_data;

  double **real_pop_center;
  double **distance_pop;

  double ***ind_pop;
  double ***ind_pop_3;
		
  int min = 999999;
  int max = 0;

  printf( "Starting time:\t\t\t\t");
  fflush(stdout);
  system("/bin/date");

  /*
   * Windows Version
  _strtime( tmpbuf );
  printf( "Starting time:\t\t\t\t%s\n", tmpbuf );
  _strdate( tmpbuf );
  printf( "Starting date:\t\t\t\t%s\n", tmpbuf );
  */

  /////////////////////////////////////////////////////////////

  if ((argc!=2) && (argc!=3) && (argc!=4) && (argc!=5) && (argc!=6))
  {
      fprintf(stderr,"Incorrect number of input parameters!\n"); 
	  fprintf(stderr,"usage:\n");
	  fprintf(stderr,"basic mode: flock data_file\n");
	  fprintf(stderr,"advanced mode 0 (specify maximum # of pops): flock data_file max_num_pop\n");
	  fprintf(stderr,"advanced mode 1 (without # of pops): flock data_file num_bin density_index\n");
	  fprintf(stderr,"advanced mode 2 (specify # of pops): flock data_file num_bin density_index number_of_pop\n");
	  fprintf(stderr,"advanced mode 3 (specify both # of pops): flock data_file num_bin density_index number_of_pop max_num_pop\n");
      exit(0);
  }	

  f_src=fopen(argv[1],"r");

  if (argc==3)
  {
	max_num_pop=atoi(argv[2]);
	printf("num of maximum pop is %ld\n",max_num_pop);
  }
  
  if (argc==4)
  {
	num_bin=atoi(argv[2]);
	printf("num_bin is %ld\n",num_bin);

	den_t_event=atoi(argv[3]);
	printf("density_index is %d\n",den_t_event);
  }

  if (argc==5)
  {
	num_bin=atoi(argv[2]);
	printf("num_bin is %ld\n",num_bin);

	den_t_event=atoi(argv[3]);
	printf("density_index is %d\n",den_t_event);

	num_pop=atoi(argv[4]);
	printf("num of pop is %ld\n",num_pop);
  }

  if (argc==6)
  {
	num_bin=atoi(argv[2]);
	printf("num_bin is %ld\n",num_bin);

	den_t_event=atoi(argv[3]);
	printf("density_index is %d\n",den_t_event);

	num_pop=atoi(argv[4]);
	printf("num of pop is %ld\n",num_pop);

	max_num_pop=atoi(argv[5]);
	printf("num of pop is %ld\n",max_num_pop);
  }


  if (num_pop==1)
  {
	  printf("it is not allowed to specify only one population\n");
	  exit(0);
  }

    

  getfileinfo(f_src, &file_Len, &num_dm, para_name_string, &time_ID); //get the filelength, number of dimensions, and num/name of parameters

  
  if (max_num_pop==0)
  {
	  max_num_pop=(long)pow(2,num_dm);
		if (max_num_pop>MAX_POP_NUM)
			max_num_pop=MAX_POP_NUM;
  }

  /************************************************* Read the data *****************************************************/
	
  rewind(f_src); //reset the file pointer	

  input_data = (double **)malloc(sizeof(double*)*file_Len);
  memset(input_data,0,sizeof(double*)*file_Len);
  for (i=0;i<file_Len;i++)
  {
     input_data[i]=(double *)malloc(sizeof(double)*num_dm);
     memset(input_data[i],0,sizeof(double)*num_dm);
  }
	
  readsource(f_src, file_Len, num_dm, input_data, time_ID); //read the data;
	
  fclose(f_src);

  normalized_data=(double **)malloc(sizeof(double*)*file_Len);
  memset(normalized_data,0,sizeof(double*)*file_Len);
  for (i=0;i<file_Len;i++)
    {
      normalized_data[i]=(double *)malloc(sizeof(double)*num_dm);
      memset(normalized_data[i],0,sizeof(double)*num_dm);
    }
	
  tran(input_data, file_Len, num_dm, NORM_METHOD, normalized_data);
	

  position=(long **)malloc(sizeof(long*)*file_Len);
  memset(position,0,sizeof(long*)*file_Len);
  for (i=0;i<file_Len;i++)
    {
      position[i]=(long*)malloc(sizeof(long)*num_dm);
      memset(position[i],0,sizeof(long)*num_dm);
    }

  all_grid_ID=(long *)malloc(sizeof(long)*file_Len);
  memset(all_grid_ID,0,sizeof(long)*file_Len);

  sorted_seq=(long*)malloc(sizeof(long)*file_Len);
  memset(sorted_seq,0,sizeof(long)*file_Len);
	
  interval=(double*)malloc(sizeof(double)*num_dm);
  memset(interval,0,sizeof(double)*num_dm);

  /************************************************* select_bin *************************************************/
	
  if (num_bin>=1)  //num_bin has been selected by user
  {
  	select_bin(normalized_data, file_Len, num_dm, MIN_GRID, MAX_GRID, position, sorted_seq, all_grid_ID, &num_nonempty, interval,num_bin);
	printf("user set bin number is %ld\n",num_bin); 
  }
  else  //num_bin has not been selected by user
  {
	num_bin=select_bin(normalized_data, file_Len, num_dm, MIN_GRID, MAX_GRID, position, sorted_seq, all_grid_ID, &num_nonempty, interval,num_bin);
	printf("selected bin number is %ld\n",num_bin);    
  }
  printf("number of non-empty grids is %ld\n",num_nonempty);
		
 

  /* Although we return sorted_seq from select_bin(), we don't use it for anything, except possibly diagnostics */
  free(sorted_seq);
	
	
  dense_grid_reverse=(long*)malloc(sizeof(long)*num_nonempty);
  memset(dense_grid_reverse,0,sizeof(long)*num_nonempty);	

  /************************************************* select_dense **********************************************/

  if (den_t_event>=1) //den_t_event must be larger or equal to 2 if the user wants to set it
	den_t_event=select_dense(file_Len, all_grid_ID, num_nonempty, &num_dense_grids, &num_dense_events, dense_grid_reverse, den_t_event);
  else
  {
	den_t_event=select_dense(file_Len, all_grid_ID, num_nonempty, &num_dense_grids, &num_dense_events, dense_grid_reverse, den_t_event);
	printf("automated selected density threshold is %d\n",den_t_event);
  }		

  printf("Number of dense grids is %ld\n",num_dense_grids);

  densegridID_To_gridevent = (long *)malloc(sizeof(long)*num_dense_grids);
  memset(densegridID_To_gridevent,0,sizeof(long)*num_dense_grids);

  for (i=0;i<num_dense_grids;i++)
    densegridID_To_gridevent[i]=-1; //initialize all densegridID_To_gridevent values to -1
	

  eventID_To_denseventID=(long *)malloc(sizeof(long)*file_Len);
  memset(eventID_To_denseventID,0,sizeof(long)*file_Len);     //eventID_To_denseventID[i]=k means event i is in a dense grid and its ID within dense events is k


  grid_To_event(file_Len, dense_grid_reverse, all_grid_ID, eventID_To_denseventID, densegridID_To_gridevent);

	
  dense_grid_seq=(long **)malloc(sizeof(long*)*num_dense_grids);
  memset(dense_grid_seq,0,sizeof(long*)*num_dense_grids);
  for (i=0;i<num_dense_grids;i++)
    {
      dense_grid_seq[i]=(long *)malloc(sizeof(long)*num_dm);
      memset(dense_grid_seq[i],0,sizeof(long)*num_dm);
    }


  /* Look up the binned data values for each dense grid */
  generate_grid_seq(num_dm, num_dense_grids, densegridID_To_gridevent, position, dense_grid_seq); 	
	
	
  /************************************************* allocate memory *********************************************/
	
  grid_clusterID=(long *)malloc(sizeof(long)*num_dense_grids);
  memset(grid_clusterID,0,sizeof(long)*num_dense_grids);

  grid_ID=(long *)malloc(sizeof(long)*num_dense_events);
  memset(grid_ID,0,sizeof(long)*num_dense_events);

  cluster_ID=(long *)malloc(sizeof(long)*num_dense_events);
  memset(cluster_ID,0,sizeof(long)*num_dense_events);


  /*********************************************** merge_grids ***********************************************/
  //long merge_grids(long file_Len, long num_dm, long num_bin, long **position, long num_dense_grids, long *dense_grid_rank, long *dense_grid_reverse,
  //			 long **dense_grid_seq, long *eventID_To_denseventID, long *densegridID_To_gridevent, long *all_grid_ID,
  //			 long *cluster_ID, long *grid_ID, long *grid_clusterID)
	
  num_clust = merge_grids(normalized_data, interval, file_Len, num_dm, num_bin, position, num_dense_grids, dense_grid_reverse, dense_grid_seq, eventID_To_denseventID, densegridID_To_gridevent, all_grid_ID, cluster_ID, grid_ID, grid_clusterID);
	
  printf("computed number of groups is %ld\n",num_clust);	

	
  /************************************** release unnecessary memory and allocate memory and compute centers **********************************/
	
	
  for (i=0;i<file_Len;i++)
    free(position[i]);
  free(position);

  for (i=0;i<num_dense_grids;i++)
    free(dense_grid_seq[i]);
  free(dense_grid_seq);

  free(dense_grid_reverse);
	
  free(densegridID_To_gridevent);
  free(all_grid_ID);
	
  // cluster_center ////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  cluster_center=(double **)malloc(sizeof(double*)*num_clust);
  memset(cluster_center,0,sizeof(double*)*num_clust);
  for (i=0;i<num_clust;i++)
  {
     cluster_center[i]=(double*)malloc(sizeof(double)*num_dm);
     memset(cluster_center[i],0,sizeof(double)*num_dm);
  }
	
  ID2Center(normalized_data,file_Len,num_dm,eventID_To_denseventID,num_clust,cluster_ID,cluster_center); //produce the centers with normalized data
	
  //printf("pass the first ID2center\n");

  /*** population_ID and grid_populationID **/
	
  cluster_populationID=(long*)malloc(sizeof(long)*num_clust);
  memset(cluster_populationID,0,sizeof(long)*num_clust);

  grid_populationID=(long*)malloc(sizeof(long)*num_dense_grids);
  memset(grid_populationID,0,sizeof(long)*num_dense_grids);

  population_ID=(long*)malloc(sizeof(long)*num_dense_events);
  memset(population_ID,0,sizeof(long)*num_dense_events);

  num_population = merge_clusters(num_clust, num_dm, interval, cluster_center, cluster_populationID, max_num_pop);

  

  for (i=0;i<num_clust;i++)
    free(cluster_center[i]);
  free(cluster_center);

  free(interval);
	
  for (i=0;i<num_dense_grids;i++)
    {
      grid_populationID[i]=cluster_populationID[grid_clusterID[i]];
    }

  for (i=0;i<num_dense_events;i++)
    {
      population_ID[i]=cluster_populationID[cluster_ID[i]];
    }

  printf("Number of merged groups before second-run partitioning is %ld\n",num_population);

  

  // population_center /////////////////////////////////////////////////////////////////////////////////////////////////////

 
  population_center=(double **)malloc(sizeof(double*)*num_population);
  memset(population_center,0,sizeof(double*)*num_population);
  for (i=0;i<num_population;i++)
  {
     population_center[i]=(double*)malloc(sizeof(double)*num_dm);
     memset(population_center[i],0,sizeof(double)*num_dm);
  }

  
	
  ID2Center(normalized_data,file_Len,num_dm,eventID_To_denseventID,num_population,population_ID,population_center); //produce population centers with normalized data
	

  // show ////////////////////////////////////////////////////////////////////////////////
  all_population_ID=(long*)malloc(sizeof(long)*file_Len);
  memset(all_population_ID,0,sizeof(long)*file_Len);

  kmeans(normalized_data, num_population, KMEANS_TERM, file_Len, num_dm, all_population_ID, population_center);

  for (i=0;i<num_population;i++)
	free(population_center[i]);

  free(population_center);
  
  ////// there needs to be another step to further partition the data
  
	
	  center_p_1=(double *)malloc(sizeof(double)*3);
	  memset(center_p_1,0,sizeof(double)*3);

	  center_p_2=(double *)malloc(sizeof(double)*3);
	  memset(center_p_2,0,sizeof(double)*3);

	  center_p_3=(double *)malloc(sizeof(double)*3);
	  memset(center_p_3,0,sizeof(double)*3);

	  size_p=(long *)malloc(sizeof(long)*num_population);
	  memset(size_p,0,sizeof(long)*num_population);
	
	  temp_size_j=(long *)malloc(sizeof(long)*num_population);
	  memset(temp_size_j,0,sizeof(long)*num_population);

	  all_computed_pop_ID=(long *)malloc(sizeof(long)*file_Len);
	  memset(all_computed_pop_ID,0,sizeof(long)*file_Len);

	  for (i=0;i<file_Len;i++)
	  {
		size_p[all_population_ID[i]]++;
		all_computed_pop_ID[i]=all_population_ID[i];
	  }

	  ind_pop=(double***)malloc(sizeof(double**)*num_population);
	  memset(ind_pop,0,sizeof(double**)*num_population);

	  ind_pop_3=(double***)malloc(sizeof(double**)*num_population);
	  memset(ind_pop_3,0,sizeof(double**)*num_population);

	  for (i=0;i<num_population;i++)
	  {
		  ind_pop[i]=(double**)malloc(sizeof(double*)*size_p[i]);
		  memset(ind_pop[i],0,sizeof(double*)*size_p[i]);

		  ind_pop_3[i]=(double**)malloc(sizeof(double*)*size_p[i]);
		  memset(ind_pop_3[i],0,sizeof(double*)*size_p[i]);

		  for (j=0;j<size_p[i];j++)
		  {
			  ind_pop[i][j]=(double*)malloc(sizeof(double)*num_dm);
			  memset(ind_pop[i][j],0,sizeof(double)*num_dm);

			  ind_pop_3[i][j]=(double*)malloc(sizeof(double)*3);
			  memset(ind_pop_3[i][j],0,sizeof(double)*3);
		  }
		  temp_size_j[i]=0;
	  }

	  for (i=0;i<file_Len;i++) //to generate ind_pop[i]
	  {
		  index_id=all_population_ID[i];
		  
		  j=temp_size_j[index_id];

		  for (t=0;t<num_dm;t++)
		  {
			ind_pop[index_id][j][t]=input_data[i][t];
		  }
		  temp_size_j[index_id]++;		
	  }

	  aver=(double*)malloc(sizeof(double)*num_dm);
	  memset(aver,0,sizeof(double)*num_dm);

	  std=(double*)malloc(sizeof(double)*num_dm);
	  memset(std,0,sizeof(double)*num_dm);

	  temp_population_center=(double**)malloc(sizeof(double*)*2);
	  memset(temp_population_center,0,sizeof(double*)*2);
	  for (i=0;i<2;i++)
	  {
		temp_population_center[i]=(double*)malloc(sizeof(double)*3);
		memset(temp_population_center[i],0,sizeof(double)*3);
	  }

	  size_p_2=(long*)malloc(sizeof(long)*2);
	  memset(size_p_2,0,sizeof(long)*2);

	  partit=(long*)malloc(sizeof(long)*num_population);
	  memset(partit,0,sizeof(long)*num_population);

	  num_computed_population=num_population;

	  temp_population_ID=(long**)malloc(sizeof(long*)*num_population);
	  memset(temp_population_ID,0,sizeof(long*)*num_population);

	  for (i=0;i<num_population;i++)
	  {
		temp_population_ID[i]=(long*)malloc(sizeof(long)*size_p[i]);
		memset(temp_population_ID[i],0,sizeof(long)*size_p[i]);
	  }
	
	  //printf("num_population=%d\n",num_population);
	  
	  for (i=0;i<num_population;i++) 
	  {
		  partit[i]=0;
		  tran(ind_pop[i], size_p[i], num_dm, 1, ind_pop[i]);//0-1 normalize ind_pop[i]
          
		  //find the 3 dimensions with the largest std for ind_pop[i]
		  d1=-1;
		  d2=-1;
		  d3=-1;

		  for (t=0;t<num_dm;t++)
		  {
			aver[t]=0;
			for (j=0;j<size_p[i];j++)
			{
				aver[t]=aver[t]+ind_pop[i][j][t];
			}
			aver[t]=aver[t]/(double)size_p[i];

			std[t]=0;
			for (j=0;j<size_p[i];j++)
			{
				std[t]=std[t]+((ind_pop[i][j][t]-aver[t])*(ind_pop[i][j][t]-aver[t]));
			}
			std[t]=sqrt(std[t]/(double)size_p[i]);
		  }

		  

		  for (j=0;j<3;j++)
		  {
			max_d_dist=0;
			for (t=0;t<num_dm;t++)
			{
				if ((t!=d1) && (t!=d2))
				{
					dist=std[t];
				}
				else
					dist=-1;

				if (dist>max_d_dist)
				{
					max_d_dist=dist;
					d_d=t;						
				}
			}

			if (j==0)
				d1=d_d;
			if (j==1)
				d2=d_d;
			if (j==2)
				d3=d_d;
		  }
		 
		 
		  for (t=0;t<size_p[i];t++)
		  {
				ind_pop_3[i][t][0]=ind_pop[i][t][d1];
				ind_pop_3[i][t][1]=ind_pop[i][t][d2];
				ind_pop_3[i][t][2]=ind_pop[i][t][d3];
		  }
		

		 temp_population_center[0][0]=(aver[d1])/2.0;
		 
		 temp_population_center[0][1]=aver[d2];
	     temp_population_center[0][2]=aver[d3];
		 
		 temp_population_center[1][0]=(aver[d1]+1.0)/2.0;

		 temp_population_center[1][1]=aver[d2];
		 temp_population_center[1][2]=aver[d3];
		 
		 

		 //run K-means with K=2
		 kmeans(ind_pop_3[i], 2, KMEANS_TERM, size_p[i], 3, temp_population_ID[i], temp_population_center);

		 for (j=0;j<2;j++)
			 size_p_2[j]=0;

		 for (j=0;j<size_p[i];j++)
			 size_p_2[temp_population_ID[i][j]]++;

		 for (j=0;j<2;j++)
			 if (size_p_2[j]==0)
				 printf("size_p_2=0 at i=%ld\n",i);


		 //check whether the 2 parts should be merged
		

		 Ei=get_avg_dist(temp_population_center[0], 0,1, temp_population_ID[i], 2,size_p[i], 3, ind_pop_3[i], 0,1,2,size_p_2);
		 Ej=get_avg_dist(temp_population_center[1], 0,1, temp_population_ID[i], 2,size_p[i], 3, ind_pop_3[i], 0,1,2,size_p_2);

		 for (t=0;t<3;t++)
		 {
			center_p_1[t]=temp_population_center[0][t]*0.25+temp_population_center[1][t]*0.75;
			center_p_2[t]=temp_population_center[0][t]*0.5+temp_population_center[1][t]*0.5;
			center_p_3[t]=temp_population_center[0][t]*0.75+temp_population_center[1][t]*0.25;
		 }
		
		 E1=get_avg_dist(center_p_1, 0,1, temp_population_ID[i], 2,size_p[i], 3, ind_pop_3[i], 0,1,2,size_p_2);
		
		 E2=get_avg_dist(center_p_2, 0,1, temp_population_ID[i], 2,size_p[i], 3, ind_pop_3[i], 0,1,2,size_p_2);

		 E3=get_avg_dist(center_p_3, 0,1, temp_population_ID[i], 2,size_p[i], 3, ind_pop_3[i], 0,1,2,size_p_2);

		 //printf("i=%d;E1=%f;E2=%f;E3=%f\n",i,E1,E2,E3);

		 if (E1<E2)
			Ep=E2;
		 else
			Ep=E1;
			
		 Ep=max(Ep,E3); //Ep is the most sparse area

		 

		 if ((Ep>Ei) && (Ep>Ej)) //the two parts should be partitioned
		 {
			 partit[i]=num_computed_population;
			 num_computed_population++;			 
		 }
		
	  } //end for (i=0;i<num_population;i++)

	  
	  for (i=0;i<num_population;i++)
		  temp_size_j[i]=0;

	  for (i=0;i<file_Len;i++)
	  {
		index_id=all_population_ID[i];
		if (partit[index_id]>0)
		{
			j=temp_size_j[index_id];
			if (temp_population_ID[index_id][j]==1)
				all_computed_pop_ID[i]=partit[index_id];
					
			temp_size_j[index_id]++;
		}
	  }

	  printf("Number of groups after partitioning is %ld\n", num_computed_population);

	
      free(size_p_2);
	  
	  free(aver);
	  free(std);
	  free(partit);

	  free(center_p_1);
	  free(center_p_2);
	  free(center_p_3);

	 

	  for (i=0;i<num_population;i++)
		  free(temp_population_ID[i]);
	  free(temp_population_ID);

	  for (i=0;i<num_population;i++)
	  {
		for (j=0;j<size_p[i];j++)
		{
			  free(ind_pop[i][j]);
			  free(ind_pop_3[i][j]);
		}
	  }

	 

	  for (i=0;i<num_population;i++)
	  {
		  free(ind_pop[i]);
		  free(ind_pop_3[i]);
	  }
	  free(ind_pop);
	  free(ind_pop_3);

	  
	 

	  for (i=0;i<2;i++)
		  free(temp_population_center[i]);

	  free(temp_population_center);

	  free(temp_size_j);
	  free(size_p);

	 

	  //update the IDs, Centers, and # of populations
	  for (i=0;i<file_Len;i++)
	  {
		  all_population_ID[i]=all_computed_pop_ID[i];
		  if (all_population_ID[i]>=num_computed_population)
			  printf("all_population_ID[%ld]=%ld\n",i,all_population_ID[i]);
	  }

	  num_population=num_computed_population;

	  free(all_computed_pop_ID);
    //end partitioning
  // since the num_population has changed, population_center needs to be redefined as below

  population_center=(double **)malloc(sizeof(double*)*num_population);
  memset(population_center,0,sizeof(double*)*num_population);
  for (i=0;i<num_population;i++)
  {
     population_center[i]=(double*)malloc(sizeof(double)*num_dm);
     memset(population_center[i],0,sizeof(double)*num_dm);
  }

  	  
  ID2Center_all(normalized_data,file_Len,num_dm,num_population,all_population_ID,population_center);
  


  ////// end of further partitioning

  ////////////////////////////////////////////////////////////
  //  to further merge populations to avoid overpartitioning
  //  Added June 20, 2010
  ////////////////////////////////////////////////////////////
  

  
  num_real_pop=num_population;
  keep_merge=1;

  if (num_pop>num_real_pop)
  {
		fprintf(stderr,"number of populations specified too large\n");
		exit(0);
  }

  center_1=(double *)malloc(sizeof(double)*num_dm);
  memset(center_1,0,sizeof(double)*num_dm);

  center_2=(double *)malloc(sizeof(double)*num_dm);
  memset(center_2,0,sizeof(double)*num_dm);

  center_3=(double *)malloc(sizeof(double)*num_dm);
  memset(center_3,0,sizeof(double)*num_dm);


    
  while (((num_real_pop>num_pop) && (num_pop!=0)) || ((keep_merge==1) && (num_pop==0) && (num_real_pop>2)))
  {
	  
	keep_merge=0;  //so, if no entering merge function, while will exit when num_pop==0

	real_pop_center=(double **)malloc(sizeof(double*)*num_real_pop);
    memset(real_pop_center,0,sizeof(double*)*num_real_pop);
    for (i=0;i<num_real_pop;i++)
    {
      real_pop_center[i]=(double*)malloc(sizeof(double)*num_dm);
      memset(real_pop_center[i],0,sizeof(double)*num_dm);
    }

	min_pop=(double **)malloc(sizeof(double*)*num_real_pop);
    memset(min_pop,0,sizeof(double*)*num_real_pop);
    for (i=0;i<num_real_pop;i++)
    {
      min_pop[i]=(double*)malloc(sizeof(double)*num_dm);
      memset(min_pop[i],0,sizeof(double)*num_dm);
    }

	max_pop=(double **)malloc(sizeof(double*)*num_real_pop);
    memset(max_pop,0,sizeof(double*)*num_real_pop);
    for (i=0;i<num_real_pop;i++)
    {
      max_pop[i]=(double*)malloc(sizeof(double)*num_dm);
      memset(max_pop[i],0,sizeof(double)*num_dm);
    }

	if (num_real_pop==num_population)
	{
		for (i=0;i<num_real_pop;i++)
			for (j=0;j<num_dm;j++)
				real_pop_center[i][j]=population_center[i][j];	
	}
	else
	{
		ID2Center_all(normalized_data,file_Len,num_dm,num_real_pop,all_population_ID,real_pop_center);
	}

	for (i=0;i<num_real_pop;i++)
	{
		for (j=0;j<num_dm;j++)
		{
			min_pop[i][j]=MAX_VALUE;
			
			max_pop[i][j]=0;
		}			
	}
	
	for (i=0;i<file_Len;i++)
	{
		index_id=all_population_ID[i];
		for (j=0;j<num_dm;j++)
		{
			if (normalized_data[i][j]<min_pop[index_id][j])
				min_pop[index_id][j]=normalized_data[i][j];
			
			if (normalized_data[i][j]>max_pop[index_id][j])
				max_pop[index_id][j]=normalized_data[i][j];
		}			
	}

/*	for (i=0;i<num_real_pop;i++)
	{
		for (j=0;j<num_dm;j++)
		{
			printf("min_pop is %f\t",min_pop[i][j]);
			//printf("max_pop is %f\n",max_pop[i][j]);
		}
		printf("\n");
	}
*/


	distance_pop=(double **)malloc(sizeof(double*)*num_real_pop);
    memset(distance_pop,0,sizeof(double*)*num_real_pop);
    for (i=0;i<num_real_pop;i++)
    {
      distance_pop[i]=(double*)malloc(sizeof(double)*num_real_pop);
      memset(distance_pop[i],0,sizeof(double)*num_real_pop);
    }

	for (i=0;i<num_real_pop-1;i++)
	{
		for (j=i+1;j<num_real_pop;j++)
		{
		  distance=0;
		
		  for (t=0;t<num_dm;t++)
		  {
			dist=real_pop_center[i][t]-real_pop_center[j][t];
			distance=distance+dist*dist;
		  }
		  distance_pop[i][j]=distance;
		  distance_pop[j][i]=distance;	      
		}
	}
 
	if (num_real_pop>2)
		num_checked_range=(long)((double)num_real_pop*(num_real_pop-1)/2.0); //to find the mergeable pair among the num_checked_range pairs of smallest distances
	else
		num_checked_range=1;

	size_c=(long *)malloc(sizeof(long)*num_real_pop);
    memset(size_c,0,sizeof(long)*num_real_pop);
	
	for (i=0;i<file_Len;i++)
	  size_c[all_population_ID[i]]++;

	for (p=0;p<num_checked_range;p++)  
	{
		current_smallest_dist=MAX_VALUE;

		for (i=0;i<num_real_pop-1;i++)
		{
			for (j=i+1;j<num_real_pop;j++)
			{
				if (distance_pop[i][j]<current_smallest_dist)
				{
					current_smallest_dist=distance_pop[i][j];
				
					temp_i=i;
					temp_j=j;
				}
			}
		}

		if (p==0)
		{
			first_i=temp_i;
			first_j=temp_j;
		}

		distance_pop[temp_i][temp_j]=MAX_VALUE+1; //make sure this pair won't be selected next time
		
		//normalize and calculate std for temp_i and temp_j
		
		//pick the largest 3 dimensions on standard deviation
		d1=-1;
		d2=-1;
		d3=-1;

		for (i=0;i<3;i++)
		{
			max_d_dist=0;
			for (j=0;j<num_dm;j++)
			{
				if ((j!=d1) && (j!=d2))
				{
					temp=min(min_pop[temp_i][j],min_pop[temp_j][j]);
					tmp=max(max_pop[temp_i][j],max_pop[temp_j][j]);
					if (tmp<=temp)
					{
						//printf("min_pop[%d][%d]=%f \t min_pop[%d][%d]=%f \t max_pop[%d][%d]=%f \t max_pop[%d][%d]=%f\n",temp_i,j,min_pop[temp_i][j],temp_j, j, min_pop[temp_j][j],temp_i,j,max_pop[temp_i][j],temp_j,j,max_pop[temp_j][j]);
						dist=0;
					}
					else
						dist=(real_pop_center[temp_i][j]-real_pop_center[temp_j][j])/(tmp-temp);
					if (dist<0)
						dist=-dist;
				}
				else
					dist=-1;

				if (dist>max_d_dist)
				{
					max_d_dist=dist;
					d_d=j;						
				}
			}

			if (i==0)
				d1=d_d;
			if (i==1)
				d2=d_d;
			if (i==2)
				d3=d_d;
		}

		//printf("d1=%d;d2=%d;d3=%d\n",d1,d2,d3);

		Ei=get_avg_dist(real_pop_center[temp_i], temp_i,temp_j, all_population_ID, num_real_pop,file_Len, num_dm, normalized_data, d1,d2,d3,size_c);
		Ej=get_avg_dist(real_pop_center[temp_j], temp_i,temp_j, all_population_ID, num_real_pop,file_Len, num_dm, normalized_data, d1,d2,d3,size_c);

		for (t=0;t<num_dm;t++)
		{
			center_1[t]=real_pop_center[temp_i][t]*0.25+real_pop_center[temp_j][t]*0.75;
			center_2[t]=real_pop_center[temp_i][t]*0.5+real_pop_center[temp_j][t]*0.5;
			center_3[t]=real_pop_center[temp_i][t]*0.75+real_pop_center[temp_j][t]*0.25;
		}
		
		E1=get_avg_dist(center_1, temp_i,temp_j, all_population_ID, num_real_pop,file_Len, num_dm, normalized_data, d1,d2,d3,size_c);
		
		E2=get_avg_dist(center_2, temp_i,temp_j, all_population_ID, num_real_pop,file_Len, num_dm, normalized_data, d1,d2,d3,size_c);

		E3=get_avg_dist(center_3, temp_i,temp_j, all_population_ID, num_real_pop,file_Len, num_dm, normalized_data, d1,d2,d3,size_c);

		if (E1<E2)
			Ep=E2;
		else
			Ep=E1;
			
		Ep=max(Ep,E3); //Ep is the most sparse area

		Ep=Ep*E_T;
		//printf("Ep=%f;Ei=%f;Ej=%f\n",Ep,Ei,Ej);

		if ((Ep<=Ei) || (Ep<=Ej))//if the most sparse area between i and j are still denser than one of them, temp_i and temp_j should be merged
		{
			keep_merge=1;
			break;		
		}//end if (Ep)
	} //end for p

	//printf("keep_merge=%d\n",keep_merge);

	for (i=0;i<num_real_pop;i++)
	{
		free(real_pop_center[i]);
		free(distance_pop[i]);
		free(min_pop[i]);
		free(max_pop[i]);
	}
	free(real_pop_center);
	free(min_pop);
	free(max_pop);
	free(distance_pop);
	free(size_c);

	//printf("temp_i=%d;temp_j=%d\n",temp_i,temp_j);

	if (keep_merge)  //found one within p loop
	{
		for (i=0;i<file_Len;i++)
		{
			if (all_population_ID[i]>temp_j)
			{
				all_population_ID[i]=all_population_ID[i]-1;
			}
			else 
			{
				if (all_population_ID[i]==temp_j)
				{
					all_population_ID[i]=temp_i;
				}
			}
		}

		num_real_pop--;
	}
	else
	{
		if ((num_pop!=0) && (num_pop<num_real_pop)) //not reach the specified num of pop
		{
			for (i=0;i<file_Len;i++)
			{
				if (all_population_ID[i]>first_j)
				{
					all_population_ID[i]=all_population_ID[i]-1;
				}
				else
				{
					if (all_population_ID[i]==first_j)
					{
						all_population_ID[i]=first_i;
					}
				}
			}

			num_real_pop--;

		

		}
	}
		//printf("num_real_pop is %d\n",num_real_pop);
	
  } //end of while

  

  free(center_1);
  free(center_2);
  free(center_3);
  printf("Final number of populations is %ld\n",num_real_pop);

  new_population_center=(double **)malloc(sizeof(double*)*num_real_pop);
  memset(new_population_center,0,sizeof(double*)*num_real_pop);
  for (i=0;i<num_real_pop;i++)
  {
      new_population_center[i]=(double*)malloc(sizeof(double)*num_dm);
      memset(new_population_center[i],0,sizeof(double)*num_dm);
  }
  
 
  ///////////////////////////////////////////
  //End of population mapping
  ///////////////////////////////////////////

  show(input_data, all_population_ID, file_Len, num_real_pop, num_dm, para_name_string);

  ID2Center_all(input_data,file_Len,num_dm,num_real_pop,all_population_ID,new_population_center);
  

  f_cid=fopen("population_id.txt","w");
  f_ctr=fopen("population_center.txt","w");
  f_out=fopen("coordinates.txt","w");
  f_results=fopen("flock_results.txt","w");

/*
  f_parameters=fopen("parameters.txt","w");
  fprintf(f_parameters,"Number_of_Bins\t%d\n",num_bin);
  fprintf(f_parameters,"Density\t%f\n",aver_index);
  fclose(f_parameters);
*/

  for (i=0;i<file_Len;i++)
	fprintf(f_cid,"%ld\n",all_population_ID[i]+1); //all_population_ID[i] changed to all_population_ID[i]+1 to start from 1 instead of 0: April 16, 2009

  /*
   * New to check for min/max to add to parameters.txt
   *
  */
  
  fprintf(f_out,"%s\n",para_name_string);
  //fprintf(f_results,"%s\tEvent\tPopulation\n",para_name_string);
  fprintf(f_results,"%s\tPopulation\n",para_name_string);
  for (i=0;i<file_Len;i++)
  {
	for (j=0;j<num_dm;j++)
	{
		if (input_data[i][j] < min) {
			min = (int)input_data[i][j];
		}
		if (input_data[i][j] > max) {
			max = (int)input_data[i][j];
		}
		if (j==num_dm-1)
		{
			fprintf(f_out,"%d\n",(int)input_data[i][j]);
			fprintf(f_results,"%d\t",(int)input_data[i][j]);
		}
		else
		{
			fprintf(f_out,"%d\t",(int)input_data[i][j]);
			fprintf(f_results,"%d\t",(int)input_data[i][j]);
		}
	}
	//fprintf(f_results,"%ld\t",i + 1);
	fprintf(f_results,"%ld\n",all_population_ID[i]+1); //all_population_ID[i] changed to all_population_ID[i]+1 to start from 1 instead of 0: April 16, 2009
  }

/*
  f_parameters=fopen("parameters.txt","w");
  fprintf(f_parameters,"Number_of_Bins\t%ld\n",num_bin);
  fprintf(f_parameters,"Density\t%d\n",den_t_event);
  fprintf(f_parameters,"Min\t%d\n",min);
  fprintf(f_parameters,"Max\t%d\n",max);
  fclose(f_parameters);
*/

  f_properties=fopen("fcs.properties","w");
  fprintf(f_properties,"Bins=%ld\n",num_bin);
  fprintf(f_properties,"Density=%d\n",den_t_event);
  fprintf(f_properties,"Min=%d\n",min);
  fprintf(f_properties,"Max=%d\n",max);
  fprintf(f_properties,"Populations=%ld\n",num_real_pop);
  fprintf(f_properties,"Events=%ld\n",file_Len);
  fprintf(f_properties,"Markers=%ld\n",num_dm);
  fclose(f_properties);

  for (i=0;i<num_real_pop;i++) {
	/* Add if we want to include population id in the output
	*/
	fprintf(f_ctr,"%ld\t",i+1);  //i changed to i+1 to start from 1 instead of 0: April 16, 2009

	for (j=0;j<num_dm;j++) {
		if (j==num_dm-1)
			fprintf(f_ctr,"%.0f\n",new_population_center[i][j]);
		else
			fprintf(f_ctr,"%.0f\t",new_population_center[i][j]);
	}
  }

  	//added April 16, 2009
	f_mfi=fopen("MFI.txt","w");

	for (i=0;i<num_real_pop;i++)
	{
		fprintf(f_mfi,"%ld\t",i+1);

		for (j=0;j<num_dm;j++)
		{
			if (j==num_dm-1)
				fprintf(f_mfi,"%.0f\n",new_population_center[i][j]);
			else
				fprintf(f_mfi,"%.0f\t",new_population_center[i][j]);
		}
	}
	fclose(f_mfi);

	//ended April 16, 2009
			
  fclose(f_cid);
  fclose(f_ctr);
  fclose(f_out);
  fclose(f_results);


  for (i=0;i<num_population;i++)
  	free(population_center[i]);
  
  free(population_center);
 
  for (i=0;i<num_real_pop;i++)
	  free(new_population_center[i]);
  
  free(new_population_center);

  for (i=0;i<file_Len;i++)
    free(normalized_data[i]);
  free(normalized_data);	
	
  free(grid_populationID);

  free(cluster_populationID);
  free(grid_clusterID);
  free(cluster_ID);

  for (i=0;i<file_Len;i++)
    free(input_data[i]);
  free(input_data);

  free(grid_ID);
  free(population_ID);
  free(all_population_ID);
  free(eventID_To_denseventID);
		
  ///////////////////////////////////////////////////////////
  printf("Ending time:\t\t\t\t");
  fflush(stdout);
  system("/bin/date");

  /*
   * Windows version
  _strtime( tmpbuf );
  printf( "Ending time:\t\t\t\t%s\n", tmpbuf );
  _strdate( tmpbuf );
  printf( "Ending date:\t\t\t\t%s\n", tmpbuf );
 */
  
  return 0;
}
