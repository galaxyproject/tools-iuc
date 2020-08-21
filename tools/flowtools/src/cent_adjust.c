/////////////////////////////////////////////////////////
//  Cent_adjust version number and modification history
//  ImmPort BISC project
//  Author: Yu "Max" Qian
//  v1.01: Oct 16, 2009
//         Line 899 of the main function:
//         Changed kmean_term=1 to kmean_term=2
//////////////////////////////////////////////////////////


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DEBUG 0
#define LINE_LEN 1024
#define FILE_NAME_LEN 128
#define PARA_NAME_LEN 64
#define MAX_VALUE 1000000000
#define CUBE 0


void getctrfileinfo(FILE *f_src_ctr, long *num_clust)
{
	int ch='\n';
	int prev='\n';
	long num_rows=0;

	while ((ch = fgetc(f_src_ctr))!= EOF )
    {
		if (ch == '\n')
        {
			++num_rows;
        }
		prev = ch;
    }
	if (prev!='\n')
		++num_rows;
	
	*num_clust=num_rows;
	//printf("center file has %ld rows\n", *num_clust);
}

/************************************* Read basic info of the source file **************************************/
void getfileinfo(FILE *f_src, long *file_Len, long *num_dm, char *name_string, int *time_ID)
{
	char src[LINE_LEN];
	char current_name[64];
	char prv;

	long num_rows=0;
	long num_columns=0;
	int ch='\n';
	int prev='\n';
	long time_pos=0;
	long i=0;
	long j=0;

	

	src[0]='\0';
	fgets(src, LINE_LEN, f_src);

	name_string[0]='\0';
	current_name[0]='\0';
	prv='\n';

	while ((src[i]==' ') || (src[i]=='\t')) //skip space and tab characters
		i++;

	while ((src[i]!='\r') && (src[i]!='\n')) //repeat until the end of the line
	{
		current_name[j]=src[i];
		
		if ((src[i]=='\t') && (prv!='\t')) //a complete word
		{
			current_name[j]='\0';
			
          /* 
           * Commented out John Campbell, June 10 2010
           * We no longer want to automatically remove Time column.
           * This column should have been removed by column selection
          if (0!=strcmp(current_name,"Time"))
            {
              num_columns++; //num_columns does not inlcude the column of Time
              time_pos++;
              strcat(name_string,current_name); 
              strcat(name_string,"\t");
            }
          else
            {
              *time_ID=time_pos;
            }
          */

           num_columns++;
           strcat(name_string,current_name);
           strcat(name_string,"\t");

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
	
	if (prv!='\t') //the last one hasn't been retrieved
	{
		current_name[j]='\0';
      /* 
       * Commented out John Campbell, June 10 2010
       * We no longer want to automatically remove Time column.
       * This column should have been removed by column selection
      if (0!=strcmp(current_name,"Time"))
        {
          num_columns++;
          strcat(name_string,current_name);
          time_pos++;
        }
      else
        {
          *time_ID=time_pos;
        }
      */

      num_columns++;
      strcat(name_string,current_name);
	}

	if (DEBUG==1)
	{
		printf("time_ID is %d\n",*time_ID);
		printf("name_string is %s\n",name_string);
	}

	// # of rows

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
	
	*file_Len=num_rows;
	*num_dm=num_columns; 

	//printf("original file size is %ld; number of dimensions is %ld\n", *file_Len, *num_dm);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/************************************* Read the source file into uncomp_data **************************************/
void readsource(FILE *f_src, long file_Len, long num_dm, double **uncomp_data, int time_ID)
{
	long time_pass=0; //to mark whether the time_ID has been passed
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
		time_pass=0;
						
		if (time_ID==-1)
		{
			for (t=0;t<num_dm;t++) //there is no time_ID
			{
				xc[0]='\0';
				j=0;
				while ((src[i]!='\r') && (src[i]!='\n') && (src[i]!=' ') && (src[i]!='\t'))
				{
					xc[j]=src[i];
					i++;
					j++;
				}
		
				xc[j]='\0';	    
				i++;

				uncomp_data[index][t]=atof(xc);
			}	
		}
		else
		{
			for (t=0;t<=num_dm;t++) //the time column needs to be skipped, so there are num_dm+1 columns
			{
				xc[0]='\0';
				j=0;
				while ((src[i]!='\r') && (src[i]!='\n') && (src[i]!=' ') && (src[i]!='\t'))
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
		}        	
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readcenter(FILE *f_src_ctr, long num_clust, long num_dm, double **cluster_center, long *IDmapping)
{
	char src[LINE_LEN];
	char xc[LINE_LEN/10];

	long i=0;
	long j=0;
	int m=0;
	int t=0;

	for (i=0;i<num_clust;i++)
	{
		src[0]='\0';
		fgets(src,LINE_LEN, f_src_ctr); 
		m=0;
		for (j=0;j<num_dm+1;j++)
		{
			xc[0]='\0';
			t=0;
			while ((src[m]!='\r') && (src[m]!='\n') && (src[m]!=' ') && (src[m]!='\t'))
			{
				xc[t]=src[m];
				m++;
				t++;
			}
			xc[t]='\0';	    
			m++;
			if (j==0)
				IDmapping[i]=atoi(xc);
			else
				cluster_center[i][j-1]=atof(xc);
			//printf("cluster_center[%d][%d]=%f\n",i,j,cluster_center[i][j]);
		}
	}
}


/**************************************** Normalization ******************************************/
void tran(double **orig_data, long clean_Len, long num_dm, long norm_used, double **matrix_to_cluster)
{
	long i=0;
	long j=0;

	double biggest=0;
	double smallest=MAX_VALUE;

	double *aver; //average of each column
	double *std; //standard deviation of each column

	aver=(double*)malloc(sizeof(double)*clean_Len);
	memset(aver,0,sizeof(double)*clean_Len);

	std=(double*)malloc(sizeof(double)*clean_Len);
	memset(std,0,sizeof(double)*clean_Len);	
		
	if (norm_used==2) //z-score normalization
	{
		for (j=0;j<num_dm;j++)
		{
			aver[j]=0;
			for (i=0;i<clean_Len;i++)
				aver[j]=aver[j]+orig_data[i][j];
			aver[j]=aver[j]/(double)clean_Len;

			std[j]=0;
			for (i=0;i<clean_Len;i++)
				std[j]=std[j]+(orig_data[i][j]-aver[j])*(orig_data[i][j]-aver[j]);
			std[j]=sqrt(std[j]/(double)clean_Len);
			
			for (i=0;i<clean_Len;i++)
				matrix_to_cluster[i][j]=(orig_data[i][j]-aver[j])/std[j];  //z-score normalization
		}
	}

	if (norm_used==1) //0-1 min-max normalization
	{
		for (j=0;j<num_dm;j++)
		{
			biggest=0;
			smallest=MAX_VALUE;
			for (i=0;i<clean_Len;i++)
			{
				if (orig_data[i][j]>biggest)
					biggest=orig_data[i][j];
				if (orig_data[i][j]<smallest)
					smallest=orig_data[i][j];
			}
			
			for (i=0;i<clean_Len;i++)
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
		for (i=0;i<clean_Len;i++)
			for (j=0;j<num_dm;j++)
				matrix_to_cluster[i][j]=orig_data[i][j];
	}

	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void assign_event(double **Matrix, long k, long dist_used, double kmean_term, long file_Len, long num_dm, long *shortest_id, double **center, int random_init)
{
	
	long i=0;

	long j=0;
	long t=0;
	long random=0;
	long random1=0;
	long random2=0;
	long times=0;
	long times_allowed=0;
	
 	long *num;  //num[i]=t means the ith cluster has t points
	
	double vvv=1.0; // the biggest variation;
	double distance=0.0;
	double xv=0.0;
	double variation=0.0;
	double EPS=0.0;
	double diff=0.0;
	double mean_dx=0;
	double mean_dy=0;
	double sum_var=0;
	double dx=0;
	double dy=0;
	double sd_x=0;
	double sd_y=0;	
	
	double *temp_center;	
	double *shortest_distance;
			
	double **sum;	
 
	temp_center = (double *)malloc(sizeof(double)*num_dm);
	memset(temp_center,0,sizeof(double)*num_dm);

	/* Choosing Centers */
	if (random_init)
	{
		for (i=0;i<k;i++)
		{	
			random1=rand()*rand();
			//srand( (unsigned)time( NULL ) );
			random2=abs((random1%5)+1);
			for (t=0;t<random2;t++)
				random2=random2*rand()+rand();
	
			random=abs(random2%file_Len);
			//printf("random=%d\n",random);
			for (j=0;j<num_dm;j++)
				center[i][j]=Matrix[random][j];			
			
		}
	}

	//printf("finish random selection\n");
	/* To compute the nearest center for every point */

	shortest_distance = (double *)malloc(sizeof(double)*file_Len);
	memset(shortest_distance,0,sizeof(double)*file_Len);

	num = (long *)malloc(sizeof(long)*k);
	memset(num,0,sizeof(long)*k);

	sum = (double **)malloc(sizeof(double*)*k);
	memset(sum,0,sizeof(double*)*k);
	for (i=0;i<k;i++)
	{
		sum[i] = (double *)malloc(sizeof(double)*num_dm);
		memset(sum[i],0,sizeof(double)*num_dm);
	}

	for (i=0;i<k;i++)
		for (j=0;j<num_dm;j++)
			sum[i][j]=0.0;        //sum[i][j] = k means the sum of the jth dimension of all points in the ith group is k 

	//printf("before recursion\n");
	if (kmean_term>=1)
		times_allowed = (long)kmean_term;
	else
		EPS = kmean_term;

	times=0;

	while (((vvv>EPS) && (kmean_term<1)) || ((times<times_allowed) && (kmean_term>=1)))
	{
		for (i=0;i<k;i++)
		{
			num[i]=0;
			for (j=0;j<num_dm;j++)
				sum[i][j]=0.0;  
		}
		
		for (i=0;i<file_Len;i++)  //for each data point i, we compute the distance between Matrix[i] and center[j]
		{		
			shortest_distance[i]=MAX_VALUE;
			for (j=0;j<k;j++)  //for each center j
			{
				
				distance=0.0;
								
				if (dist_used==0)  //Euclidean distance
				{
					for (t=0;t<num_dm;t++) //for each dimension
					{
						diff=center[j][t]-Matrix[i][t];
						if (diff<0)
							diff=-diff;

						if (CUBE)
							distance = distance+(diff*diff*diff);	
						else
							distance = distance+(diff*diff);
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
				//printf("mean_dx=%f\n",mean_dx);

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
				

				if (distance<shortest_distance[i])
				{
					shortest_distance[i]=distance;						
					shortest_id[i]=j;
				}			
			}//end for j
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
	free(shortest_distance);
		
}

//////////////////////////////////////////////////////
/*************************** Show *****************************/
void show(double **Matrix, long *cluster_id, long file_Len, long k, long num_dm, char *name_string, long *IDmapping)
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
		fprintf(fprof_id,"%ld\t",IDmapping[i]);
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
		fprintf(fpcnt_id,"%ld\t%.2f\n",IDmapping[t],(double)size_c[t]*100.0/(double)file_Len);										
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
/******************************************************** Main Function **************************************************/
int main (int argc, char **argv)
{
	//inputs
	FILE *f_src; //source file pointer
	FILE *f_src_ctr; //source center file
	//outputs
	FILE *f_cid; //cluster-id file pointer
	FILE *f_mfi; //added April 16, 2009
	

	char name_string[LINE_LEN]; //for name use

	int time_id=-1;

	long file_Len=0;
	long num_clust=0;
	long num_dm=0;
	long norm_used=0;
	long dist_used=0;
	long i=0;
	long j=0;

	long *cluster_id;
	long *IDmapping; //this is to keep the original populationID of the center.txt

	double kmean_term=0;
	
	double **cluster_center;
	double **orig_data;
	double **normalized_data;
		
/*
	_strtime( tmpbuf );
    printf( "Starting time:\t\t\t\t%s\n", tmpbuf );
	_strdate( tmpbuf );
    printf( "Starting date:\t\t\t\t%s\n", tmpbuf );
*/

	if (argc!=3)
	{
		printf("usage: cent_adjust input_center input_data_file\n");       
		exit(0);
	}	
	

	f_src_ctr=fopen(argv[1],"r");	
	
	//read source data
	f_src=fopen(argv[2],"r");
	
	getfileinfo(f_src, &file_Len, &num_dm, name_string, &time_id); //get the filelength, number of dimensions, and num/name of parameters

	rewind(f_src); //reset data file pointer	

	orig_data = (double **)malloc(sizeof(double*)*file_Len);
	memset(orig_data,0,sizeof(double*)*file_Len);
	for (i=0;i<file_Len;i++)
	{
		orig_data[i]=(double *)malloc(sizeof(double)*num_dm);
		memset(orig_data[i],0,sizeof(double)*num_dm);
	}
	
	readsource(f_src, file_Len, num_dm, orig_data, time_id); //read the data;
	
	fclose(f_src);
	/////////////////////////////////////////////////////////////////////////////
	getctrfileinfo(f_src_ctr, &num_clust); //get how many populations
	norm_used=0;
	dist_used=0;
	kmean_term=2;  //modified on Oct 16, 2009: changed kmean_term=1 to kmean_term=2

	rewind(f_src_ctr); //reset center file pointer

	//read population center
	cluster_center=(double **)malloc(sizeof(double*)*num_clust);
	memset(cluster_center,0,sizeof(double*)*num_clust);
	for (i=0;i<num_clust;i++)
	{
		cluster_center[i]=(double*)malloc(sizeof(double)*num_dm);
		memset(cluster_center[i],0,sizeof(double)*num_dm);
	}
	for (i=0;i<num_clust;i++)
		for (j=0;j<num_dm;j++)
			cluster_center[i][j]=0;

	IDmapping=(long *)malloc(sizeof(long)*num_clust);
	memset(IDmapping,0,sizeof(long)*num_clust);

	readcenter(f_src_ctr,num_clust,num_dm,cluster_center,IDmapping); //read population center
    fclose(f_src_ctr);

	/////////////////////////////////////////////////////////////////////////////
	normalized_data=(double **)malloc(sizeof(double*)*file_Len);
	memset(normalized_data,0,sizeof(double*)*file_Len);
	for (i=0;i<file_Len;i++)
	{
		normalized_data[i]=(double *)malloc(sizeof(double)*num_dm);
		memset(normalized_data[i],0,sizeof(double)*num_dm);
	}
	
	tran(orig_data, file_Len, num_dm, norm_used, normalized_data);
	/************************************************* Compute number of clusters *************************************************/
	
	cluster_id=(long*)malloc(sizeof(long)*file_Len);
	memset(cluster_id,0,sizeof(long)*file_Len);

	assign_event(normalized_data,num_clust,dist_used,kmean_term,file_Len,num_dm,cluster_id,cluster_center,0);

	
	//show(orig_data,cluster_id,file_Len,num_clust,num_dm,show_data,num_disp,name_string); 
	show(orig_data, cluster_id, file_Len, num_clust, num_dm, name_string, IDmapping);

	f_cid=fopen("population_id.txt","w");

	for (i=0;i<file_Len;i++)
		fprintf(f_cid,"%ld\n",IDmapping[cluster_id[i]]);
		

	fclose(f_cid);
 
	//added April 16, 2009
	f_mfi=fopen("MFI.txt","w");

	for (i=0;i<num_clust;i++)
	{
		fprintf(f_mfi,"%ld\t",IDmapping[i]);

		for (j=0;j<num_dm;j++)
		{
			if (j==num_dm-1)
				fprintf(f_mfi,"%.0f\n",cluster_center[i][j]);
			else
				fprintf(f_mfi,"%.0f\t",cluster_center[i][j]);
		}
	}
	fclose(f_mfi);

	//ended April 16, 2009

	for (i=0;i<num_clust;i++)
		free(cluster_center[i]);
	free(cluster_center);
		

	/********************************************** Release memory ******************************************/
  
	for (i=0;i<file_Len;i++)
	{
		free(orig_data[i]);		
		free(normalized_data[i]);
	}
	
	free(orig_data);
	free(normalized_data);
	free(cluster_id);
	free(IDmapping);

/*
	_strtime( tmpbuf );
    printf( "Ending time:\t\t\t\t%s\n", tmpbuf );
	_strdate( tmpbuf );
    printf( "Ending date:\t\t\t\t%s\n", tmpbuf );
*/

}
