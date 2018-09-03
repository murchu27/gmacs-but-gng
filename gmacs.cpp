//include directives
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include "gmacs.h"

//initialise parameters for reading motifs/sequences
FILE* out;
FILE* outR;
int DN;
char file[FN];
double **kfv_f, **kfv_r;

//initialise parameters for kfv forming & clustering
int pop=POP, gen=GEN;
double mut=MUTATE, trim=TRIM, rep=REP;
int trim_i=0, mut_i=0;
int min_motif=MIN_MOTIF;
int kmer=KMER, lmer=LMER, seql=SEQ_L;
int kfv_size=pow(4,kmer);

void print_msg(){
	printf("\n\n\t*******************************************************\n");
	printf("\t**                                                   **\n");
	printf("\t**              Welcome to GMACS (v1.2)              **\n");
	printf("\t**    Motif/Sequence Clustering and Alignment Tool   **\n");
	printf("\t**                                                   **\n");
	printf("\t**                                                   **\n");
	printf("\t**   Comments/Bugs to pilib.obroin@einstein.yu.edu   **\n");
	printf("\t**                                                   **\n");
	printf("\t*******************************************************\n\n");
}

void usage(){
	print_msg();
	printf("\tUsage:\n\n");
	printf("\t-h\t\tDisplay this help message.\n");
	printf("\t-f [file]\tMotif input filename.\n");
	printf("\t-p [pop]\tSet population size to pop. Min:%d (Default:%d)\n", POP, POP);
	printf("\t-g [gen]\tSet generation number to gen. Min:%d (Default:%d)\n", GEN, GEN);
	printf("\t-m [mutate]\tSet mutation rate to mutate. Valid range (0 -> 0.2) (Default:%.2lf)\n", MUTATE);
	printf("\t-t [trim]\tSet trim threshold to trim. Valid range (0 -> 0.4) (Default:%.1lf)\n", TRIM);
	printf("\t-s [seql]\tSet sequence length to seql. (Default:%.1lf)\n", SEQ_L);
	printf("\n\n\n");
}

//method for reading in files to newgmacs, and determining whether they are motif files, or sequence files
void get_args(int arg_C, char* arg_V[])
{

	//if no command line options provided, print a usage message, and exit
    if(arg_C==1)
    {
		usage();
        exit(0);
    }

	//otherwise, iterate over each of the arguments provided
    for(int j=1; j<arg_C; j++)
    {
        int count=0; //counts number of valid arguments
        char valid_args[6][3]= {"-h", "-f", "-p", "-g", "-m", "-t"}; //specifies options that are valid
        for(int i=0; i<6; i++)
            if((strcmp(arg_V[j], valid_args[i])==0)||(strcmp(arg_V[j-1], valid_args[i])==0))
                count++; //increment if the current argument is a valid option 
        if(count==0)
        {
            printf("\nInvalid option '%s'...Aborting\nTry 'gmacs -h' for usage.\n\n", arg_V[j]);
            exit(0);
        }
        else
        {
            if(strcmp(arg_V[j], "-h")==0)
            {//print usage message, then exit
                usage();
                exit(0);
            }

            if(strcmp(arg_V[j], "-f")==0)
            {//user wishes to specify a file for gmacs to work with
                if(arg_C>=j+2 && arg_V[j+1]!=NULL && strcmp(arg_V[j+1], " ")!=0)
                {//ensure that there is another argument after '-f', and that it is not empty
                	//copy this argument into 'file' variable, then call read_file 
                    strcpy(file, arg_V[j+1]);
                    read_file(file);
                }
                else
                {
                	//if no file is specified, cannot proceed
                    printf("\nNo motif file specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                    exit(0);
                }
            }

            if((strcmp(arg_V[j], "-p"))==0)
            {//user wishes to specify population size
                if(arg_C>=j+2 && arg_V[j+1]!=NULL && strcmp(arg_V[j+1], " ")!=0)
                {//ensure that there is another argument after '-p', and that it is not empty
                    //give value of that argument to 'pop' variable
                    pop=atoi(arg_V[j+1]);
                    if(pop<POP)
                    {//invalid population, cannot proceed
                        printf("\nInvalid population size specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                        exit(0);
                    }
                }
                else
                {
                    //no population specified, cannot proceed
                    printf("\nNo population size specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                    exit(0);

                }
            }

            if((strcmp(arg_V[j], "-g"))==0)
            {//user wishes to specify generation number
                if(arg_C>=j+2 && arg_V[j+1]!=NULL && strcmp(arg_V[j+1], " ")!=0)
                {//ensure that there is another argument after '-g', and that it is not empty
                    //give value of that argument to 'gen' variable
                    gen=atoi(arg_V[j+1]);
                    if(gen<GEN)
                    {//invalid generation, cannot proceed
                        printf("\nInvalid generation number specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                        exit(0);
                    }
                }
                else
                {
                    //no generation specified, cannot proceed
                    printf("\nNo generation number specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                    exit(0);
                }
            }

            if((strcmp(arg_V[j], "-m"))==0)
            {//user wishes to specify mutation rate
                if(arg_C>=j+2 && arg_V[j+1]!=NULL && strcmp(arg_V[j+1], " ")!=0)
                {//ensure that there is another argument after '-m', and that it is not empty
                    //give value of that argument to 'mut' variable
                    mut=atof(arg_V[j+1]);
                    mut_i=(floor(mut*100));
                    mut=(double)mut_i/100;
                    if(mut_i<0||mut_i>20)
                    {//invalid mutation, cannot proceed
                        printf("\nInvalid mutation rate specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                        exit(0);
                    }
                }
                else
                {
                    //no mutation specified, cannot proceed
                    printf("\nNo mutation rate specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                    exit(0);
                }
            }

            if((strcmp(arg_V[j], "-t"))==0)
            {//user wishes to specify trim threshold
                if(arg_C>=j+2 && arg_V[j+1]!=NULL && strcmp(arg_V[j+1], " ")!=0)
                {//ensure that there is another argument after '-t', and that it is not empty
                    //give value of that argument to 'trim' variable
                    trim=atof(arg_V[j+1]);
                    trim_i=(floor(trim*100));
                    trim=(double)trim_i/100;
                    if(trim_i<0||trim_i>40)
                    {//invalid trim, cannot proceed
                        printf("\nInvalid trim threshold specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                        exit(0);
                    }
                }
                else
                {
                    //no trim specified, cannot proceed
                    printf("\nNo trim threshold specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                    exit(0);
                }
            }
            
            if((strcmp(arg_V[j], "-s"))==0)
            {//user wishes to specify sequence length
                if(arg_C>=j+2 && arg_V[j+1]!=NULL && strcmp(arg_V[j+1], " ")!=0)
                {//ensure that there is another argument after '-s', and that it is not empty
                    //give value of that argument to 'seql' variable
                    seql=atoi(arg_V[j+1]);
                    if(seql<2)
                    {//invalid sequence length, cannot proceed
                        printf("\nInvalid sequence length specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                        exit(0);
                    }
                }
                else
                {
                    //no sequence length specified, cannot proceed
                    printf("\nNo sequence length specified...Aborting\nTry 'gmacs -h' for usage.\n\n");
                    exit(0);

                }
            }
        }
    }
}

int main(int arg_C, char* arg_V[])
{
	//open file for reporting kfvs/clusters, check for errors
	out=fopen("report.txt", "w");
    if (out==NULL)
        exit(0);

	outR=fopen("reportR.txt", "w");
    if (outR==NULL)
    	exit(0);
	
	//pass command line arguments to get_args for handling
    get_args(arg_C, arg_V);
	
	//get_args then calls other methods which populate 'kfv_f' and 'kfv_r' 
	//it also prints the contents of these kfv arrays to the report files, so we are finished with them, and can close them
    fclose(out);
    fclose(outR);
    
    
    //now, we have our kfvs!

    return (0);

}

