#include "gmacs.h"

//method for determining type of file specified by user; a motif file, or a peak sequence file
void read_file(char* file)
{

    /* 
    
    mot_file format:
    >AGL3_MADS
    0.000	0.969	0.010	0.021
    0.031	0.773	0.000	0.196
    ...(4xl frequency values)
    >ARNT_bHLH
    0.200	0.800	0.000	0.000
    0.950	0.000	0.050	0.000
    ...

    seq_file format:
    >chr1_33274558_33277058
    AACCTGGGCCATAAAGCTAATTGCCTCCATTCAGCCTGTAATCTAGATGC
    CTCGGCCCAGGGCTCACTCTACCCCTCAAACCTTTCCAGTCTCTAGACCC
    ...(50x50 chars per chromosome)
    >chr1_150232838_150235338
    GGGCGGAGCTCACCTTGGCCGAGGCGCGGCGGACGCTGGGCGAGCTGGGC
    ...
	
	*/

	//regardless of file, we will need the following
    int c=0; //for reading in chars from file
    int count=0; //number of mot_file motifs OR number of seq_file sequences
    FILE* fp; //FILE object for opening the file containing either motifs or sequences
    char name[FN]; //a c-string to store the names of motifs/sequences

	//open file specified by user, check that it is a valid file
    fp=fopen(file, "r");
    if(fp == NULL)
    {//exit if there is an issue with file
        printf("\nError reading from file '%s'...Aborting.\nCheck file exists.\n\n", file);
        exit(0);
    }

    //count motifs/sequences
    //since each motif/sequence name is printed with a leading '>' (and there is no other reason to use a '>'...), 
    //we can count the number of motifs/sequences by counting the number of '>' characters
    do
    {//simply iterate over the characters in the file and count how many times '>' appears
        c = fgetc(fp);
        if(c=='>')
            count++;
    }
    while(c!=EOF);

	//if '>' does not appear, the file is not a valid motif/sequence file, or there is some other problem with reading
    if(count==0)
    {
        printf("\nError reading from file '%s' or no motifs/sequences found...Aborting.\nPlease see readme for acceptable formats.\n\n", file);
        exit(0);
    }

	//copy value of count to global variable 'DN'
    DN = count;

	//exit if there are 2 or fewer motifs/sequences
    if (DN<=2)
    {
        printf("\n2 or fewer motifs/sequences found.\n\n");
        exit(0);
    }

    fseek(fp, 0, SEEK_SET); //goes back to start of file
    do //get to start of name
    {
        c = fgetc(fp);
    }
    while(c!='>' && c!=EOF);


    //get name, verify no errors
    if(fgets(name, FN, fp)==NULL)
        printf("\nError reading from file '%s'...Aborting.\nPlease see readme for acceptable formats.\n\n", file);


    //skip control chars
    do
    {
        c = fgetc(fp);
    }
    while(c==32 || c==10);


    //check first character to determine if this is a motif or a sequence, and pass file to relevant method
    //if nucleotides are involved, it's a sequence file; otherwise, assume it is a motif file
    if (c=='A'||c=='C'||c=='G'||c=='T')//sequence
        seq_read(fp);
    else//motif
        mot_read(fp);

}

//method for reading motif data from a given file, and storing it in motif objects
void mot_read(FILE* mot)
{
	/*
	
	have ints for:
	 ---reading characters ('c') 
	 ---measuring the length (in nucleotides) of motifs ('len') 
	 ---keeping track of the motif object that we are working with ('index')
	
	*/
	 
	int c=0, len=0, index=0;
    char name[FN]; //names of motifs
    double pr[4]; //the probabilities of nucs at each position
	
	//create array of motif objects, one for each motif (i.e., DN objects)
	motif* motifs = (motif*) malloc(DN*sizeof(struct motif));
	//motifs = (motif*) malloc(DN*sizeof(struct motif));

    //now we need to read in PSSMs for each motif
    fseek(mot, 0, SEEK_SET); //goes back to start of file

    //motif name should start with '>'
    c = fgetc(mot);
    while(c!=EOF)
    {
        //keep getting chars until start of motif name is reached
        while(c!='>' && c!=EOF)
            c = fgetc(mot);

        //get name of motif, ensuring that a name exists (and exiting if it does not) 
        if(fgets(name, FN, mot)==NULL){
        	printf("\nError reading from file '%s'...Aborting.\nPlease see readme for acceptable formats.\n\n", file);
        	exit(0);
        }
        
       	//strip newline char, replace with terminating null character
        char *pos;
        if ((pos=strchr(name, '\n')) != NULL)
            *pos = '\0';
            
        //copy into name attribute of current motif
        strcpy (motifs[index].name, name);

        //skip control chars to get to PSSM
        do{
        	c = fgetc(mot);
        }while(c==32 || c==10);

		//this is a new motif; set len to zero before starting
        len = 0;
        
        //iterate until the end of this motif (at which point, this motif's PSSM will be filled)
        while(c!='>' && c!=EOF)
        {
            //seek back one char for the one just removed with fgetc() when skipping control chars
            fseek(mot, -1, SEEK_CUR);
            
            //copy probability scores into array, alert if not all scores are read properly
            if(fscanf(mot, "%lf %lf %lf %lf", &pr[0], &pr[1], &pr[2], &pr[3])==0){
                printf("\nError reading from file '%s'...Aborting.\nPlease see readme for acceptable formats.\n\n", file);
                exit(0);
			}
			
			//probability scores read successfully, so now copy them into the PSSM of the current motif
            for(int j=0; j<4; j++)
                motifs[index].pssm[len][j] = pr[j];
            len++; //move on to next line of PSSM

            //ignore control chars
            do{
	            c = fgetc(mot);
            }while(c==32 || c==10);
        }

		//ensure motif is not larger than permitted
        if(len>MAX_L){
			printf("\nError reading from file '%s'; motif '%s' is too long!\nMaximum length motif permitted is %i...Aborting.\n\n", file, motifs[index].name, MAX_L);
            exit(0);
        }

		//copy len into length attribute of current motif
        motifs[index].length=len;
        
        //this motif is complete; move on to the next
        index++;
    }

	//array of motifs is complete, pass along to mot_to_kfv for conversion to kfv's
    mot_to_kfv(motifs);
}

//method for reading sequence data from a given file, and storing it in sequence objects
void seq_read(FILE* seq)
{
    /*
	
	have ints for:
	 ---reading characters ('c') 
	 ---keeping track of the sequence object that we are working with ('index')
	 ---keeping track of the char we are inserting into the 'seq' attribute of the current sequence object ('len')
	
	*/
    
    int c = 0, index = 0, len = 0;
    char name[FN]; //names of sequences

	//create array of sequence objects, one for each sequence (i.e., DN objects)
	sequence* sequences;
    sequences = (sequence*) malloc(DN*sizeof(struct sequence));

    //read in sequences
    fseek(seq, 0, SEEK_SET); //goes back to start of file
    
    //sequence name should start with '>'
    c = fgetc(seq);
	do{
        while(c!='>' && c!=EOF)
            c = fgetc(seq);

        //get name of sequence, ensuring that a name exists (and exiting if it does not) 
        if(fgets(name, FN, seq)==NULL){
        	printf("\nError reading from file '%s'...Aborting.\nPlease see readme for acceptable formats.\n\n", file);
        	exit(0);
        }
        
       	//strip newline char, replace with terminating null character
        char *pos;
        if ((pos=strchr(name, '\n')) != NULL)
            *pos = '\0';
        
        //copy into name attribute of current sequence
        strcpy (sequences[index].name, name);

        //skip control chars to get to nucleotides
        do{
            c = fgetc(seq);
        }while(c==32 || c==10);

		//this is a new sequence; set len to zero before starting 
        len = 0;

		//iterate until the end of the sequence is reached (at which point the seq attribute of this sequence object will be filled)
        do{
            //ensure that c is a nucleotide
            if (!(c==65||c==67||c==71||c==84||c==97||c==99||c==103||c==116)){
                printf("\nError reading from file '%s'...Aborting.\nPlease see readme for acceptable formats.\n\n", file);
                exit(0);
            }

			//add it to this sequence's seq attribute
			sequences[index].seq[len] = (char)c;
            len++;//move to next char of seq

            //ignore control chars
            do{
                c = fgetc(seq);
            }while(c==32 || c==10);
            
        }while(c!='>' && c!=EOF);
        
        //this sequence is complete; move on to the next
        index++;
    
    }while(c!=EOF);

	//array of sequences is complete, pass along to seq_to_kfv for conversion to kfv's
    seq_to_kfv(sequences);
}

//NOT FAMILIAR WITH tHIS METHOD, REFER TO PILIB
void trim_motif(struct motif* motif_x)
{

    int st_pos=0, end_pos=0, new_len=0;

    /*first, calculate the information content*/
    get_ic(motif_x);

    st_pos=0;
    end_pos=(motif_x->length)-1;
    new_len=motif_x->length;

    /*do not trim beyond min_motif*/
    while(((motif_x->ic[st_pos]<trim_i)||(motif_x->ic[end_pos]<trim_i))&&new_len>min_motif)
    {

        if(motif_x->ic[st_pos]<motif_x->ic[end_pos])
        {
            st_pos++;
            new_len--;
        }
        else
        {
            end_pos--;
            new_len--;
        }
    }

    for(int i=0; i<new_len; i++)
    {
        motif_x->ic[i]=motif_x->ic[i+st_pos];
        for(int j=0; j<4; j++)
            motif_x->pssm[i][j]=motif_x->pssm[i+st_pos][j];
    }

    motif_x->length=new_len;
}

//NOT FAMILIAR WITH THIS METHOD, REFER TO PILIB
void get_ic(struct motif* motif_x)
{

    double sum=0.0;

    for(int y=0; y<MAX_L; y++)
        motif_x->ic[y]=0.0;

    /*get IC for each col*/
    for(int i=0; i<motif_x->length; i++)
    {
        sum=0.0;
        for(int j=0; j<4; j++)
            if(motif_x->pssm[i][j]>0.0)
                sum+=motif_x->pssm[i][j]*log2(motif_x->pssm[i][j]);

        motif_x->ic[i]=2+sum;
    }
}

//method for converting a given array of motif objects into kfv's
void mot_to_kfv(motif* motifs)
{
    /*
	
	 *	ii,jj,kk,ll: keeps track of pssm score for the individual nucleotides that make up the kmer whose frequency is currently being calculated
	 *	score: product of pssm scores for individual nucleotides at current step (i.e., probability of current kmer at current step)
	 *	sum: sum of 'score' values over all steps (i.e., probability of current kmer over all steps)
	 *	index: keeps track of the current kmer in the kfv_f and kfv_r arrays
	 * 	len: holds value of motifs[x].length for current motif
	 *	rev_motif: used so that frequency of kmers in reverse direction can be calculated
	 *  sym: used for printing out the symbols of the current motif when reporting to output file
	
	*/
	
    double ii=0.0, jj=0.0, kk=0.0, ll=0.0;
    double score=0.0, sum=0.0;
    int index=0, len=0;
    motif rev_motif;    
    char sym[4] = {'A','C','G','T'};	

    //first, allocate memory for frequency vectors
    kfv_f = (double**) malloc((DN)*sizeof(double*));
    kfv_r = (double**) malloc((DN)*sizeof(double*));
	for(int x=0; x<DN; x++)
    {
        kfv_f[x] = (double*) malloc(kfv_size*sizeof(double));
        kfv_r[x] = (double*) malloc(kfv_size*sizeof(double));
        for(int y=0; y<kfv_size; y++){
            kfv_f[x][y] = 0.0;
            kfv_r[x][y] = 0.0;
        }
    }

	//iterate over each motif
    for(int x=0; x<DN; x++)
    {
		//each motif should be trimmed        
        trim_motif(&motifs[x]);
        
        //give len new value of motifs[x].length
        len=motifs[x].length;

        //set up rev_motif for this motif
        //rev_motif's contents should be made by reading the current motif in reverse, and by taking the complement of each nucleotide
        //e.g., normal motif: ACCTG, reverse motif: CAGGT
        for(int i=0; i<len; i++)
            for(int j=0; j<4; j++)
                rev_motif.pssm[i][j]=motifs[x].pssm[(len-1)-i][3-j];

		//before beginning to report kfv scores, output name of current motif for simplicity
        fprintf(out, "%i. %s\n", x+1, motifs[x].name);
        fprintf(out, "        FOR      REV\n");
        
        //calculate kfv for current motif (CURRENTLY ONLY k=4 SUPPORTED)
        //i tracks the first nucleotide of the kmer, j the second, etc.
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                for(int k=0; k<4; k++)
                {
                    for(int l=0; l<4; l++)
                    {
                    	//iterates over kmers in order AAAA,AAAC,...,AACA,...,TTTT
                    	//calculate correct index of kfv_f and kfv_r for this kmer
                    	index = i*pow(4,3) + j*pow(4,2) + k*4 + l;
                    	
                    	//print out symbols for current kmer
                    	fprintf(out, "%c%c%c%c   ", sym[i], sym[j], sym[k], sym[l]);
                    	
                    	//new kmer, reset sum to zero
                        sum=0;
                        
                        //check for probability of current kmer occurring at each step of the motif (stop at (len-kmer) to avoid memory access issues) 
                        for(int pos=0; pos<=(len-kmer); pos++)
                        {
                            ii=motifs[x].pssm[pos][i];
                            jj=motifs[x].pssm[pos+1][j];
                            kk=motifs[x].pssm[pos+2][k];
                            ll=motifs[x].pssm[pos+3][l];
                            score=ii*jj*kk*ll; //probability of current kmer at current step
                            sum+=score; //add to total probability of current kmer over entire motif
                        }
                        
                        //at end, store the value of sum in the kfv_f array
                        kfv_f[x][index]=sum;
                       	
                       	//also report it to output file
						fprintf(out, "%.5f    ", sum);

						//now do the same for the reverse motif, storing sum in the kfv_r array at the end
                        sum=0;
                        for(int pos=0; pos<=(len-kmer); pos++)
                        {
                            ii=rev_motif.pssm[pos][i];
                            jj=rev_motif.pssm[pos+1][j];
                            kk=rev_motif.pssm[pos+2][k];
                            ll=rev_motif.pssm[pos+3][l];
                            score=ii*jj*kk*ll;
                            sum+=score;
                        }
                        kfv_r[x][index]=sum;
						fprintf(out, "%.5f\n", sum); //move to next line for next kmer
                    }
                }
            }
        }
        //move onto the next motif
        fprintf(out, "\n\n");
    }//all motifs finished!
    
    //output file for R import
    fprintf(outR, "Name    Direction    ");
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                for(int l=0; l<4; l++)
                    fprintf(outR, "%c%c%c%c	   ", sym[i], sym[j], sym[k], sym[l]);
	fprintf(outR, "\n");
	
	for(int x=0; x<DN; x++){
		fprintf(outR, "%s    1   ", motifs[x].name);
		for(int y=0; y<kfv_size;y++)
            fprintf(outR, "%.5f ", kfv_f[x][y]);
        fprintf(outR, "\n");
        
        fprintf(outR, "%s    2   ", motifs[x].name);
		for(int y=0; y<kfv_size;y++)
            fprintf(outR, "%.5f ", kfv_r[x][y]);
        fprintf(outR, "\n");
    }
}

//method for converting a given array of sequence objects into kfv's
void seq_to_kfv(sequence* sequences)
{
    /*
	
	 *	strcnt: number of full length strings (of length 'lmer')
	 *	r: length of remainder string
	 *	kmertot: total number of kmers (excluding remainder string)
	 *	f_index/r_index: used to find the index of kfv_f/kfv_r corresponding to the kmer frequency score that should be changed  
	 * 	n: used for converting nucleotides to quaternary (0,1,2,3) to help find f_index/r_index
	 *	sym: used for printing out the symbols of the current motif when reporting to output file
	
	*/
    
    int strcnt = seql/lmer;
    int r = seql-(strcnt*lmer);
    int kmertot = (lmer-kmer+1)*strcnt;
    
    if (r>=kmer) //then there are more kmers in the remainder string, so add those too
        kmertot += (r-kmer+1);

    int f_index, r_index, n;
	char sym[4] = {'A','C','G','T'};

	//first, allocate memory for frequency vectors
    kfv_f = (double**) malloc((DN)*sizeof(double*));
    kfv_r = (double**) malloc((DN)*sizeof(double*));
	for(int x=0; x<DN; x++)
    {
        kfv_f[x] = (double*) malloc(kfv_size*sizeof(double));
        kfv_r[x] = (double*) malloc(kfv_size*sizeof(double));
        for(int y=0; y<kfv_size; y++){
            kfv_f[x][y] = 0.0;
            kfv_r[x][y] = 0.0;
        }
    }

    for(int x=0; x<DN; x++) //iterate over each sequence
    {
        
        //print out name of sequence to report file for simplicity
		fprintf(out, "%i. %s\n", x+1, sequences[x].name);
        fprintf(out, "         FOR        REV\n");
        
        for(int l=0; l<strcnt; l++) //iterate over each lmer in the sequence
        {
            for(int k=0; k<(lmer-kmer+1); k++) //iterate over kmers in each lmer
            {
                /*for(int b=0; b<kmer; b++) //iterate over chars in each kmer
                {
                	nucs[k] = sequences[x].seq[(l*lmer)+k+b];
                }//now nucs holds the kmer at this step
                //we need to add it to the kfv*/
                
                //new kmer; reset f_index and r_index 
                f_index = 0;
				r_index = 0;       
				
				//each time a kmer occurs in the sequence, we update the frequency score corresponding to that kmer in the kfv arrays
				//so, we read the kmer that is present at each step, and determine the index that points to the kmer in the kfv array
				//check each base in the current kmer, starting with the last
                for(int b=kmer-1; b>=0; b--)
                {
                    
                    switch (toupper(sequences[x].seq[(l*lmer)+k+b]))
                    {
                    case 'A':
                        n = 0;
                        break;
                    case 'C':
                        n = 1;
                        break;
                    case 'G':
                        n = 2;
                        break;
                    case 'T':
                        n = 3;
                        break;
                    }
					
					//in kfv array, frequency scores are stored for each kmer in order AAAA, AAAC, ..., AACA, ..., TTTT
					//e.g., AGTA lies at index (0 * 1)+(3 * 4)+(2 * 16)+(0 * 64)
                    f_index += n*int(pow(4,kmer-1-b));
                    r_index += n*int(pow(4,b));
                }

				//now add frequency of a kmer occurring across total number of kmers
                kfv_f[x][f_index] += (1.0/kmertot);
                kfv_r[x][r_index] += (1.0/kmertot);

            }//proceed to next kmer
        }//proceed to next lmer


		//if there are any kmers in the remainder string, we handle those the same way as the earlier kmers
        if (r>=kmer)
        {
            for(int k=0; k<(r-kmer+1); k++)
            {
                f_index = 0;
                r_index = 0;
                for(int b=kmer-1; b>=0; b--)
                {
                    switch (toupper(sequences[x].seq[(strcnt*lmer)+k+b]))
                    {
                    case 'A':
                        n = 0;
                        break;
                    case 'C':
                        n = 1;
                        break;
                    case 'G':
                        n = 2;
                        break;
                    case 'T':
                        n = 3;
                        break;
                    }

                    f_index += n*int(pow(4,kmer-1-b));
                    r_index += n*int(pow(4,b));
                }
				kfv_f[x][f_index] += (1.0/kmertot);
                kfv_r[x][r_index] += (1.0/kmertot);
            }
        }
        
        //now report on the scores that were determined
        f_index = 0; //use f_index as an indexing variable here to save memory
        for(int i=0; i<4; i++)
            for(int j=0; j<4; j++)
                for(int k=0; k<4; k++)
                    for(int l=0; l<4; l++)
                    {
                        fprintf(out, "%c%c%c%c   %.5f    %.5f\n", sym[i], sym[j], sym[k], sym[l], kfv_f[x][f_index], kfv_r[x][f_index]);
                        f_index++;
                    }
        //finished with this sequence, move to the next
        fprintf(out, "\n\n");
    }
    //all sequences finished!
    	
    //output file for R import
    fprintf(outR, "Name    Direction    ");
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                for(int l=0; l<4; l++)
                    fprintf(outR, "%c%c%c%c	   ", sym[i], sym[j], sym[k], sym[l]);
	fprintf(outR, "\n");
	
	for(int x=0; x<DN; x++){
		fprintf(outR, "%s    1   ", sequences[x].name);
		for(int y=0; y<kfv_size;y++)
            fprintf(outR, "%.5f ", kfv_f[x][y]);
        fprintf(outR, "\n");
        
        fprintf(outR, "%s    2   ", sequences[x].name);
		for(int y=0; y<kfv_size;y++)
            fprintf(outR, "%.5f ", kfv_r[x][y]);
        fprintf(outR, "\n");
    }
}
