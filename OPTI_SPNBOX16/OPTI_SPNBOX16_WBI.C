/**************    optimized white-box implementation of SPNBOX-16   ************/
/*******            2019.12                           ***************/
/**********     junl1212@163.com              ***********************/

#include <stdlib.h>
#include <sys/time.h>
#include <givaro/gfq.h>
#include <ctime>
#include <stdio.h>
#include <sbox.h>
#include <invsbox.h>

#define round 10
#define block 65536 // 1MB data, block 65536
#define parallel block/8    //deal with 8 blocks at one time
#define loops 1 //loops 100
#define u8 uint8_t
#define u16 uint16_t

#include <stdint.h>     //for int8_t
#include <string.h>     //for memcmp
#include <wmmintrin.h>  // intrinsics for AES-NI
#include <emmintrin.h>  // intrinsics for SSE2
#include <tmmintrin.h>
#include <immintrin.h>   //intrinsics for AVX2
#include <avx2intrin.h>

static __m128i allzero =_mm_set_epi16 (0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000);
static __m128i MSB8_m= allzero;//used in finite field arithmetic
static __m128i xmm[8][8]; //used in finite field arithmetic
static __m128i row[8];     //used in finite field arithmetic
static __m128i input[block],temp[block];   // message and copy message

using namespace Givaro;
int modulus[] = { 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }; // x^16 + x^5 + x^3 + x + 1  (0x2b)
GFqDom<int32_t> GF216(2, 16, modulus);

/*__m128i xtime(__m128i x)   //   core function: 0x02 * parameter
{// whether MSB of x is 1 or 0? if it's 1, then compare is 0xffff, otherwise it's 0x0000
    __m128i compare = _mm_cmpgt_epi16(MSB8_m,x);//MSB8_m   >   x ? 0xffff : 0x0000
    __m128i count =_mm_set_epi16(0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0001);//count
    x=_mm_sll_epi16(x,count);  //left shift by 1
    __m128i poly=_mm_set_epi16(0x002b,0x002b,0x002b,0x002b,0x002b,0x002b,0x002b,0x002b);//reduction poly
    compare=_mm_and_si128 (poly,compare);//determine whether XOR poly, if compare is 0xffff, with XOR ; otherwise no XOR 
    x=_mm_xor_si128(compare,x);
    return x;
}
__m128i mul(u16 operand, __m128i x)
{
	if (operand==0x01)
		return  x;                                               				// 0x01 * parameter
	if (operand==0x03)
	   return _mm_xor_si128(xtime(x),x); 										 // 0x03 * parameter
	if (operand==0x04)
	   return xtime(xtime(x)); 													// 0x04 * parameter
	if (operand==0x05)
	   return _mm_xor_si128(xtime(xtime(x)),x); 								// 0x05 * parameter
	if (operand==0x06)
	   return _mm_xor_si128(xtime(xtime(x)),xtime(x));							 // 0x06 * parameter
	if (operand==0x08)
	   return xtime(xtime(xtime(x)));											// 0x08 * parameter
	if (operand==0x0b)
	   return _mm_xor_si128(_mm_xor_si128(xtime(xtime(xtime(x))),xtime(x)),x); // 0x0b * parameter
	if (operand==0x07)
	   return _mm_xor_si128(_mm_xor_si128(xtime(xtime(x)),xtime(x)),x);        // 0x07 * parameter
}*/
__m128i mul(u16 operand, __m128i x)   // FINITE FIELD MULTIPLICATION
{
		for (u16 k = 0; k < 8; k++)
		{
			GFqDom<int32_t>::Element a, b, c;
			GF216.init(a, operand);   // initialize
			GF216.init(b, ((u16 *)&x)[k]);
			GF216.mul(c, a, b);   // field multiplication
			int32_t c_int;
		    GF216.convert(c_int, c);
			((u16 *)&x)[k]=c_int;
		}
		return x;		
}

/*********************   print __m128i     ****************************/
void print128_num(__m128i var)
{
	  u16 *val=(u16*) &var;
      for (int i=0;i<8;i++)
      {
            printf("%04X",val[i]);
      }
      //printf("\n");
}

/****************     MDS matrix comes from block cipher khazad    ***********************/
u16 PM[8][8] = { //khazad matrix = had(1,3,4,5,6,8,b,7)
	{0x01,0x03,0x04,0x05,0x06,0x08,0x0b,0x07},
	{0x03,0x01,0x05,0x04,0x08,0x06,0x07,0x0b},
	{0x04,0x05,0x01,0x03,0x0b,0x07,0x06,0x08},
	{0x05,0x04,0x03,0x01,0x07,0x0b,0x08,0x06},
	{0x06,0x08,0x0b,0x07,0x01,0x03,0x04,0x05},  
	{0x08,0x06,0x07,0x0b,0x03,0x01,0x05,0x04},
	{0x0b,0x07,0x06,0x08,0x04,0x05,0x01,0x03},
	{0x07,0x0b,0x08,0x06,0x05,0x04,0x03,0x01}
}; /*************************************************************************************/


/*******************   (inv)nonlinear layer of SPNBOX-16 :LUT  ****************************/
void nonlinear(__m128i *m)
{
	*m=_mm_set_epi16(sbox[((u16 *)&(*m))[7]],sbox[((u16 *)&(*m))[6]],sbox[((u16 *)&(*m))[5]],sbox[((u16 *)&(*m))[4]],sbox[((u16 *)&(*m))[3]],sbox[((u16 *)&(*m))[2]],sbox[((u16 *)&(*m))[1]],sbox[((u16 *)&(*m))[0]]);
	
}
void invnonlinear(__m128i *m)
{
	*m=_mm_set_epi16(invsbox[((u16 *)&(*m))[7]],invsbox[((u16 *)&(*m))[6]],invsbox[((u16 *)&(*m))[5]],invsbox[((u16 *)&(*m))[4]],invsbox[((u16 *)&(*m))[3]],invsbox[((u16 *)&(*m))[2]],invsbox[((u16 *)&(*m))[1]],invsbox[((u16 *)&(*m))[0]]);
}
/*******************   encryption of SPNBOX-16  ****************************/
void encryptionwhite(__m128i input[block])
{
	for (int i = 0; i < round; i++)
	{
		for (int j = 0; j < parallel; j++)
		{
			/*****************     nonlinear layer        ******************************/
			nonlinear(&input[j*8+0]);nonlinear(&input[j*8+1]);nonlinear(&input[j*8+2]);nonlinear(&input[j*8+3]);nonlinear(&input[j*8+4]);nonlinear(&input[j*8+5]);nonlinear(&input[j*8+6]);nonlinear(&input[j*8+7]);
            
			 /*****************   linear layer :  collect several blocks in one register  (for parallel)  and  xor corresponding places in every register      *********/
            for (int a=0;a<8;a++) row[a]=allzero;   
            for (int a=0;a<8;a++)
            {
                for (int b=0;b<8;b++)
                {
                            __m128i temp=_mm_set_epi16(
								((u16*)&input[j*8+7])[b],
								((u16*)&input[j*8+6])[b],
								((u16*)&input[j*8+5])[b],
								((u16*)&input[j*8+4])[b],
								((u16*)&input[j*8+3])[b],
								((u16*)&input[j*8+2])[b],
								((u16*)&input[j*8+1])[b],
								((u16*)&input[j*8+0])[b]
								);
                            xmm[a][b]=mul(PM[a][b],temp);
                            row[a]=_mm_xor_si128(row[a],xmm[a][b]);//row[a]=_mm_xor_si128(xmm[a][7],_mm_xor_si128(xmm[a][6],_mm_xor_si128(xmm[a][5],_mm_xor_si128(xmm[a][4],_mm_xor_si128(xmm[a][3],_mm_xor_si128(xmm[a][2],_mm_xor_si128(xmm[a][1],xmm[a][0])))))));
                }
    
            }
	        for (int t=0;t<8;t++)
	        {
				input[j*8+t]=_mm_set_epi16(((u16*)&(row[7]))[t],((u16*)&(row[6]))[t],((u16*)&(row[5]))[t],((u16*)&(row[4]))[t],((u16*)&(row[3]))[t],((u16*)&(row[2]))[t],((u16*)&(row[1]))[t],((u16*)&(row[0]))[t]);
	        }
            
   /******************    affine layer *********************/
            __m128i roundconstant=_mm_set_epi16((u16)(8*i+8),(u16)(8*i+7),(u16)(8*i+6),(u16)(8*i+5),(u16)(8*i+4),(u16)(8*i+3),(u16)(8*i+2),(u16)(8*i+1));
            for (int t = 0; t < 8; t++)  //affine layer
            {
               input[j*8+t] ^= roundconstant;
            }
		}
	}
}
/*******************   decryption of SPNBOX-16  ****************************/
void decryptionwhite(__m128i input[block])
{	
	for (int i = round - 1; i >= 0; i--)
	{
		for (int j = 0; j < parallel; j++)
		{
			/******************     inverse    affine layer *********************/
            __m128i roundconstant=_mm_set_epi16((u16)(8*i+8),(u16)(8*i+7),(u16)(8*i+6),(u16)(8*i+5),(u16)(8*i+4),(u16)(8*i+3),(u16)(8*i+2),(u16)(8*i+1));
            for (int t = 0; t < 8; t++)  //affine layer
            {
               input[j*8+t] ^= roundconstant;
            }
			/********* inverse  linear layer (the same as linear layer):  collect several blocks in one register  (for parallel)  and  xor corresponding places in every register      *********/
            for (int a=0;a<8;a++) row[a]=allzero;   
            for (int a=0;a<8;a++)
            {
                for (int b=0;b<8;b++)
                {
                            __m128i temp=_mm_set_epi16(
								((u16*)&input[j*8+7])[b],
								((u16*)&input[j*8+6])[b],
								((u16*)&input[j*8+5])[b],
								((u16*)&input[j*8+4])[b],
								((u16*)&input[j*8+3])[b],
								((u16*)&input[j*8+2])[b],
								((u16*)&input[j*8+1])[b],
								((u16*)&input[j*8+0])[b]
								);
                            xmm[a][b]=mul(PM[a][b],temp);
                            row[a]=_mm_xor_si128(row[a],xmm[a][b]);//row[a]=_mm_xor_si128(xmm[a][7],_mm_xor_si128(xmm[a][6],_mm_xor_si128(xmm[a][5],_mm_xor_si128(xmm[a][4],_mm_xor_si128(xmm[a][3],_mm_xor_si128(xmm[a][2],_mm_xor_si128(xmm[a][1],xmm[a][0])))))));
                }
    
            }
	        for (int t=0;t<8;t++)
	        {
				input[j*8+t]=_mm_set_epi16(((u16*)&(row[7]))[t],((u16*)&(row[6]))[t],((u16*)&(row[5]))[t],((u16*)&(row[4]))[t],((u16*)&(row[3]))[t],((u16*)&(row[2]))[t],((u16*)&(row[1]))[t],((u16*)&(row[0]))[t]);
	        }
			 /******************    inverse        nonlinear layer          **************************************/
            invnonlinear(&input[j*8+0]);invnonlinear(&input[j*8+1]);invnonlinear(&input[j*8+2]);invnonlinear(&input[j*8+3]);invnonlinear(&input[j*8+4]);invnonlinear(&input[j*8+5]);invnonlinear(&input[j*8+6]);invnonlinear(&input[j*8+7]);			
		}
	}	
}
/*********************************************************************************************

BELOW are some useful functions!

************************************************************************************************/

/***********************calculate average of list elements ***************************/
double Average(double list[], int lenlist)
{
	double ave, sum = 0;
	for (int i = 0; i < lenlist; i++) {
		sum += list[i];
	}
	ave = sum / lenlist;
	return ave;
}
/***********************  verify decryption   **********************************/
void verifydecryption(__m128i arr1[block], __m128i arr2[block])
{
	int i = memcmp(arr1, arr2, 16*block);
	if (i == 0)
		printf("VERIFY DECRYPTION CORRECT!");
	else
		printf("VERIFY DECRYPTION WRONG!");
}
void verifydecryption2(__m128i arr1[block], __m128i arr2[block])   // check how many blocks are wrong
{
	int counter=0;
	for (int i=0; i < block; i++)
	{
		if (((u16*)&(arr1[i]))[0]!=((u16*)&(arr2[i]))[0] || ((u16*)&(arr1[i]))[1]!=((u16*)&(arr2[i]))[1]  || ((u16*)&(arr1[i]))[2]!=((u16*)&(arr2[i]))[2] || ((u16*)&(arr1[i]))[3]!=((u16*)&(arr2[i]))[3] || ((u16*)&(arr1[i]))[4]!=((u16*)&(arr2[i]))[4] || ((u16*)&(arr1[i]))[5]!=((u16*)&(arr2[i]))[5] || ((u16*)&(arr1[i]))[6]!=((u16*)&(arr2[i]))[6] || ((u16*)&(arr1[i]))[7]!=((u16*)&(arr2[i]))[7])
		{
		    counter++;
			//printf("the %d-th block is wrong!!!\n",i);
			//printf("wrong blocks of arr1 are %04X,%04X,%04X,%04X,%04X,%04X,%04X,%04X",arr1[i][0],arr1[i][1],arr1[i][2],arr1[i][3],arr1[i][4],arr1[i][5],arr1[i][6],arr1[i][7]);
			//printf("\n");
			//printf("wrong blocks of arr2 are %04X,%04X,%04X,%04X,%04X,%04X,%04X,%04X",arr2[i][0],arr2[i][1],arr2[i][2],arr2[i][3],arr2[i][4],arr2[i][5],arr2[i][6],arr2[i][7]);
			//printf("\n");
		}
	}
	printf("there are %d wrong blocks in total",counter);
}
void verifydecryption3(__m128i arr1[block], __m128i arr2[block])// check how many bytes are wrong
{
	int counter=0;
	for (int i=0; i < block; i++)
	{
		for (int j=0;j<8;j++)
		{
			if (((u16*)&(arr1[i]))[j]!=((u16*)&(arr2[i]))[j])
			{
				counter++;
				printf("%d-th byte in %d-th block is wrong\n",j,i);
			}
		}
	}
	printf("there are %d wrong bytes in total",counter);
} 
void printmessage(__m128i arr[block])
{
	for (int i = 0; i < block; i++)
	{
		print128_num(arr[i]);
	}
		
}
/********************** generate random message (plaintext)  ************************************/
void generatemessage(__m128i arr1[block])
{
	for (int i = 0; i < block; i++)
	{
		arr1[i] = _mm_set_epi16(rand(),rand(),rand(),rand(),rand(),rand(),rand(),rand());
	}
}

/***************************************                     main function           ***************************************************************/
int main() 
{
	/*srand(time(0));
	double elapsed_secs[loops];


	for (int k = 0; k < loops; k++)
	{	
		printf("********* loop %d results below  ************\n", k);

		generatemessage(input);
		memcpy(temp, input, 16*block);

		printf("PLAINTEXT is: \n");
		printmessage(input);		
	
		clock_t begin = clock(); // be careful about the begin time 
		encryptionwhite(input);//encrypt message
		clock_t end = clock();
		elapsed_secs[k] = double(end - begin) / CLOCKS_PER_SEC;

		printf("\n");
		printf("CIPHERTEXT is: \n");
		printmessage(input);

		decryptionwhite(input);//decrypt message

		printf("\n");
		printf("VERIFY PLAINTEXT is: \n");
		printmessage(input);

		printf("\n");
		verifydecryption(input,temp);	//verify decryption result
		printf("\n");
	}
	double avetime = Average(elapsed_secs,loops);
	printf("\n");
	printf("average time for %d loops encryption is %f s\n",loops, avetime);
	*/



/********  BELOW IS TRACEABLE TEST     *************/
    srand(time(0));
	double elapsed_secs[loops];


	for (int k = 0; k < loops; k++)
	{	
		printf("********* loop %d results below  ************\n", k);

		FILE *fp;  // read random ciphertext
        fp=fopen("pertu_ciphertext01.txt","r");
		rewind(fp);
		for (int i=0;i<65536;i++)
		{
		  //fscanf(fp,"%x,",&input[i][0]);fscanf(fp,"%x,",&input[i][1]);fscanf(fp,"%x,",&input[i][2]);fscanf(fp,"%x,",&input[i][3]);fscanf(fp,"%x,",&input[i][4]);fscanf(fp,"%x,",&input[i][5]);fscanf(fp,"%x,",&input[i][6]);fscanf(fp,"%x,",&input[i][7]);
		   fscanf(fp,"%x,",&(((u16*)&input[i])[0]));fscanf(fp,"%x,",&(((u16*)&input[i])[1]));fscanf(fp,"%x,",&(((u16*)&input[i])[2]));fscanf(fp,"%x,",&(((u16*)&input[i])[3]));fscanf(fp,"%x,",&(((u16*)&input[i])[4]));fscanf(fp,"%x,",&(((u16*)&input[i])[5]));fscanf(fp,"%x,",&(((u16*)&input[i])[6]));fscanf(fp,"%x,",&(((u16*)&input[i])[7]));
		}	  	
        fclose(fp);		
	
		clock_t begin = clock(); // be careful about the begin time 
	    decryptionwhite(input);//decrypt message
		clock_t end = clock();
		elapsed_secs[k] = double(end - begin) / CLOCKS_PER_SEC;

		fp=fopen("testtraitor_plaintext01.txt","w+"); // write new plaintext of the random ciphertext
        for (int i=0;i<block;i++) 
        {
           // fprintf(fp,"%04X,%04X,%04X,%04X,%04X,%04X,%04X,%04X,",
             //   input[i][0], input[i][1], input[i][2], input[i][3], input[i][4], input[i][5], input[i][6], input[i][7]);
				fprintf(fp,"%04X,%04X,%04X,%04X,%04X,%04X,%04X,%04X,",((u16*)&input[i])[0], ((u16*)&input[i])[1], ((u16*)&input[i])[2],((u16*)&input[i])[3],((u16*)&input[i])[4],((u16*)&input[i])[5],((u16*)&input[i])[6],((u16*)&input[i])[7]);
				
        }   
        fclose(fp);	

		printf("\n\n\n\n\n\n\n\n\n");
		fp=fopen("newpertu_plaintext01.txt","r");// read the original perturbation plaintext
        for (int i=0;i<block;i++) 
        {
            //fscanf(fp,"%x,",&temp[i][0]);fscanf(fp,"%x,",&temp[i][1]);fscanf(fp,"%x,",&temp[i][2]);fscanf(fp,"%x,",&temp[i][3]);fscanf(fp,"%x,",&temp[i][4]);fscanf(fp,"%x,",&temp[i][5]);fscanf(fp,"%x,",&temp[i][6]);fscanf(fp,"%x,",&temp[i][7]);	
			fscanf(fp,"%x,",&(((u16*)&temp[i])[0]));fscanf(fp,"%x,",&(((u16*)&temp[i])[1]));fscanf(fp,"%x,",&(((u16*)&temp[i])[2]));fscanf(fp,"%x,",&(((u16*)&temp[i])[3]));fscanf(fp,"%x,",&(((u16*)&temp[i])[4]));fscanf(fp,"%x,",&(((u16*)&temp[i])[5]));fscanf(fp,"%x,",&(((u16*)&temp[i])[6]));fscanf(fp,"%x,",&(((u16*)&temp[i])[7]));
        }   
        fclose(fp);

		printf("\n\n\n\n\n\n\n\n\n");
		//verifydecryption(input,temp);
		verifydecryption2(input,temp);
        //verifydecryption3(input,temp);

		

	
	}
	double avetime = Average(elapsed_secs,loops);
	printf("\n");
	printf("average time for %d loops encryption is %f s\n",loops, avetime);

	return 0;
	
}
