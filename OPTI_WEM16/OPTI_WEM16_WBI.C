/**************    optimized white-box implementation of WEM-16   ************/
/*******            2019.12                           ***************/
/**********     junl1212@163.com              ***********************/

#include <stdlib.h>
#include<sys/time.h>
#include <ctime>
#include <stdio.h>
#include <string.h>
//#include <iostream>
//#include <new> 
#include <sbox.h>
#include <newinvsbox10.h>



#define round 12
#define roundAES 5
#define block 65536// 1MB data, block 65536, 1GB data block 65536*1024
#define loops 1 //loops 100
#define u8 uint8_t
#define u16 uint16_t
#define u32 uint32_t

#ifndef __AES_NI_H__
#define __AES_NI_H__

#include <stdint.h>     //for int8_t
#include <string.h>     //for memcmp
#include <wmmintrin.h>  //intrinsics for AES-NI
#include <stdio.h>
#define u8 uint8_t
#define u16 uint16_t

//macros
#define DO_ENC_BLOCK(m,k) \
    do{\
        m = _mm_xor_si128           (m, k[ 0]); \
        m = _mm_aesenc_si128    (m, k[1]); \
        m = _mm_aesenc_si128    (m, k[2]); \
        m = _mm_aesenc_si128    (m, k[3]); \
        m = _mm_aesenc_si128    (m, k[4]); \
        m = _mm_aesenclast_si128    (m, k[5]);\ 
    }while(0)

#define DO_DEC_BLOCK(m,k) \
    do{\
        m = _mm_xor_si128           (m, k[5+0]); \
        m = _mm_aesdec_si128    (m, k[6]); \
        m = _mm_aesdec_si128    (m, k[7]); \
        m = _mm_aesdec_si128    (m, k[8]); \
        m = _mm_aesdec_si128    (m, k[9]); \
        m = _mm_aesdeclast_si128    (m, k[0]);\
    }while(0)

#define AES_128_key_exp(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))
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
/****************************    AES key schedule     ****************************************/
static __m128i aes_128_key_expansion(__m128i key, __m128i keygened)
{
    keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    return _mm_xor_si128(key, keygened);
}

//public API

static void aes128_load_key_enc_only(u8 * masterkey, __m128i *roundkey)
{
    roundkey[0] = _mm_loadu_si128((const __m128i*) masterkey);
    roundkey[1]  = AES_128_key_exp(roundkey[0], 0x01);
    roundkey[2]  = AES_128_key_exp(roundkey[1], 0x02);
    roundkey[3]  = AES_128_key_exp(roundkey[2], 0x04);
    roundkey[4]  = AES_128_key_exp(roundkey[3], 0x08);
    roundkey[5]  = AES_128_key_exp(roundkey[4], 0x10);

    /*for (int i=0;i<=5;i++)
    {    
        printf("%d-th round key is: ",i);
        print128_num(roundkey[i]);
    }*/
}

static void aes128_load_key(u8 * masterkey, __m128i *roundkey)
{
    aes128_load_key_enc_only(masterkey, roundkey);
    roundkey[6] = _mm_aesimc_si128(roundkey[4]);
    roundkey[7] = _mm_aesimc_si128(roundkey[3]);
    roundkey[8] = _mm_aesimc_si128(roundkey[2]);
    roundkey[9] = _mm_aesimc_si128(roundkey[1]);   

    /*for (int i=6;i<=9;i++)
    {    
        printf("%d-th round key is: ",i);
        print128_num(roundkey[i]);
    }*/
}

static void aes128_enc(__m128i *roundkey, u8 *plainText)
{
    __m128i m = _mm_loadu_si128((__m128i *) plainText);

    DO_ENC_BLOCK(m,roundkey);

    _mm_storeu_si128((__m128i *) plainText, m);
}

static void aes128_dec(__m128i *roundkey, u8 *cipherText)
{
    __m128i m = _mm_loadu_si128((__m128i *) cipherText);
 
    DO_DEC_BLOCK(m,roundkey);

    _mm_storeu_si128((__m128i *) cipherText, m);
}
#endif


/*******************   (inv)nonlinear layer of WEM-16 :LUT  ****************************/
void nonlinear(u16 arr[8])
{
	for (u16 i = 0; i < 8; i++) {
		arr[i]= sbox[arr[i]];
	}
}
void invnonlinear(u16 arr[8])
{
	for (u16 i = 0; i < 8; i++) {	
		arr[i] = invsbox[arr[i]];	
	}
}
void transfer168(u16 arr1[8], u8 arr2[16])
{
	for (int i = 0; i < 16; i=i+2)
	{
		arr2[i] = (u8)((arr1[i / 2]) >> 8);
		arr2[i + 1] = (u8)(arr1[i / 2]);
	}

}
void transfer816(u8 arr1[16], u16 arr2[8])
{
	for (int i = 0; i < 8; i ++)
	{
		arr2[i] = ((arr1[2 * i]) << 8) ^ arr1[2 * i + 1];
	}
}
/*******************   (inv)linear layer of WEM-16:5-round AES ****************************/
void linear(u16 arr[8], __m128i rk[10])
{
	u8 arr2[16];
	transfer168(arr, arr2);
    aes128_enc(rk,arr2);
	transfer816(arr2,arr);
}
void invlinear(u16 arr[8], __m128i rk[10])
{
	u8 arr2[16];
	transfer168(arr, arr2);
    aes128_dec(rk,arr2);
	transfer816(arr2, arr);
}
/*******************   encryption of WEM-16  ****************************/
void encryptionwhite(u16 input[block][8], __m128i rk[10])
{
	for (int i = 0; i < round; i++)
	{
		for (int j = 0; j < block; j++)
		{
			nonlinear(input[j]);  // nonlinear layer
			linear(input[j],rk);    //linear layer	
		}
	}
	for (int j = 0; j < block; j++)
	{
		nonlinear(input[j]);  // last nonlinear layer
	}
}
/*******************   decryption of WEM-16  ****************************/
void decryptionwhite(u16 input[block][8], __m128i rk[10])
{	
	for (int j = 0; j < block; j++)
	{
		invnonlinear(input[j]);  // nonlinear layer
	}
	for (int i = round - 1; i >= 0; i--)
	{
		for (int j = 0; j < block; j++)
		{
			invlinear(input[j],rk);    //invlinear layer 
			invnonlinear(input[j]);  // nonlinear layer
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
void printmessage(u16 arr[block][8])
{
	for (int i=0;i<block;i++) 
	   printf("%04X%04X%04X%04X%04X%04X%04X%04X",
		arr[i][0], arr[i][1], arr[i][2], arr[i][3], arr[i][4], arr[i][5], arr[i][6], arr[i][7]);
}
/***********************  verify decryption   **********************************/
void verifydecryption(u16 arr1[block][8],u16 arr2[block][8])
{
	int i = memcmp(arr1, arr2, 16*block);
	if (i == 0)
		printf("VERIFY DECRYPTION CORRECT!");
	else
		printf("VERIFY DECRYPTION WRONG!");
}
void verifydecryption2(u16 arr1[block][8],u16 arr2[block][8])   // check how many blocks are wrong
{
	int counter=0;
	for (int i=0; i < block; i++)
	{
		if (arr1[i][0]!=arr2[i][0] || arr1[i][1]!=arr2[i][1]  || arr1[i][2]!=arr2[i][2] || arr1[i][3]!=arr2[i][3] || arr1[i][4]!=arr2[i][4] || arr1[i][5]!=arr2[i][5] || arr1[i][6]!=arr2[i][6] || arr1[i][7]!=arr2[i][7])
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
void verifydecryption3(u16 arr1[block][8],u16 arr2[block][8])// check how many bytes are wrong
{
	int counter=0;
	for (int i=0; i < block; i++)
	{
		for (int j=0;j<8;j++)
		{
			if (arr1[i][j]!=arr2[i][j])
			{
				counter++;
				printf("%d-th byte in %d-th block is wrong\n",j,i);
			}
		}
	}
	printf("there are %d wrong bytes in total",counter);
} 
/********************** generate random message (plaintext)  ************************************/
void generatemessage(u16 arr1[block][8])
{
	for (int i = 0; i < block; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			arr1[i][j] = rand();
		}
	}
}

int main(int argc, char** argv) 
{

	
	srand(time(0));
	double elapsed_secs[loops];
	u16 input[block][8];// = new u16[65536*1024];
	u16 temp[block][8];// = new u16[block][8];
	u16 randomcipher[block][8];

    u8 masterkey[] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
    __m128i roundkey[10];
    printf("round keys following:\n");
    aes128_load_key(masterkey, roundkey);
	
	for (int k = 0; k < loops; k++)
	{	
		printf("********* %d-th loop results below  ************\n", k);

		/* generatemessage(input);	
		FILE *fp;
        fp=fopen("random_ciphertext.txt","w+");
        for (int i=0;i<block;i++) 
        {
            fprintf(fp,"%04X,%04X,%04X,%04X,%04X,%04X,%04X,%04X,",
                input[i][0], input[i][1], input[i][2], input[i][3], input[i][4], input[i][5], input[i][6], input[i][7]);
        }  
	
        fclose(fp);		
	 
		clock_t begin = clock(); // be careful about the begin time 
		decryptionwhite(input,roundkey);//decrypt message
		clock_t end = clock();
		elapsed_secs[k] = double(end - begin) / CLOCKS_PER_SEC;
        fp=fopen("random_plaintext.txt","w+");
        for (int i=0;i<block;i++) 
        {
            fprintf(fp,"%04X,%04X,%04X,%04X,%04X,%04X,%04X,%04X,",
                input[i][0], input[i][1], input[i][2], input[i][3], input[i][4], input[i][5], input[i][6], input[i][7]);
        }   
        fclose(fp);	

		*/
	
		FILE *fp;  // read random ciphertext
        fp=fopen("random_ciphertext.txt","r");
		rewind(fp);
		for (int i=0;i<65536;i++)
		{
		  fscanf(fp,"%x,",&input[i][0]);fscanf(fp,"%x,",&input[i][1]);fscanf(fp,"%x,",&input[i][2]);fscanf(fp,"%x,",&input[i][3]);fscanf(fp,"%x,",&input[i][4]);fscanf(fp,"%x,",&input[i][5]);fscanf(fp,"%x,",&input[i][6]);fscanf(fp,"%x,",&input[i][7]);
		}	  	
        fclose(fp);	
		
		clock_t begin = clock(); // be careful about the begin time 
		decryptionwhite(input,roundkey);//decrypt message
		clock_t end = clock();
		elapsed_secs[k] = double(end - begin) / CLOCKS_PER_SEC;

        fp=fopen("random_plaintext10.txt","w+"); // write new plaintext of the random ciphertext
        for (int i=0;i<block;i++) 
        {
            fprintf(fp,"%04X,%04X,%04X,%04X,%04X,%04X,%04X,%04X,",
                input[i][0], input[i][1], input[i][2], input[i][3], input[i][4], input[i][5], input[i][6], input[i][7]);
        }   
        fclose(fp);	

		printf("\n\n\n\n\n\n\n\n\n");
		fp=fopen("random_plaintext.txt","r");// read the original perturbation plaintext
        for (int i=0;i<block;i++) 
        {
            fscanf(fp,"%x,",&temp[i][0]);fscanf(fp,"%x,",&temp[i][1]);fscanf(fp,"%x,",&temp[i][2]);fscanf(fp,"%x,",&temp[i][3]);fscanf(fp,"%x,",&temp[i][4]);fscanf(fp,"%x,",&temp[i][5]);fscanf(fp,"%x,",&temp[i][6]);fscanf(fp,"%x,",&temp[i][7]);	
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
