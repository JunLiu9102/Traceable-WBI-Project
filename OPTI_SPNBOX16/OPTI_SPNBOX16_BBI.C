/**************    optimized black-box implementation of SPNBOX-16    ************/
/*******            2019.12                          ***************/
/**********     junl1212@163.com            ***********************/
/*********    encryption only (decryption is not covenient because AES16 decryption is not covenient by using AES-NI   ************/

#include <stdlib.h>
#include <sys/time.h>
#include <givaro/gfq.h>
#include <ctime>
#include <stdio.h>
#include<time.h>
#include<errno.h>
#include<vector>


#include <stdint.h>     //for int8_t
#include <string.h>     //for memcmp
#include <wmmintrin.h>  // intrinsics for AES-NI
#include <emmintrin.h>  // intrinsics for SSE2
#include <tmmintrin.h>
#include <immintrin.h>   //intrinsics for AVX2
#include <avx2intrin.h>

#define u8 uint8_t
#define u16 uint16_t   
#define u32 uint32_t
#define u64 uint64_t
#define roundaes16 32
#define round 10  
#define block 65536// 1MB data = 65536 blocks
#define parallel block/8     //deal with 8 blocks at one time
#define loops 100 //loops =100

#ifndef __AES_NI_H__
#define __AES_NI_H__

using namespace Givaro;
int modulus[] = { 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }; // x^16 + x^5 + x^3 + x + 1  (0x2b)
GFqDom<int32_t> GF216(2, 16, modulus);

static __m128i invSRindex =_mm_set_epi8 (0x03,0x06,0x09,0x0c,0x0f,0x02,0x05,0x08,0x0b,0x0e,0x01,0x04,0x07,0x0a,0x0d,0x00); // inverse shiftrows index in AES
static __m128i SRindex    =_mm_set_epi8 (0x0B,0x06,0x01,0x0c,0x07,0x02,0x0D,0x08,0x03,0x0e,0x09,0x04,0x0F,0x0a,0x05,0x00); // shiftrows indes in AES
static __m128i allzero    =_mm_set_epi8 (0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00); 
static __m128i rk[roundaes16 +1];   //round key array
static __m128i MSB8_m= _mm_set_epi16(0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000);//used in finite field arithmetic
static __m128i xmm[8][8]; //used in finite field arithmetic
static __m128i row[8];     //used in finite field arithmetic
static __m128i input[block];   // message and copy message

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

#endif
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

/**************       AES16  ***************************/

/***************    key generation of AES-16 :KDF  *********************/
void AES16KeySchedule(__m128i rk[roundaes16 +1])   
{
    for (u8 i=0;i<roundaes16 +1;i++)
    {
        //printf("%d-th round key is :",i);
		rk[i] = _mm_set_epi8 (0x00,0x00,i+1,i,0x00,0x00,i+1,i,0x00,0x00,i+1,i,0x00,0x00,i+1,i);
		//print128_num(rk[i]);
    }
}

/*******************   (inv)nonlinear layer of SPNBOX16  ****************************/
///////    parallel: 2 AESENC instructions can compute one block at one time    ///////////////////
void nonlinear(__m128i *m, __m128i rk[roundaes16 +1]) // 
{
	__m128i m1=_mm_set_epi8 (0x00,0x00,((u8 *)&(*m))[6],((u8 *)&(*m))[7],0x00,0x00,((u8 *)&(*m))[4],((u8 *)&(*m))[5],0x00,0x00,((u8 *)&(*m))[2],((u8 *)&(*m))[3],0x00,0x00,((u8 *)&(*m))[0],((u8 *)&(*m))[1]);
    m1 = _mm_xor_si128    (m1, rk[0]); 
    for (int i=1;i<=roundaes16;i++)
    {
        m1 = _mm_set_epi8 (0x52,0x52,((u8 *)&(m1))[13],((u8 *)&(m1))[12],0x52,0x52,((u8 *)&(m1))[9],((u8 *)&(m1))[8],0x52,0x52,((u8 *)&(m1))[5],((u8 *)&(m1))[4],0x52,0x52,((u8 *)&(m1))[1],((u8 *)&(m1))[0]);//we lose data here
		m1 = _mm_shuffle_epi8        (m1, invSRindex);    //  first do inverse shiftrows
		m1 = _mm_aesenc_si128        (m1, rk[i]);    //then do a full encryption (take advantage of _mm_aesenc_si128 instruction!!)
    }   

    __m128i m2=_mm_set_epi8 (0x00,0x00,((u8 *)&(*m))[14],((u8 *)&(*m))[15],0x00,0x00,((u8 *)&(*m))[12],((u8 *)&(*m))[13],0x00,0x00,((u8 *)&(*m))[10],((u8 *)&(*m))[11],0x00,0x00,((u8 *)&(*m))[8],((u8 *)&(*m))[9]);
    m2 = _mm_xor_si128    (m2, rk[0]); 
    for (int i=1;i<=roundaes16;i++)
    {
        m2 = _mm_set_epi8 (0x52,0x52,((u8 *)&(m2))[13],((u8 *)&(m2))[12],0x52,0x52,((u8 *)&(m2))[9],((u8 *)&(m2))[8],0x52,0x52,((u8 *)&(m2))[5],((u8 *)&(m2))[4],0x52,0x52,((u8 *)&(m2))[1],((u8 *)&(m2))[0]);//we lose data here
		m2 = _mm_shuffle_epi8        (m2, invSRindex);    //  first do inverse shiftrows			
		m2 = _mm_aesenc_si128        (m2, rk[i]);    //then do a full encryption (take advantage of _mm_aesenc_si128 instruction!!)
    }   
	*m=_mm_set_epi16(
		(((u8 *)&(m2))[12] << 8)^(((u8 *)&(m2))[13]),(((u8 *)&(m2))[8] << 8)^(((u8 *)&(m2))[9]),(((u8 *)&(m2))[4] << 8)^(((u8 *)&(m2))[5]),(((u8 *)&(m2))[0] << 8)^(((u8 *)&(m2))[1]),
	    (((u8 *)&(m1))[12] << 8)^(((u8 *)&(m1))[13]),(((u8 *)&(m1))[8] << 8)^(((u8 *)&(m1))[9]),(((u8 *)&(m1))[4] << 8)^(((u8 *)&(m1))[5]),(((u8 *)&(m1))[0] << 8)^(((u8 *)&(m1))[1])
		);
}

/*******************   encryption of SPNBOX16  ****************************/
void encryptionblack(__m128i input[block], __m128i rk[roundaes16 +1])
{
	
    for (int i = 0; i < round; i++)
    {
        for (int j = 0; j < parallel; j++)   // modify here: parallel compute 8 blocks one time
        {
   /******************            nonlinear layer          **************************************/
            nonlinear(&input[j*8+0],rk);nonlinear(&input[j*8+1],rk);nonlinear(&input[j*8+2],rk);nonlinear(&input[j*8+3],rk);nonlinear(&input[j*8+4],rk);nonlinear(&input[j*8+5],rk);nonlinear(&input[j*8+6],rk);nonlinear(&input[j*8+7],rk);
        
   /*****************   linear layer :  collect several blocks in one register  (for parallel)  and  xor corresponding places in every register      *********/ 
            for (int a=0;a<8;a++)
            {
				row[a]=allzero; 
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
            {input[j*8+t] ^= roundconstant;}
        }
    }
}

/*********************************************************************************************

BELOW are some useful functions!

************************************************************************************************/
double Average(double list[], int lenlist)
{
	double ave, sum = 0;
	for (int i = 0; i < lenlist; i++) {
		sum += list[i];
	}
	ave = sum / lenlist;
	return ave;
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

int main(int argc, char** argv)
{

	srand(time(0));
	double elapsed_secs[loops];
	AES16KeySchedule(rk);

	for (int k = 0; k < loops; k++)
	{
		printf("********* loop %d results below  ************\n", k);

		generatemessage(input);
		printf("PLAINTEXT is: \n");
		printmessage(input);
	    printf("\n");

		clock_t begin = clock();
		encryptionblack(input,rk); //encrypt message
		clock_t end = clock();
		elapsed_secs[k] = double(end - begin) / CLOCKS_PER_SEC;
        
		printf("\n");
		printf("CIPHERTEXT is: \n");
		printmessage(input);
		printf("\n");
	}
	double avetime = Average(elapsed_secs, loops);
	printf("average time is %f s\n", avetime);

	return 0;
}







