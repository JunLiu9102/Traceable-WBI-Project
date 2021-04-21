/************************ optimized  S-box generation of SPNbox-16     ***********/
/*************************            2019.12    *************************/
/*************************   junl1212@163.com******************************/

#define u8 uint8_t
#define u16 uint16_t   //
#define u32 uint32_t
#define roundaes16 32

#ifndef __AES_NI_H__
#define __AES_NI_H__

#include <stdint.h>     //for int8_t
#include <string.h>     //for memcmp
#include <wmmintrin.h>  //for intrinsics for AES-NI
#include <emmintrin.h>  //for intrinsics for SSE
#include <tmmintrin.h>
#include <stdio.h>
#define u8 uint8_t
#define u16 uint16_t

static __m128i invSRindex=_mm_set_epi8 (0x03,0x06,0x09,0x0c,0x0f,0x02,0x05,0x08,0x0b,0x0e,0x01,0x04,0x07,0x0a,0x0d,0x00);
static __m128i SRindex     =_mm_set_epi8 (0x0B,0x06,0x01,0x0c,0x07,0x02,0x0D,0x08,0x03,0x0e,0x09,0x04,0x0F,0x0a,0x05,0x00);
static __m128i allzero        =_mm_set_epi8 (0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00);
static u8 cache[roundaes16+1][2]; // to store some data in encryption (the data lost in encryption)
static u16 sbox[65536];
static u16 invsbox[65536];
//macros


/*********************   print __m128i     ****************************/
void print128_num(__m128i var)
{
      u8 *val=(u8*) &var;
      for (int i=0;i<16;i++)
      {
            printf("0x%02X,",val[i]);
      }
      printf("\n");
}

//public API
/**************       AES16  ***************************/
static void AES16Encrypt(__m128i *m, __m128i *c, __m128i k[roundaes16 + 1])
{
            *c = _mm_xor_si128    (*m, k[0]); 
            for (int i=1;i<=roundaes16;i++)
            {
                  cache[i-1][0] = ((u8 *)&(*c))[3]; cache[i-1][1]=((u8 *)&(*c))[2];
                  *c = _mm_set_epi8 (0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x52,0x52,((u8 *)&(*c))[1],((u8 *)&(*c))[0]);
                  //we lose data here
                  *c = _mm_shuffle_epi8(*c, invSRindex);    //SB
                  *c = _mm_aesenclast_si128    (*c, allzero); 
                  *c = _mm_aesdeclast_si128    (*c, allzero);    //MC 
                  *c = _mm_aesenc_si128        (*c, allzero); 
                  *c = _mm_xor_si128           (*c, k[i]);  //ARK

            }   
}
static void AES16Decrypt(__m128i *c, __m128i *m, __m128i k[roundaes16 + 1])
{
             *m= _mm_xor_si128    (*c, k[roundaes16]); 
            for (int i= roundaes16-1;i>=0;i--)
            {
                  *m = _mm_aesenclast_si128  (*m,allzero);// invMC
                  *m = _mm_aesdec_si128      (*m,allzero);
                  *m = _mm_shuffle_epi8      (*m, SRindex);  //invSB
                  *m = _mm_aesdeclast_si128  (*m,allzero);
                  *m = _mm_set_epi8 (0x63,0x63,0x63,0x63,0x63,0x63,0x63,0x63,0x63,0x63,0x63,0x63,cache[i][0],cache[i][1],((u8 *)&(*m))[1],((u8 *)&(*m))[0]);
                  *m = _mm_xor_si128         (*m, k[i]);    //ARK
            }   
}
#endif

void generateroundkey(__m128i rk[roundaes16 +1])
{
      for (u8 i=0;i<roundaes16 +1;i++)
      {
            printf("%d-th round key is :",i);
            rk[i] = _mm_set_epi8 (0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,i+1,i);
            printf("0x%02X%02X\n",((u8 *)&rk[i])[0],((u8 *)&rk[i])[1]);
      }
}

void generatesbox(__m128i P[256][256], __m128i C[256][256], __m128i rk[roundaes16 + 1])
{
	for (u8 i=0x00; i<=0xff; i++)
      {
           for (u8 j=0x00;  j<=0xff;  j++)
           {
                 P[i][j] = _mm_set_epi8 (0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,j,i);
                 AES16Encrypt(&P[i][j],&C[i][j],rk);
                 //AES16Decrypt(&C[i][j],&P[i][j],rk);  //verify
                 sbox[256*i+j]= ((((u8 *)&C[i][j])[0])<<8) ^ (((u8 *)&C[i][j])[1]); 
                 printf("0x%02X%02X,",((u8 *)&C[i][j])[0],((u8 *)&C[i][j])[1]);       
                 if (j==0xff)
                    break;
           }
           if (i==0xff)
           break;  
      }     
}
void generateinvsbox(u16 arr1[65536], u16 arr2[65536])  //return index in sboxgen[65536],
/*why not AES16 decryption: because AES-NI is not appropriate for AES-16 decryption (even encryption is convenient),
AES16Decrypt defined above is only appropriate for AES16Encrypt defined above, but it can not be used alone! */
{
	u16 i, j;
	for ( i = 0x0000; i <= 0xffff; i++)
	{
		for ( j = 0; j <= 0xffff; j++)
		{
			if (arr1[j] == i)
				break;
		}
	      arr2[i] = j;
		printf("0x%04X,",arr2[i]);
		if (i == 0xffff)
			break;
	}
}

int main()
{
      __m128i rk[roundaes16 +1];  __m128i plain[256][256];   __m128i cipher[256][256]; 
      generateroundkey(rk);
      printf("sbox following:\n");
      generatesbox(plain, cipher, rk);
      printf("\n\n\n\n\n\n");
	
	printf("inverse sbox following:\n");
	generateinvsbox(sbox,invsbox);
	printf("\n");

	return 0;

}
