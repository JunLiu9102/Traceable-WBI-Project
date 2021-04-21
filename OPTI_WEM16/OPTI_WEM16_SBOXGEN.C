/************************optimized  S-box generation of WEM-16     ***********/
/*************************            2019.12    *************************/
/*************************   junl1212@163.com******************************/
/*********************  ***************************/

#ifndef __AES_NI_H__
#define __AES_NI_H__

#include <stdint.h>     //for int8_t
#include <string.h>     //for memcmp
#include <wmmintrin.h>  //for intrinsics for AES-NI
#include <stdio.h>
#define u8 uint8_t
#define u16 uint16_t


//compile using gcc and following arguments: -g;-o xxx;-Wall;-msse2;-msse;-march=native;-maes

//internal stuff

//macros
#define DO_ENC_BLOCK(m,k) \
    do{\
        m = _mm_xor_si128           (m, k[ 0]); \
        m = _mm_aesenc_si128    (m, k[ 1]); \
        m = _mm_aesenc_si128    (m, k[ 2]); \
        m = _mm_aesenc_si128    (m, k[ 3]); \
        m = _mm_aesenc_si128    (m, k[ 4]); \
        m = _mm_aesenc_si128    (m, k[ 5]); \
        m = _mm_aesenc_si128    (m, k[ 6]); \
        m = _mm_aesenc_si128    (m, k[ 7]); \
        m = _mm_aesenc_si128    (m, k[ 8]); \
        m = _mm_aesenc_si128    (m, k[ 9]); \
        m = _mm_aesenclast_si128(m, k[10]);\
    }while(0)

#define DO_DEC_BLOCK(m,k) \
    do{\
        m = _mm_xor_si128       (m, k[10+0]); \
        m = _mm_aesdec_si128    (m, k[10+1]); \
        m = _mm_aesdec_si128    (m, k[10+2]); \
        m = _mm_aesdec_si128    (m, k[10+3]); \
        m = _mm_aesdec_si128    (m, k[10+4]); \
        m = _mm_aesdec_si128    (m, k[10+5]); \
        m = _mm_aesdec_si128    (m, k[10+6]); \
        m = _mm_aesdec_si128    (m, k[10+7]); \
        m = _mm_aesdec_si128    (m, k[10+8]); \
        m = _mm_aesdec_si128    (m, k[10+9]); \
        m = _mm_aesdeclast_si128(m, k[0]);\
    }while(0)

#define AES_128_key_exp(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))
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
    roundkey[6]  = AES_128_key_exp(roundkey[5], 0x20);
    roundkey[7]  = AES_128_key_exp(roundkey[6], 0x40);
    roundkey[8]  = AES_128_key_exp(roundkey[7], 0x80);
    roundkey[9]  = AES_128_key_exp(roundkey[8], 0x1B);
    roundkey[10] = AES_128_key_exp(roundkey[9], 0x36);
    for (int i=0;i<=10;i++)
    {    
        printf("%d-th round key is: ",i);
        print128_num(roundkey[i]);
    }
}

static void aes128_load_key(u8 * masterkey, __m128i *roundkey)
{
    aes128_load_key_enc_only(masterkey, roundkey);

    // generate decryption keys in reverse order.
    // k[10] is shared by last encryption and first decryption rounds
    // k[0] is shared by first encryption round and last decryption round (and is the original user key)
    // For some implementation reasons, decryption key schedule is NOT the encryption key schedule in reverse order
    roundkey[11] = _mm_aesimc_si128(roundkey[9]);
    roundkey[12] = _mm_aesimc_si128(roundkey[8]);
    roundkey[13] = _mm_aesimc_si128(roundkey[7]);
    roundkey[14] = _mm_aesimc_si128(roundkey[6]);
    roundkey[15] = _mm_aesimc_si128(roundkey[5]);
    roundkey[16] = _mm_aesimc_si128(roundkey[4]);
    roundkey[17] = _mm_aesimc_si128(roundkey[3]);
    roundkey[18] = _mm_aesimc_si128(roundkey[2]);
    roundkey[19] = _mm_aesimc_si128(roundkey[1]);
}

static void aes128_enc(__m128i *roundkey, u8 *plainText, u8 *cipherText)
{
    __m128i m = _mm_loadu_si128((__m128i *) plainText);
    __m128i temp = m;

    DO_ENC_BLOCK(m,roundkey);
    m=_mm_xor_si128(m,temp);  // CTR

    _mm_storeu_si128((__m128i *) cipherText, m);
}

static void aes128_dec(__m128i *roundkey, u8 *cipherText, u8 *plainText)
{
    __m128i m = _mm_loadu_si128((__m128i *) cipherText);
    __m128i temp =_mm_loadu_si128((__m128i *) plainText);
   m=_mm_xor_si128(m,temp);  // CTR

    DO_DEC_BLOCK(m,roundkey);

    _mm_storeu_si128((__m128i *) plainText, m);
}
#endif



void generatesbox(u8 P[32][256][16], u8 C[32][256][16], __m128i rk[20])
{
	for (u8 j = 0x00; j < 0x20; j++)  // 
	{
		for (u8 i = 0x00; i <= 0xff; i++)
		{
			for (u8 k = 0; k < 14; k++)
			{
				P[j][i][k] = 0x00;
			}
			P[j][i][14] = j;
			P[j][i][15] = i;   // initilize counter and plaintext (lexicographical)
	
			aes128_enc(rk, P[j][i],C[j][i]);

			/*aes128_dec(rk, C[j][i],P[j][i]);  //verify encryption is correct
			for (u8 k = 0; k < 16; k++)  
			{				
				printf("%02X,", P[j][i][k]);
			}
			printf("\n");
			*/
			if (i == 0xff)   // be careful about dead loop
				break;
		}	
	}
}

void generateinvsbox(u16 arr1[65536], u16 arr2[65536])  //return index in sboxgen[65536]
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




void transfer(u8 arr1[32][256][16], u16 arr2[65536])
{
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j <256; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				arr2[2048 * i + 8 * j + k] = (arr1[i][j][2*k] << 8) ^ (arr1[i][j][2*k+1]);
			}
		}
	}
}

void FisherYates(u16 arr[65536],u16 arr2[65536])  //shuffle in arr2 (arr controls)
{
	for (int i = 0; i < 65536; i++)
	{
		arr2[i] = (u16)i;
	}
	for (int i = 65535; i >= 0; i--)
	{
	
		int j = ((int)(arr[i])) % (i+1);
		u16 temp = arr2[j];    //exchange arr2[j] and arr2[i]
 		arr2[j] = arr2[i];
		arr2[i] = temp;
	}
	for (int i = 0; i < 65536; i++)
	{
		printf("0x%04X,",arr2[i]);

	}
}



int main()
{

    u8 P[32][256][16]; //8192 AES-128 plaintexts
    u8 C[32][256][16]; //ciphertext  
    u16 FYarray[65536];//Fisher-Yates random sequence
    u16 sboxgen[65536],invsboxgen[65536];   //generated sbox

    u8 masterkey[] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
    __m128i roundkey[20];
    printf("round keys following:\n");
    aes128_load_key(masterkey, roundkey);

    printf("sbox following:\n");
	generatesbox(P, C, roundkey);
    transfer(C, FYarray);  //transfer sequence of random bits
	FisherYates(FYarray,sboxgen); //Fisher-Yates shuffle algorithm

	printf("\n\n\n");
	//checkpermutation(sboxgen); //verify whether the generated sbox is a PRP
    printf("inverse sbox following:\n");
	generateinvsbox(sboxgen, invsboxgen);
	//checkpermutation(invsboxgen);
	printf("\n\n\n");
	return 0;
}


