#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
#include <openssl/rand.h>
#include <openssl/sha.h>
#include <openssl/bn.h>
#include <math.h>

#define MAXTHREADS 1024		// TODO to be changed if more needed

const unsigned int ct_size = 512;	// Paillier cipher size
unsigned int m = 8000;		// size of bloom filter
unsigned int n = 1000;		// nbr of vehicules
unsigned int k = 4;			// nbr of hash functions
unsigned int p_size = 10;	// modulus bitsize
unsigned int p = 1021;		// modulus
unsigned int car_t = 4;	// nbr of parallel threads
unsigned int server_t = 4;	// nbr of parallel threads
unsigned int ct_count = 0;	// will be set later
unsigned int bytes_per_ct = 0;	// will be set later
unsigned int car_elems_per_thread = 0;	// will be =ct_count/car_t
unsigned int server_elems_per_thread = 0;	// will be =ct_count/server_t
unsigned char *sec, *enc_sec, *dec_sec;	// used by threads

typedef struct {
	BIGNUM *n;			/* public */
	BIGNUM *n2;			/* public */
	BIGNUM *phi;		/* private */
	BIGNUM *phi_inv;	/* private */
} Paillier_key;

typedef struct {
	BIGNUM *c;
} Paillier_cipher;

Paillier_key pk;

void init_paillier(const char *, const char *, Paillier_key *, BN_CTX *ctx);
void encrypt_paillier(Paillier_cipher *, BIGNUM *, Paillier_key *, BN_CTX *ctx);
void decrypt_paillier(BIGNUM *, Paillier_cipher *, Paillier_key *, BN_CTX *ctx);
void add_paillier(Paillier_cipher *, Paillier_cipher *, Paillier_cipher *, Paillier_key *, BN_CTX *ctx);
void mult_paillier(Paillier_cipher *, BIGNUM *, Paillier_cipher *, Paillier_key *, BN_CTX *ctx);
double print_time(struct timeval *start, struct timeval *end);
typedef enum {big, little} endianess_t;
static int bn2binpad(const BIGNUM *a, unsigned char *to, int tolen, endianess_t endianess);
void get_indices(size_t *idx, size_t k, size_t m);
void *sub_enc(void *vargp);
void *sub_dec(void *vargp);

void expand_16(unsigned char *in, uint16_t *out, size_t out_size, unsigned int unit_size) {
	unsigned int mask = pow(2, unit_size) - 1;
	for (size_t i=0; i<out_size; i++)
		out[i] = (uint16_t)((*(unsigned int*)(&in[i*unit_size/8])) >> i*unit_size%8) & mask;
}

void expand_64(unsigned char *in, unsigned long *out, size_t out_size, unsigned int unit_size) {
	unsigned long mask = pow(2, unit_size) - 1;
	for (size_t i=0; i<out_size; i++)
		out[i] = (unsigned long)((*(unsigned long*)(&in[i*unit_size/8])) >> i*unit_size%8) & mask;
}

void uint16_t2char(uint16_t *in, unsigned char *out, size_t in_size, unsigned int unit_size) {
	BIGNUM *bn_out = BN_new();
	BN_zero(bn_out);
	size_t i = in_size;
	do {
		i--;
		BN_lshift(bn_out, bn_out, unit_size);
		BN_add_word(bn_out, in[i]);
	} while (i>0);
	bn2binpad(bn_out, out, in_size*(unit_size)/8, little);
}

void ulong2char(unsigned long *in, unsigned char *out, size_t in_size, unsigned int unit_size) {
	BIGNUM *bn_out = BN_new();
	BN_zero(bn_out);
	size_t i = in_size;
	do {
		i--;
		BN_lshift(bn_out, bn_out, unit_size);
		BN_add_word(bn_out, in[i]);
	} while (i>0);
	bn2binpad(bn_out, out, in_size*(unit_size)/8, big);
}

int main(int argc, char **argv)
{
	struct timeval start, end, s, e;
	unsigned char *Q, *q;

	if (RAND_load_file("/dev/urandom", 1024) < 64) {
		printf("PNRG not seeded!\n");
		abort();
	}
	if (argc != 7) {
		printf("Usage: %s <m> <n> <k> <p_size> <car_t> <server_t>\n", argv[0]);
		exit(0);
	}
	m = atoi(argv[1]);
	if (m%64 != 0) {
		printf("For simulation purposes, m has to be multiple of 64!\n");
		abort();
	}
	n = atoi(argv[2]);
	unsigned int n_size = log2(n);
	if (n > pow(2, n_size))
		n_size++;
	k = atoi(argv[3]);
	p_size = atoi(argv[4]);
	if (p_size<7 || p_size>16) {
		printf("p_size should be in {7, 8, 9, 10, 11, 12, 13, 14, 15, 16}!\n");
		abort();
	}
	car_t = atoi(argv[5]);
	server_t = atoi(argv[6]);
	pthread_t tid[MAXTHREADS];
	unsigned int myid[MAXTHREADS];
	unsigned int tab_primes[] = {113, 251, 509, 1021, 2039, 4093, 8191, 16381, 32749, 65521};
	p = tab_primes[p_size-7];

	size_t i=0, j=0;
	BN_CTX *ctx = BN_CTX_new();
	init_paillier("pub.txt", "priv.txt", &pk, ctx);

	// Construct query:
	printf("Constructing query...\n");
	gettimeofday(&start,NULL);

	// 1.computing q:query and sec:secret
	Q = (unsigned char *) malloc(m*2*sizeof(unsigned char));
	uint16_t *Q_16 = (uint16_t *) Q;
	RAND_bytes(Q, m*2);
	for (i=0; i<m; i++)
		Q_16[i] %= p;
	sec = (unsigned char *) malloc(m*(p_size+n_size)/8*sizeof(unsigned char));
	uint16_t2char(Q_16, sec, m, p_size+n_size);

	// 2.changing k indices in q
	unsigned int temp=0;
	size_t idx[k];
	get_indices(idx, k, m);
	for (i=0; i<k; i++) {
		RAND_bytes(&temp, sizeof(unsigned int));
		temp = temp%(p-1) + 1;	// sould never add 0
		Q_16[idx[i]] = (Q_16[idx[i]]+temp)%p;
	}
	q = (unsigned char *) malloc(m*p_size/8*sizeof(unsigned char));
	uint16_t2char(Q_16, q, m, p_size);

	gettimeofday(&end,NULL);
	printf("%lf\n", print_time(&start, &end));

	// 3.encrypting sec into enc_sec
	printf("Encrypting sec...\n");
	gettimeofday(&start,NULL);
	unsigned int blocks_per_ct = 2047*8/m;
	bytes_per_ct = blocks_per_ct*m/64;
	ct_count = m*(p_size+n_size)/8/bytes_per_ct;
	if (m*(p_size+n_size)/8 % bytes_per_ct != 0)
		ct_count++;
	car_elems_per_thread = ct_count/car_t;
	//TODO check when the last ct contains less than blocks_per_ct blocks 
	printf("ct_count=%u\n", ct_count);
	enc_sec = (unsigned char *) malloc(ct_count*ct_size*sizeof(unsigned char));
	int v=0;
	for (v=0; v<car_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, sub_enc, &myid[v]);
	}
	for (v=0; v<car_t; v++)
		pthread_join(tid[v], NULL);

	gettimeofday(&end,NULL);
	printf("%lf\n", print_time(&start, &end));


	// process query:
	printf("Processing query...\n");
	// 0.preparing old query
	// for simulation we are adding 0 to q and enc(0) to enc(sec)
	unsigned int old_exp_q_16[m];
	for (i=0; i<m; i++)
		old_exp_q_16[i] = p;
	Paillier_cipher old_ctx[ct_count];
	BIGNUM *bn_0 = BN_new();
	BN_zero(bn_0);
	for (i=0; i<ct_count; i++) {
		old_ctx[i].c = BN_new();
		encrypt_paillier(&old_ctx[i], bn_0, &pk, ctx);
	}
	
	gettimeofday(&start,NULL);

	// 1.expand rec_q to exp_q_16
	uint16_t exp_q_16[m];
	expand_16(q, exp_q_16, m, p_size);

	// 2.perform plain addition of exp_q_16
	for (i=0; i<m; i++)
		exp_q_16[i] = (old_exp_q_16[i]+exp_q_16[i])%p;

	//3.Paillier addition of enc_sec
	Paillier_cipher ct;
	ct.c = BN_new();
	for (i=0; i<ct_count; i++) {
		BN_bin2bn(&enc_sec[i*ct_size], ct_size, ct.c);
		add_paillier(&ct, &ct, &old_ctx[i], &pk, ctx);
		bn2binpad(ct.c, &enc_sec[i*ct_size], ct_size, big);
	}

	// 4.shrink exp_q_16 to agg_q and ct_list to enc_sec
	unsigned char *agg_q = (unsigned char *) malloc(m*p_size/8*sizeof(unsigned char));
	uint16_t2char(exp_q_16, agg_q, m, p_size);

	for (i=0; i<ct_count; i++)

	gettimeofday(&end,NULL);
	printf("%lf\n", print_time(&start, &end));

	//Correctness check. Could be removed////////////////////////////////////////
	for (j=0; j<m; j++){
		if (Q_16[j] != exp_q_16[j]) {
			printf("ERROR: Q_16[%zu]=%hx != exp_q_16[%zu]=%hx\n", j, Q_16[j], j, exp_q_16[j]);
			exit(0);
		}
	}
	////////////////////////////////////////////////////////////////////////////

	// extract result:
	printf("Result decoding...\n");
	gettimeofday(&start,NULL);
	server_elems_per_thread = ct_count/server_t;
	// 1.decrypt all received ciphertexts (enc_sec) to dec_sec
	dec_sec = (unsigned char *) malloc(m*(p_size+n_size)/8*sizeof(unsigned char));
	for (v=0; v<server_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, sub_dec, &myid[v]);
	}
	for (v=0; v<server_t; v++)
		pthread_join(tid[v], NULL);

	// 2.expand dec_sec to exp_dec_sec_long
	unsigned long exp_dec_sec_long[m];
	expand_64(dec_sec, exp_dec_sec_long, m, p_size+n_size);
	for (i=0; i<m; i++)
		exp_dec_sec_long[i] %= p;

	// 3.expand agg_q to agg_q_16
	uint16_t agg_q_16[m];
	expand_16(agg_q, agg_q_16, m, p_size);

	// 4.computing bloom_filter_long = exp_dec_sec_long - agg_q_16
	unsigned long *bloom_filter_long = exp_dec_sec_long;
	for (i=0; i<m; i++)
		bloom_filter_long[i] = bloom_filter_long[i] - agg_q_16[i] != 0;

	gettimeofday(&end,NULL);
	printf("%lf\n", print_time(&start, &end));

	//TODO to delete/////////////////////////////////////////////////////////////
	printf("bloom_filter_long=");
	for (j=0; j<300; j++)
		printf("%lx", bloom_filter_long[j]);
	printf("\n");
	if (!memcmp(sec, dec_sec, m*(p_size+n_size)/8))
		printf("Decoding successful!! sec==dec_sec\n");
	else
		printf("Decoding failed :( sec!=dec_sec\n");
	if (!memcmp(agg_q, q, m*p_size/8))
		printf("Decoding successful!! agg_q==q\n");
	else
		printf("Decoding failed :( agg_q!=q\n");
	///////////////////////////////////////////////////////////////////////////////

	return 0;
}

void *sub_enc(void *vargp) {
	unsigned int myid = *((unsigned int *)vargp);
	unsigned int start = myid * car_elems_per_thread;
	unsigned int end = (myid != car_t - 1) ? start + car_elems_per_thread : ct_count;
	BN_CTX *ctx = BN_CTX_new();
	// printf("thread %u started: start=%d end=%d ct_count=%u car_elems_per_thread=%d bytes_per_ct=%d\n", myid, start, end, ct_count, car_elems_per_thread, bytes_per_ct);
	BIGNUM *pt = BN_new();
	Paillier_cipher ct;
	ct.c = BN_new();
	for (unsigned int i=start; i<end; i++) {
		BN_bin2bn(&sec[i*bytes_per_ct], bytes_per_ct, pt);
		encrypt_paillier(&ct, pt, &pk, ctx);
		bn2binpad(ct.c, &enc_sec[i*ct_size], ct_size, big);
	}
}

void *sub_dec(void *vargp) {
	unsigned int myid = *((unsigned int *)vargp);
	unsigned int start = myid * server_elems_per_thread;
	unsigned int end = (myid != server_t - 1) ? start + server_elems_per_thread : ct_count;
	BN_CTX *ctx = BN_CTX_new();
	BIGNUM *pt2 = BN_new();
	Paillier_cipher ct2;
	ct2.c = BN_new();
	for (unsigned int i=start; i<end; i++) {
		BN_bin2bn(&enc_sec[i*ct_size], ct_size, ct2.c);
		decrypt_paillier(pt2, &ct2, &pk, ctx);
		bn2binpad(pt2, &dec_sec[i*bytes_per_ct], bytes_per_ct, big);
	}
}

// void *sub_add(void *vargp) {
// 	unsigned int myid = *((unsigned int *)vargp);
// 	unsigned int start = myid * server_elems_per_thread;
// 	unsigned int end = (myid != server_t - 1) ? start + server_elems_per_thread : ct_count;
// 	BN_CTX *ctx = BN_CTX_new();
// 	Paillier_cipher ct;
// 	ct.c = BN_new();

// }

void get_indices(size_t *idx, size_t k, size_t m) {
	unsigned char *key = "secret_seed1";
	size_t key_length = strlen(key);
	unsigned char hash[SHA256_DIGEST_LENGTH];
	BIGNUM *bn_temp = BN_new();
	for(int l=0; l<k; l++) {
		SHA256(key, key_length, hash);
		BN_bin2bn(hash, SHA256_DIGEST_LENGTH, bn_temp);
		idx[l] = BN_mod_word(bn_temp, m);
		key = hash;
		key_length = SHA256_DIGEST_LENGTH;
	}
}

double print_time(struct timeval *start, struct timeval *end) {
	double usec;
	usec = (end->tv_sec*1000000 + end->tv_usec) - (start->tv_sec*1000000 + start->tv_usec);
	return usec/1000.0;
}

void init_paillier(const char *pub, const char *priv, Paillier_key *pk, BN_CTX *ctx) {
	FILE *fp;
	char str[2000], str2[2000];
	pk->n = BN_new();
	pk->n2 = BN_new();
	pk->phi = BN_new();
	pk->phi_inv = BN_new();
	fp = fopen(pub, "r");
	if (!fp) {
		printf("Error: Cannot open public key file %s\n", pub);
		exit(1);
	}
	fscanf(fp, "n\n%s\n", str);
	BN_dec2bn(&pk->n, str);
	fclose(fp);
	BN_sqr(pk->n2, pk->n, ctx);
	fp = fopen(priv, "r");
	if (!fp) {
		printf("Error: Cannot open private key file %s\n", priv);
		exit(1);
	}
	fscanf(fp, "phi\n%s\nphi_inv%s\n", str, str2);
	BN_dec2bn(&pk->phi, str);
	BN_dec2bn(&pk->phi_inv, str2);
	fclose(fp);
}

void encrypt_paillier(Paillier_cipher *pc, BIGNUM *m, Paillier_key *pk, BN_CTX *ctx) {
	BIGNUM *r, *aux, *aux2;
	r = BN_new();
	aux = BN_new();
	aux2 = BN_new();
	BN_rand_range(r, pk->n);
	BN_mod_exp(aux, r, pk->n, pk->n2, ctx);
	BN_mod_mul(aux2, m, pk->n, pk->n2, ctx);
	BN_add_word(aux2, 1);
	BN_mod_mul(pc->c, aux2, aux, pk->n2, ctx);
}

void decrypt_paillier(BIGNUM *m, Paillier_cipher *pc, Paillier_key *pk, BN_CTX *ctx) {
	BIGNUM *aux, *aux2;
	aux = BN_new();
	aux2 = BN_new();
	BN_mod_exp(aux, pc->c, pk->phi, pk->n2, ctx);
	BN_sub_word(aux, 1);
	BN_div(aux2, m, aux, pk->n, ctx);
	BN_mod_mul(m, aux2, pk->phi_inv, pk->n, ctx);
}

void add_paillier(Paillier_cipher *pc, Paillier_cipher *pc1, Paillier_cipher *pc2, Paillier_key *pk, BN_CTX *ctx) {
	BN_mod_mul(pc->c, pc1->c, pc2->c, pk->n2, ctx);
}


void mult_paillier(Paillier_cipher *pc, BIGNUM *m, Paillier_cipher *pc2, Paillier_key *pk, BN_CTX *ctx) {
	BN_mod_exp(pc->c, pc2->c, m, pk->n2, ctx);
}

static int bn2binpad(const BIGNUM *a, unsigned char *to, int tolen, endianess_t endianess) {
    int n;
    size_t i, lasti, j, atop, mask;
    BN_ULONG l;

    /*
     * In case |a| is fixed-top, BN_num_bytes can return bogus length,
     * but it's assumed that fixed-top inputs ought to be "nominated"
     * even for padded output, so it works out...
     */
    n = BN_num_bytes(a);
    if (tolen == -1) {
        tolen = n;
    } else if (tolen < n) {     /* uncommon/unlike case */
        BIGNUM temp = *a;

        bn_correct_top(&temp);
        n = BN_num_bytes(&temp);
        if (tolen < n)
            return -1;
    }

    /* Swipe through whole available data and don't give away padded zero. */
    atop = a->dmax * BN_BYTES;
    if (atop == 0) {
        OPENSSL_cleanse(to, tolen);
        return tolen;
    }

    lasti = atop - 1;
    atop = a->top * BN_BYTES;
    if (endianess == big)
        to += tolen; /* start from the end of the buffer */
    for (i = 0, j = 0; j < (size_t)tolen; j++) {
        unsigned char val;
        l = a->d[i / BN_BYTES];
        mask = 0 - ((j - atop) >> (8 * sizeof(i) - 1));
        val = (unsigned char)(l >> (8 * (i % BN_BYTES)) & mask);
        if (endianess == big)
            *--to = val;
        else
            *to++ = val;
        i += (i - lasti) >> (8 * sizeof(i) - 1); /* stay on last limb */
    }

    return tolen;
}