#include <immintrin.h>
#include "polynome.h"
#include "vect.h"

void vect_mod_add(Uint *res1, Uint *tab1, Uint *tab2)
{
    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i a = _mm256_loadu_si256((__m256i *)tab1);
    __m256i b = _mm256_loadu_si256((__m256i *)tab2);
    __m256i x = _mm256_add_epi32(a, b);
    __m256i result = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
    // __m256i y = _mm256_sub_epi32(x, p); // Première version avec les mask, plus lente
    // __m256i mask = _mm256_cmpgt_epi32(p, x);
    // __m256i result = _mm256_blendv_epi8(y, x, mask);
    _mm256_storeu_si256((__m256i *)res1, result);
}

void vect_mod_sub(Uint *res1, Uint *tab1, Uint *tab2)
{
    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i a = _mm256_loadu_si256((__m256i *)tab1);
    __m256i b = _mm256_loadu_si256((__m256i *)tab2);
    __m256i x = _mm256_sub_epi32(a, b);
    __m256i result = _mm256_min_epu32(x, _mm256_add_epi32(x, p));
    // __m256i y = _mm256_add_epi32(x, p); // Première version avec les mask, plus lente
    // __m256i mask = _mm256_cmpgt_epi32(p, y);
    // __m256i result = _mm256_blendv_epi8(x, y, mask);
    _mm256_storeu_si256((__m256i *)res1, result);
}

void vect_mod_mult(Uint *res, Uint *tab1, Uint *tab2)
{
    __m256i x = _mm256_loadu_si256((__m256i *)tab1);
    __m256i y = _mm256_loadu_si256((__m256i *)tab2);
    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i q = _mm256_set1_epi32(NB_Q);

    __m256i a_low = _mm256_mul_epu32(x, y);
    __m256i b_low = _mm256_srli_epi64(a_low, NB_S);
    __m256i c_low = _mm256_srli_epi64(_mm256_mul_epu32(b_low, q), NB_T);

    __m256i a_high = _mm256_mul_epu32(_mm256_srli_si256(x, 4), _mm256_srli_si256(y, 4));
    __m256i b_high = _mm256_srli_epi64(a_high, NB_S);
    __m256i c_high = _mm256_srli_epi64(_mm256_mul_epu32(b_high, q), NB_T);

    __m256i d_low = _mm256_sub_epi64(a_low, _mm256_mul_epu32(c_low, p));
    __m256i d_high = _mm256_sub_epi64(a_high, _mm256_mul_epu32(c_high, p));
    __m256i d = _mm256_or_si256(d_low, _mm256_slli_si256(d_high, 4));

    __m256i result = _mm256_min_epu32(d, _mm256_sub_epi32(d, p));
    _mm256_storeu_si256((__m256i *)res, result);
}

void vect_mod_add_sub_eval(Uint *res_add, Uint *res_sub, Uint *tab1, Uint *tab2)
{
    __m256i x, result;                               // Initialisation de variables
    __m256i p = _mm256_set1_epi32(NB_P);             // Chargement du tableau avec que des p
    __m256i a = _mm256_loadu_si256((__m256i *)tab1); // Chargement de 8 cases successives de tab1
    __m256i b = _mm256_loadu_si256((__m256i *)tab2); // Chargement de 8 cases successives de tab2

    // addition modulo p
    x = _mm256_add_epi32(a, b);                           // x = a + b
    result = _mm256_min_epu32(x, _mm256_sub_epi32(x, p)); // min_pos(x, x - p) = (a + b) % p
    _mm256_storeu_si256((__m256i *)res_add, result);      // Stockage du résultat dans res_add

    // soustraction modulo p
    x = _mm256_sub_epi32(a, b);                           // x = a - b
    result = _mm256_min_epu32(x, _mm256_add_epi32(x, p)); // min_pos(x, x + p) = (a - b) % p
    _mm256_storeu_si256((__m256i *)res_sub, result);      // Stockage du résultat dans res_sub
}

void vect_mod_mult_eval(Uint *res, Uint *tab1, Uint *tab2, Uint i, Uint pas)
{
    __m256i x = _mm256_loadu_si256((__m256i *)tab1);
    __m256i tmp = _mm256_set1_epi32(i);
    __m256i u = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    tmp = _mm256_add_epi32(tmp, u);
    __m256i v_pas = _mm256_set1_epi32(pas);
    tmp = _mm256_mullo_epi32(tmp, v_pas);
    __m256i y = _mm256_i32gather_epi32((const int *)tab2, tmp, 4);

    __m256i p = _mm256_set1_epi32(NB_P);
    __m256i q = _mm256_set1_epi32(NB_Q);

    __m256i a_low = _mm256_mul_epu32(x, y);
    __m256i b_low = _mm256_srli_epi64(a_low, NB_S);
    __m256i c_low = _mm256_srli_epi64(_mm256_mul_epu32(b_low, q), NB_T);

    __m256i a_high = _mm256_mul_epu32(_mm256_srli_si256(x, 4), _mm256_srli_si256(y, 4));
    __m256i b_high = _mm256_srli_epi64(a_high, NB_S);
    __m256i c_high = _mm256_srli_epi64(_mm256_mul_epu32(b_high, q), NB_T);

    __m256i d_low = _mm256_sub_epi64(a_low, _mm256_mul_epu32(c_low, p));
    __m256i d_high = _mm256_sub_epi64(a_high, _mm256_mul_epu32(c_high, p));
    __m256i d = _mm256_or_si256(d_low, _mm256_slli_si256(d_high, 4));

    __m256i result = _mm256_min_epu32(d, _mm256_sub_epi32(d, p));
    _mm256_storeu_si256((__m256i *)res, result);
}

void vect_mod_div_n(Uint *tab, Uint n)
{
    Uint n_inv = inv(n);
    __m256i y = _mm256_set1_epi32(n_inv);
    for (int i = 0; i < n; i += 8)
    {
        __m256i x = _mm256_loadu_si256((__m256i *)&tab[i]);
        __m256i p = _mm256_set1_epi32(NB_P);
        __m256i q = _mm256_set1_epi32(NB_Q);

        __m256i a_low = _mm256_mul_epu32(x, y);
        __m256i b_low = _mm256_srli_epi64(a_low, NB_S);
        __m256i c_low = _mm256_srli_epi64(_mm256_mul_epu32(b_low, q), NB_T);

        __m256i a_high = _mm256_mul_epu32(_mm256_srli_si256(x, 4), _mm256_srli_si256(y, 4));
        __m256i b_high = _mm256_srli_epi64(a_high, NB_S);
        __m256i c_high = _mm256_srli_epi64(_mm256_mul_epu32(b_high, q), NB_T);

        __m256i d_low = _mm256_sub_epi64(a_low, _mm256_mul_epu32(c_low, p));
        __m256i d_high = _mm256_sub_epi64(a_high, _mm256_mul_epu32(c_high, p));
        __m256i d = _mm256_or_si256(d_low, _mm256_slli_si256(d_high, 4));

        __m256i result = _mm256_min_epu32(d, _mm256_sub_epi32(d, p));
        _mm256_storeu_si256((__m256i *)&tab[i], result);
    }
}
