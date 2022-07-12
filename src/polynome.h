#ifndef POLYNOME_H
#define POLYNOME_H

#include <math.h>
#include <immintrin.h>

// #define NB_P 2013265921 // Racine = 2, Ordre = (P-1)/2, R = 31 bits
#define NB_P 754974721 // Racine = 11, Ordre = P-1, R = 30 bits
                       // 754974721 = 1 + 2^24 * 3^2 * 5 => DegMax(Poly) = 2^24

#define NB_R log(NB_P) / log(2) // 2^(R-1) < P <=  2^R
// Choix de T et S tels que T >= R et S+T < n+R-1, ici n=32
#define NB_T 30
#define NB_S 27
#define NB_Q (int)(pow(2, NB_T + NB_S) / NB_P) // Q = floor(2^(S+T)/P)

typedef unsigned int Uint;

typedef struct _Poly
{
    Uint *coeffs;
    Uint deg;
} Poly;

inline Uint mod_add(Uint a, Uint b)
{
    Uint res = a + b;
    if (res < NB_P)
    {
        return res;
    }
    return res - NB_P;
}

inline Uint mod_sub(Uint a, Uint b)
{
    if (a < b)
    {
        return NB_P - (b - a);
    }
    return a - b;
}

inline Uint mod_mult(Uint a, Uint b)
{
    return ((unsigned long)a * b) % NB_P;
}

void afficher_poly(Poly P, Uint max);
Poly gen_poly(Uint deg);
void liberer_poly(Poly P);
int compare_poly(Poly P, Poly Q);
Poly prod_poly_naif(Poly P, Poly Q);
Poly copy_poly(Poly P, Uint deb, Uint fin);
Poly poly_somme(Poly P, Poly Q, Uint deb);
Poly oppose_poly(Poly P);
Poly prod_poly_karatsuba(Poly P, Poly Q);
Uint mod_pow(Uint x, Uint n);
Uint horner(Poly P, Uint x);
Uint *get_racines(Uint racine, Uint n);
Uint *eval_malloc(Poly P, Uint *racines);
Uint *eval(Uint *coeffs, Uint n, Uint *tmp_coeffs, Uint *racines, Uint pas_rac);
Uint *vect_eval(Uint *coeffs, Uint n, Uint *tmp_coeffs, Uint *racines, Uint pas_rac, Uint *tmp_sub);
Uint *vect_eval_V2(Uint *coeffs, Uint n, Uint *tmp_coeffs, Uint *racines, Uint *tmp_sub);
Uint inv(Uint a);
Poly creer_poly_fft(Uint deg, Uint n);
Poly gen_poly_fft(Uint deg, Uint n);
Poly FFT(Poly P, Poly Q, Uint n, Uint racine_principale);
Poly vect_FFT(Poly P, Poly Q, Uint n, Uint racine_principale);
Poly vect_FFT_V2(Poly P, Poly Q, Uint n, Uint racine_principale);

#endif
