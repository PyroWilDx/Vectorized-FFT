#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "polynome.h"
#include "vect.h"

void afficher_poly(Poly P, Uint max)
{
    printf("%d", (P.coeffs)[0]);
    for (int i = 1; i <= max; i++)
    {
        if ((P.coeffs)[i] != 0)
        {
            printf(" + %dx^%d", (P.coeffs)[i], i);
        }
    }
    printf("\n");
}

Poly creer_poly(Uint deg)
{
    Poly P;
    P.deg = deg;
    P.coeffs = (Uint *)malloc(sizeof(Uint) * (deg + 1));
    for (Uint i = 0; i <= P.deg; i++)
    {
        P.coeffs[i] = 0;
    }
    return P;
}

Poly gen_poly(Uint deg)
{
    Poly P = creer_poly(deg);
    for (int i = 0; i <= deg; i++)
    {
        (P.coeffs)[i] = (rand() % (NB_P - 1)) + 1;
    }
    return P;
}

void liberer_poly(Poly P)
{
    free(P.coeffs);
}

int compare_poly(Poly P, Poly Q)
{
    if (P.deg != Q.deg)
    {
        printf("Degré différents : %d vs %d\n", P.deg, Q.deg);
        return 0;
    }
    for (int i = 0; i <= P.deg; i++)
    {
        if ((P.coeffs)[i] != (Q.coeffs)[i])
        {
            printf("Problème : Coeff1 = %d, Coeff2 = %d, pow = x^%d\n", (P.coeffs)[i], (Q.coeffs)[i], i);
            return 0;
        }
    }
    return 1;
}

Poly prod_poly_naif(Poly P, Poly Q)
{
    Poly R = creer_poly(P.deg + Q.deg);
    for (int i = 0; i <= P.deg; i++)
    {
        for (int j = 0; j <= Q.deg; j++)
        {
            Uint tmp = ((long)P.coeffs[i] * Q.coeffs[j]) % NB_P;
            (R.coeffs)[i + j] = ((R.coeffs)[i + j] + tmp) % NB_P;
        }
    }
    return R;
}

Poly copy_poly(Poly P, Uint deb, Uint fin)
{
    Poly R;
    R.coeffs = (Uint *)malloc(sizeof(Uint) * (fin - deb + 1));
    if (fin - deb < P.deg)
    {
        R.deg = fin - deb;
    }
    else
    {
        R.deg = P.deg;
    }
    for (int i = deb; i <= fin; i++)
    {
        (R.coeffs)[i - deb] = (P.coeffs)[i];
    }
    return R;
}

Poly somme_poly(Poly P, Poly Q, Uint deb)
{
    Poly R = copy_poly(P, 0, P.deg);
    for (int i = 0; i <= Q.deg; i++)
    {
        (R.coeffs)[deb + i] = ((R.coeffs)[deb + i] + (Q.coeffs)[i]) % NB_P;
    }
    return R;
}

Poly oppose_poly(Poly P)
{
    Poly R = copy_poly(P, 0, P.deg);
    for (int i = 0; i <= R.deg; i++)
    {
        (R.coeffs)[i] = -(R.coeffs)[i] + NB_P;
    }
    return R;
}

Poly prod_poly_karatsuba(Poly P, Poly Q)
{
    if (P.deg <= 130)
    {
        return prod_poly_naif(P, Q);
    }
    Uint k = (int)ceil((P.deg + 1) / 2.0);
    Poly P1 = copy_poly(P, 0, k - 1);
    Poly P2 = copy_poly(P, k, P.deg);
    Poly Q1 = copy_poly(Q, 0, k - 1);
    Poly Q2 = copy_poly(Q, k, Q.deg);
    Poly E1 = prod_poly_karatsuba(P1, Q1);
    Poly E2 = prod_poly_karatsuba(P2, Q2);
    Poly E3 = prod_poly_karatsuba(somme_poly(P1, P2, 0), somme_poly(Q1, Q2, 0));

    Poly R;
    R.deg = P.deg + Q.deg;
    R.coeffs = (Uint *)calloc(sizeof(Uint), R.deg + 1);

    R = somme_poly(R, E1, 0);
    R = somme_poly(R, E2, 2 * ((int)((P.deg + 2) / 2.0)));
    R = somme_poly(R, somme_poly(E3, oppose_poly(somme_poly(E1, E2, 0)), 0), k);

    return R;
}

Uint horner(Poly P, Uint x)
{
    Uint res = P.coeffs[P.deg];
    for (int i = P.deg - 1; i >= 0; i--)
    {
        res = mod_add(mod_mult(res, x), P.coeffs[i]);
    }
    return res;
}

Uint mod_pow(Uint x, Uint n)
{
    if (n == 0)
    {
        return 1;
    }
    if (n == 1)
    {
        return x;
    }
    Uint res;
    if (n % 2 == 0)
    {
        res = mod_pow(x, n / 2);
        return mod_mult(res, res);
    }
    res = mod_pow(x, n - 1);
    return mod_mult(res, x);
}

Uint *get_racines(Uint racine, Uint n)
{
    Uint *racines = (Uint *)malloc(sizeof(Uint) * n);
    racines[0] = 1;
    for (int i = 1; i < n; i++)
    {
        racines[i] = mod_mult(racines[i - 1], racine);
    }
    return racines;
}

Uint *eval_malloc(Poly P, Uint *racines)
{
    if (P.deg == 0)
    {
        Uint *tmp = (Uint *)malloc(sizeof(Uint));
        tmp[0] = (P.coeffs)[0];
        return tmp;
    }
    Uint k = (P.deg + 1) / 2;
    Poly R0 = creer_poly(k - 1); // creer_poly() fait un malloc
    Poly R1 = creer_poly(k - 1);

    Uint tmp;
    Uint *racines_bis = (Uint *)malloc(sizeof(Uint) * k);
    for (int i = 0; i < k; i++)
    {
        (R0.coeffs)[i] = mod_add((P.coeffs)[i], (P.coeffs)[i + k]);
        tmp = mod_sub((P.coeffs)[i], (P.coeffs)[i + k]);
        (R1.coeffs)[i] = mod_mult(tmp, racines[i]);
        racines_bis[i] = racines[2 * i];
    }
    Uint *r0 = eval_malloc(R0, racines_bis);
    Uint *r1 = eval_malloc(R1, racines_bis);
    Uint *res = (Uint *)malloc(sizeof(Uint) * 2 * k);
    for (int i = 0; i < k; i++)
    {
        res[2 * i] = r0[i];
        res[2 * i + 1] = r1[i];
    }
    return res;
}

Uint *eval(Uint *coeffs, Uint n, Uint *tmp_coeffs, Uint *racines, Uint pas_rac)
{
    if (n == 1)
        return &coeffs[0];

    Uint k = n / 2;

    Uint tmp;
    for (int i = 0; i < k; i++)
    {
        tmp_coeffs[i] = mod_add(coeffs[i], coeffs[i + k]);
        tmp = mod_sub(coeffs[i], coeffs[i + k]);
        tmp_coeffs[i + k] = mod_mult(tmp, racines[i * pas_rac]);
    }
    Uint *r0 = eval(tmp_coeffs, k, coeffs, racines, pas_rac * 2);
    Uint *r1 = eval(&tmp_coeffs[k], k, &coeffs[k], racines, pas_rac * 2);
    for (int i = 0; i < k; i++)
    {
        coeffs[2 * i] = r0[i];
        coeffs[2 * i + 1] = r1[i];
    }
    return coeffs;
}

Uint *vect_eval(Uint *coeffs, Uint n, Uint *tmp_coeffs, Uint *racines, Uint pas_rac, Uint *tmp_sub)
{
    if (n == 1)
        return &coeffs[0];

    Uint tmp;
    Uint k = n / 2;

    if (k >= 8)
    {
        for (Uint i = 0; i < k; i += 8)
        {
            vect_mod_add_sub_eval(&tmp_coeffs[i], &tmp_sub[i], &coeffs[i], &coeffs[i + k]);
            vect_mod_mult_eval(&tmp_coeffs[i + k], &tmp_sub[i], racines, i, pas_rac);
        }
    }
    else
    {
        for (Uint i = 0; i < k; i++)
        {
            tmp_coeffs[i] = mod_add(coeffs[i], coeffs[i + k]);
            tmp = mod_sub(coeffs[i], coeffs[i + k]);
            tmp_coeffs[i + k] = mod_mult(tmp, racines[i * pas_rac]);
        }
    }

    tmp = pas_rac * 2;
    Uint *r0 = vect_eval(tmp_coeffs, k, coeffs, racines, tmp, tmp_sub);
    Uint *r1 = vect_eval(&tmp_coeffs[k], k, &coeffs[k], racines, tmp, tmp_sub);

    for (Uint i = 0; i < k; i++)
    {
        tmp = 2 * i;
        coeffs[tmp] = r0[i];
        coeffs[tmp + 1] = r1[i];
    }
    return coeffs;
}

Uint *vect_eval_V2(Uint *coeffs, Uint n, Uint *tmp_coeffs, Uint *racines, Uint *tmp_sub)
{
    if (n == 1)
        return &coeffs[0];

    Uint tmp;
    Uint k = n / 2;

    if (k >= 8)
    {
        for (Uint i = 0; i < k; i += 8)
        {
            vect_mod_add_sub_eval(&tmp_coeffs[i], &tmp_sub[i], &coeffs[i], &coeffs[i + k]);
            vect_mod_mult(&tmp_coeffs[i + k], &tmp_sub[i], &racines[i]);
        }
    }
    else
    {
        for (Uint i = 0; i < k; i++)
        {
            tmp_coeffs[i] = mod_add(coeffs[i], coeffs[i + k]);
            tmp = mod_sub(coeffs[i], coeffs[i + k]);
            tmp_coeffs[i + k] = mod_mult(tmp, racines[i]);
        }
    }

    Uint *racines_bis = (Uint *)malloc(sizeof(Uint) * k);
    for (Uint i = 0; i < k; i++)
    {
        racines_bis[i] = racines[2 * i];
    }

    Uint *r0 = vect_eval_V2(tmp_coeffs, k, coeffs, racines_bis, tmp_sub);
    Uint *r1 = vect_eval_V2(&tmp_coeffs[k], k, &coeffs[k], racines_bis, tmp_sub);

    for (Uint i = 0; i < k; i++)
    {
        tmp = 2 * i;
        coeffs[tmp] = r0[i];
        coeffs[tmp + 1] = r1[i];
    }
    return coeffs;
}

Uint inv(Uint a)
{
    int r = a;
    int r_ = NB_P;
    int u = 1;
    int v = 0;
    int u_ = 0;
    int v_ = 1;
    while (r_ != 0)
    {
        int q = r / r_;
        int r2 = r;
        int u2 = u;
        int v2 = v;
        r = r_;
        u = u_;
        v = v_;
        r_ = r2 - q * r_;
        u_ = u2 - q * u_;
        v_ = v2 - q * v_;
    }
    if (u < 0)
    {
        u = u + NB_P;
    }
    return u;
}

Poly creer_poly_fft(Uint deg, Uint n)
{
    Poly P;
    P.deg = deg;
    P.coeffs = (Uint *)malloc(sizeof(Uint) * n);
    return P;
}

Poly gen_poly_fft(Uint deg, Uint n)
{
    Poly P = creer_poly_fft(deg, n);
    for (int i = 0; i <= deg; i++)
    {
        (P.coeffs)[i] = rand() % NB_P;
    }
    for (int i = deg + 1; i < n; i++)
    {
        (P.coeffs)[i] = 0;
    }
    return P;
}

Poly FFT(Poly P, Poly Q, Uint n, Uint racine_principale)
{

    // Étape 1 : Pré-Calcul
    Uint *racines = get_racines(racine_principale, n);

    // Étape 2 : Évaluation de P et Q
    Uint *tmp_coeffs = (Uint *)malloc(sizeof(Uint) * n);
    Uint *eval_P = eval(P.coeffs, n, tmp_coeffs, racines, 1);
    Uint *eval_Q = eval(Q.coeffs, n, tmp_coeffs, racines, 1);

    // Étape 3 : Produit point à point
    Uint *eval_R = (Uint *)malloc(sizeof(Uint) * n);
    for (Uint i = 0; i < n; i++)
    {
        eval_R[i] = mod_mult(eval_P[i], eval_Q[i]);
    }

    // Étape 4 : Interpolation
    Uint racine_p_inv = inv(racine_principale);
    Uint *racines_inv = get_racines(racine_p_inv, n);
    Poly R = creer_poly_fft(P.deg + Q.deg, n);
    R.coeffs = eval(eval_R, n, tmp_coeffs, racines_inv, 1);
    Uint n_inv = inv(n);
    for (int i = 0; i < n; i++)
    {
        R.coeffs[i] = mod_mult(R.coeffs[i], n_inv);
    }

    return R;
}

Poly vect_FFT(Poly P, Poly Q, Uint n, Uint racine_principale)
{

    // Étape 1 : Pré-Calcul
    Uint *racines = get_racines(racine_principale, n);

    // Étape 2 : Évaluation de P et Q
    Uint *tmp_coeffs = (Uint *)malloc(sizeof(Uint) * n);
    Uint *tmp_sub = (Uint *)malloc(sizeof(Uint) * n);
    Uint *eval_P = vect_eval(P.coeffs, n, tmp_coeffs, racines, 1, tmp_sub);
    Uint *eval_Q = vect_eval(Q.coeffs, n, tmp_coeffs, racines, 1, tmp_sub);

    // Étape 3 : Produit point à point
    Uint *eval_R = (Uint *)malloc(sizeof(Uint) * n);
    for (Uint i = 0; i < n; i += 8)
    {
        vect_mod_mult(&eval_R[i], &eval_P[i], &eval_Q[i]);
    }

    // Étape 4 : Interpolation
    Uint racine_p_inv = inv(racine_principale);
    Uint *racines_inv = get_racines(racine_p_inv, n);
    Poly R = creer_poly_fft(P.deg + Q.deg, n);
    R.coeffs = vect_eval(eval_R, n, tmp_coeffs, racines_inv, 1, tmp_sub);
    vect_mod_div_n(R.coeffs, n);

    return R;
}

Poly vect_FFT_V2(Poly P, Poly Q, Uint n, Uint racine_principale)
{

    // Étape 1 : Pré-Calcul
    Uint *racines = get_racines(racine_principale, n);

    // Étape 2 : Évaluation de P et Q
    Uint *tmp_coeffs = (Uint *)malloc(sizeof(Uint) * n);
    Uint *tmp_sub = (Uint *)malloc(sizeof(Uint) * n);
    Uint *eval_P = vect_eval_V2(P.coeffs, n, tmp_coeffs, racines, tmp_sub);
    Uint *eval_Q = vect_eval_V2(Q.coeffs, n, tmp_coeffs, racines, tmp_sub);

    // Étape 3 : Produit point à point
    Uint *eval_R = (Uint *)malloc(sizeof(Uint) * n);
    for (Uint i = 0; i < n; i += 8)
    {
        vect_mod_mult(&eval_R[i], &eval_P[i], &eval_Q[i]);
    }

    // Étape 4 : Interpolation
    Uint racine_p_inv = inv(racine_principale);
    Uint *racines_inv = get_racines(racine_p_inv, n);
    Poly R = creer_poly_fft(P.deg + Q.deg, n);
    R.coeffs = vect_eval_V2(eval_R, n, tmp_coeffs, racines_inv, tmp_sub);
    vect_mod_div_n(R.coeffs, n);

    return R;
}
