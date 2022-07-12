#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include "polynome.h"

#define MAX_PUISSANCE 24

clock_t temps_initial;
clock_t temps_final;
double temps_cpu;
double temps_tot_eval_malloc = 0;
double temps_tot_eval = 0;
double temps_tot_vect_eval = 0;

Poly test_naif(Poly P, Poly Q, int deg, double *temps) {
    temps_initial = clock();
    Poly R = prod_poly_naif(P, Q);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    *temps = temps_cpu;
    printf("Degré = %d, Naif : %f\n", deg, temps_cpu);
    return R;
}

Poly test_karatsuba(Poly P, Poly Q, int deg, double *temps) {
    temps_initial = clock();
    Poly R = prod_poly_karatsuba(P, Q);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    *temps = temps_cpu;
    printf("Degré = %d Karatsuba : %f\n", deg, temps_cpu);
    return R;
}

void compare_naif_karatsuba(int v) {
    Poly P, Q, R1, R2;
    double temps_naif, temps_karatsuba;
    for (int i = 32; i < 10000000; i *= 2) {
        P = gen_poly(i);
        Q = gen_poly(i);

        R1 = test_naif(P, Q, i, &temps_naif);
        R2 = test_karatsuba(P, Q, i, &temps_karatsuba);
        if (v == 1) assert(compare_poly(R1, R2));

        liberer_poly(P);
        liberer_poly(Q);
        liberer_poly(R1);
        liberer_poly(R2);

        printf("\n");
    }
}

void verif(Uint *res, Poly P, Uint *racines) {
    for (int i = 0; i < P.deg+1; i++) {
        int horner_x = horner(P, racines[i]);
        if (res[i] != horner_x) printf("i = %d, eval[i] = %d, horner[i] = %d\n", i, res[i], horner_x);
        assert(res[i] == horner_x);
    }
    // afficher_poly(P);
}

void test_eval_malloc(int deg, Uint racine_p, int aff, int v) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine_p, P.deg+1);
    temps_initial = clock();
    Uint *res = eval_malloc(P, racines);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    temps_tot_eval_malloc += temps_cpu;

    if (aff == 1) printf("Degré = %d | Temps eval() = %f\n", P.deg, temps_cpu);
    if (v == 1) verif(res, P_cpy, racines);
    
    liberer_poly(P);
    liberer_poly(P_cpy);
    free(racines);
}

void test_eval(int deg, Uint racine_p, int aff, int v) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine_p, P.deg+1);
    Uint *tmp_coeffs = (Uint *) malloc(sizeof(Uint)*(deg+1));
    temps_initial = clock();
    Uint *res = eval(P.coeffs, P.deg+1, tmp_coeffs, racines, 1);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    temps_tot_eval += temps_cpu;

    if (aff == 1) printf("Degré = %d | Temps eval() = %f\n", P.deg, temps_cpu);
    if (v == 1) verif(res, P_cpy, racines);

    liberer_poly(P);
    liberer_poly(P_cpy);
    free(racines);
    free(tmp_coeffs);
}

void test_vect_eval(Uint deg, Uint racine_p, int aff, int v) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine_p, P.deg+1);
    Uint *tmp_coeffs = (Uint *) malloc(sizeof(Uint)*(deg+1));
    Uint *tmp_sub = (Uint *) malloc(sizeof(Uint)*(deg+1));
    temps_initial = clock();
    Uint *res = vect_eval(P.coeffs, P.deg+1, tmp_coeffs, racines, 1, tmp_sub);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    temps_tot_vect_eval += temps_cpu;

    if (aff == 1) printf("Degré = %d | Temps vect_eval() = %f\n", P.deg, temps_cpu);
    if (v == 1) verif(res, P_cpy, racines);

    liberer_poly(P);
    liberer_poly(P_cpy);
    free(racines);
    free(tmp_coeffs);
    free(tmp_sub);
}

void test_vect_eval_V2(Uint deg, Uint racine_p, int aff, int v) {
    Poly P = gen_poly(deg);
    Poly P_cpy = copy_poly(P, 0, P.deg);
    Uint *racines = get_racines(racine_p, P.deg+1);
    Uint *tmp_coeffs = (Uint *) malloc(sizeof(Uint)*(deg+1));
    Uint *tmp_sub = (Uint *) malloc(sizeof(Uint)*(deg+1));
    temps_initial = clock();
    Uint *res = vect_eval_V2(P.coeffs, P.deg+1, tmp_coeffs, racines, tmp_sub);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
    temps_tot_vect_eval += temps_cpu;

    if (aff == 1) printf("Degré = %d | Temps vect_eval_V2() = %f\n", P.deg, temps_cpu);
    if (v == 1) verif(res, P_cpy, racines);

    liberer_poly(P);
    liberer_poly(P_cpy);
    free(racines);
    free(tmp_coeffs);
    free(tmp_sub);
}

void compare_ffts(Uint racine, Uint ordre_racine, int v) {
    for (int n = 32768; n <= pow(2, MAX_PUISSANCE); n *= 2) {
        Uint racine_principale = mod_pow(racine, ordre_racine/n);

        Poly P = gen_poly_fft(n/2-1, n);
        Poly vect_P = copy_poly(P, 0, n);
        Poly vect_P_V2 = copy_poly(P, 0, n);
        Poly Q = gen_poly_fft(n/2-1, n);
        Poly vect_Q = copy_poly(Q, 0, n);
        Poly vect_Q_V2 = copy_poly(Q, 0, n);

        temps_initial = clock();
        Poly R = FFT(P, Q, n, racine_principale);
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("n = %d, Temps FFT() : %f\n", n, temps_cpu);

        temps_initial = clock();
        Poly vect_R = vect_FFT(vect_P, vect_Q, n, racine_principale);
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("n = %d, Temps vect_FFT() : %f\n", n, temps_cpu);

        temps_initial = clock();
        Poly vect_R_V2 = vect_FFT_V2(vect_P_V2, vect_Q_V2, n, racine_principale);
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("n = %d, Temps vect_FFT_V2() : %f\n", n, temps_cpu);

        if (v == 1) assert(compare_poly(R, vect_R) && compare_poly(vect_R, vect_R_V2));

        printf("\n");
    }
}

void compare_TOUT(Uint racine, Uint ordre_racine, int v) {
    for (int n = pow(2, 16); n <= pow(2, MAX_PUISSANCE); n *= 2) {
        Uint racine_principale = mod_pow(racine, ordre_racine/n);

        Poly P_fft = gen_poly_fft(n/2-1, n);
        Poly P_naif = copy_poly(P_fft, 0, P_fft.deg);
        Poly P_kara = copy_poly(P_fft, 0, P_fft.deg);
        Poly P_vect = copy_poly(P_fft, 0, n);
        Poly P_vect_V2 = copy_poly(P_fft, 0, n);
        Poly Q_fft = gen_poly_fft(n/2-1, n);
        Poly Q_naif = copy_poly(Q_fft, 0, Q_fft.deg);
        Poly Q_kara = copy_poly(Q_fft, 0, Q_fft.deg);
        Poly Q_vect = copy_poly(Q_fft, 0, n);
        Poly Q_vect_V2 = copy_poly(Q_fft, 0, n);

        temps_initial = clock();
        Poly R_naif;
        if (n <= pow(2, 20)) {
            R_naif = prod_poly_naif(P_naif, Q_naif);
        }
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("n = %d, Temps naif() : %f\n", n, temps_cpu);

        temps_initial = clock();
        Poly R_kara;
        if (n <= pow(2, 21)) {
            R_kara = prod_poly_karatsuba(P_kara, Q_kara);
        }
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("n = %d, Temps karatsuba() : %f\n", n, temps_cpu);

        temps_initial = clock();
        Poly R_fft = FFT(P_fft, Q_fft, n, racine_principale);
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("n = %d, Temps FFT() : %f\n", n, temps_cpu);

        temps_initial = clock();
        Poly R_vect = vect_FFT(P_vect, Q_vect, n, racine_principale);
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("n = %d, Temps vect_FFT() : %f\n", n, temps_cpu);

        temps_initial = clock();
        Poly R_vect_V2 = vect_FFT_V2(P_vect_V2, Q_vect_V2, n, racine_principale);
        temps_final = clock();
        temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;
        printf("n = %d, Temps vect_FFT_V2() : %f\n", n, temps_cpu);

        if (v == 1) {
            if (n <= pow(2, 20)) {
                assert(compare_poly(R_naif, R_kara));
            }
            if (n <= pow(2, 21)) {
                assert(compare_poly(R_kara, R_fft));
            }
            assert(compare_poly(R_fft, R_vect));
            assert(compare_poly(R_vect, R_vect_V2));
        }
        printf("\n");
    }
}

int main() {
    srand(time(NULL));

    printf("Saisir un numéro : \n");
    printf("0 : Tests Naif vs Karatsuba\n");
    printf("1 : Tests des temps moyens de eval_malloc() vs eval()\n");
    printf("2 : Tests eval() vs vect_eval() vs vect_eval_V2()\n");
    printf("3 : Tests validité FFT (compraison résultats de Naif et de FFT)\n");
    printf("4 : Tests FFT() vs vect_FFT() vs vect_FFT_V2()\n");
    printf("5 : Tests comparaison de toutes les fonctions\n");
    printf("6 : Tests de la fonction inverse d'un entier inv() sur Z/pZ\n");

    int num;
    int _ = scanf("%d", &num);
    if (_ == 0) {
        return 0;
    }

    // Ancien nombre premier 2013265921
    // Uint racine = 2;
    // Uint ordre_racine = (NB_P-1)/2;
    Uint racine = 11;
    Uint ordre_racine = NB_P-1;

    Uint deg = 131071; // 2^k - 1
    Uint n = deg+1; // 2^k
    Uint racine_principale = mod_pow(racine, ordre_racine/n);

    if (num == 0) { // TESTS NAIF VS KARATSUBA
        compare_naif_karatsuba(0);
    }

    if (num == 1) { // TESTS EVAL_MALLOC() VS EVAL()
        Uint nb_tours = 1;
        for (Uint i = 32768; i <= pow(2, MAX_PUISSANCE); i = i*2) {
            Uint deg = i-1;
            Uint racine_principale = mod_pow(racine, ordre_racine/(deg+1));
            for (Uint j = 0; j < nb_tours; j++) {
                test_eval_malloc(deg, racine_principale, 0, 0);
                test_eval(deg, racine_principale, 0, 0);
            }
            printf(">>>> Degré = %d <<<<\n", deg);
            printf("=====> Temps moyen eval_malloc : %f\n", temps_tot_eval_malloc/nb_tours);
            printf("=====> Temps moyen eval : %f\n", temps_tot_eval/nb_tours);
        }
    }

    if (num == 2) { // TESTS EVAL() VS VECT_EVAL() VS VECT_EVAL_V2()
        for (Uint i = 32768; i <= pow(2, MAX_PUISSANCE); i = i*2) {
            Uint deg = i-1;
            Uint racine_principale = mod_pow(racine, ordre_racine/(deg+1));
            test_eval(deg, racine_principale, 1, 0);
            test_vect_eval(deg, racine_principale, 1, 0);
            test_vect_eval_V2(deg, racine_principale, 1, 0);
            printf("\n");
        }
    }

    if (num == 3) { // TEST VALIDATION FFT
        Poly P = gen_poly_fft(n/2-1, n);
        Poly P_cpy = copy_poly(P, 0, P.deg);
        Poly vect_P = copy_poly(P, 0, n);
        Poly vect_P_V2 = copy_poly(P, 0, n);
        Poly Q = gen_poly_fft(n/2-1, n);
        Poly Q_cpy = copy_poly(Q, 0, Q.deg);
        Poly vect_Q = copy_poly(Q, 0, n);
        Poly vect_Q_V2 = copy_poly(Q, 0, n);
        Poly R = FFT(P, Q, n, racine_principale);
        Poly vect_R = vect_FFT(vect_P, vect_Q, n, racine_principale);
        Poly vect_R_V2 = vect_FFT_V2(vect_P_V2, vect_Q_V2, n, racine_principale);
        Poly R_cpy = prod_poly_naif(P_cpy, Q_cpy);
        if (compare_poly(R, R_cpy) == 1 && compare_poly(vect_R, R_cpy) == 1 
            && compare_poly(vect_R_V2, R_cpy) == 1) {
            printf("OK\n");
        } else {
            printf("PAS OK\n");
        }
    }

    if (num == 4) { // TESTS FFT VS VECT_FFT VS VECT_FFT_V2
        compare_ffts(racine, ordre_racine, 0);
    }

    if (num == 5) { // TESTS NAIF VS KARATSUBA VS FFT VS VECT_FFT VS VECT_FFT_V2
        compare_TOUT(racine, ordre_racine, 0);
    }

    if (num == 6) { // TESTS INVERSE RACINE PRINCIPALE
        for (Uint i = 2; i <= pow(2, MAX_PUISSANCE); i = i*2) {
            Uint racine_principale = mod_pow(racine, ordre_racine/(deg+1));
            int racine_p_inv = inv(racine_principale);
            printf("n = %d, r = %d, r_inv = %d, r*r_inv mod p = %ld\n", i, racine_principale, 
                        racine_p_inv, ((long) racine_principale*racine_p_inv) % NB_P);
        }
    }

    printf("\n");
    return 0;
}
