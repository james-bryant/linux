/*

  fp_trig.c: floating-point math routines for the Linux-m68k
  floating point emulator.

  Copyright (c) 1998-1999 David Huggins-Daines / Roman Zippel.

  I hereby give permission, free of charge, to copy, modify, and
  redistribute this software, in source or binary form, provided that
  the above copyright notice and the following disclaimer are included
  in all such copies.

  THIS SOFTWARE IS PROVIDED "AS IS", WITH ABSOLUTELY NO WARRANTY, REAL
  OR IMPLIED.

*/

/*
 * This file contains a modified version of parts of Motorola's
 * floating-point package for the 68040 (arch/m68k/fpsp040): fp_fsin(),
 * fp_fcos(), fp_ftan(), fp_fasin(), fp_facos(), fp_fatan(), fp_fsinh(),
 * fp_fcosh() and fp_fsincos0() to fp_fsincos7(), with the functions and
 * the constants above them that they use, are the package's algorithms
 * written in C, and the constants are the package's in another form.
 * The package comes with this notice (arch/m68k/fpsp040/README):
 *
 *	MOTOROLA MICROPROCESSOR & MEMORY TECHNOLOGY GROUP
 *	M68000 Hi-Performance Microprocessor Division
 *	M68040 Software Package
 *
 *	M68040 Software Package Copyright (c) 1993, 1994 Motorola Inc.
 *	All rights reserved.
 *
 *	THE SOFTWARE is provided on an "AS IS" basis and without warranty.
 *	To the maximum extent permitted by applicable law,
 *	MOTOROLA DISCLAIMS ALL WARRANTIES WHETHER EXPRESS OR IMPLIED,
 *	INCLUDING IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 *	PARTICULAR PURPOSE and any warranty against infringement with
 *	regard to the SOFTWARE (INCLUDING ANY MODIFIED VERSIONS THEREOF)
 *	and any accompanying written materials.
 *
 *	To the maximum extent permitted by applicable law,
 *	IN NO EVENT SHALL MOTOROLA BE LIABLE FOR ANY DAMAGES WHATSOEVER
 *	(INCLUDING WITHOUT LIMITATION, DAMAGES FOR LOSS OF BUSINESS
 *	PROFITS, BUSINESS INTERRUPTION, LOSS OF BUSINESS INFORMATION, OR
 *	OTHER PECUNIARY LOSS) ARISING OF THE USE OR INABILITY TO USE THE
 *	SOFTWARE.  Motorola assumes no responsibility for the maintenance
 *	and support of the SOFTWARE.
 *
 *	You are hereby granted a copyright license to use, modify, and
 *	distribute the SOFTWARE so long as this entire notice is retained
 *	without alteration in any modified and/or redistributed versions,
 *	and that such modified versions are clearly identified as such.
 *	No licenses are granted by implication, estoppel or otherwise
 *	under any patents or trademarks of Motorola, Inc.
 *
 * fp_ftanh() and fp_fatanh() are not taken from the package.
 */

#include "fp_emu.h"
#include "fp_trans.h"
#include "fp_trig.h"

/*
 * The trigonometric instructions follow ssin and stan of Motorola's
 * floating-point package for the 68040 (arch/m68k/fpsp040/ssin.S and
 * stan.S; its notice is at the head of this file) step by step, so that
 * they compute the same digits.  The package states an error below one
 * unit in the last place for the sine and the cosine of an argument
 * below 15 pi, and below three units for the tangent.  It does not
 * keep these bounds, even against the function of the argument as its
 * own pi/2 reduces it: results almost two units off were found for the
 * sine and the cosine, and three and a half for the tangent.
 *
 * Pi/2 has 66 bits in the package.  The sine, the cosine and the
 * tangent of an argument next to a multiple of pi/2 and those of a
 * large argument are those of the argument as these 66 bits reduce it,
 * and no more accurate than that.
 *
 * The constants are the package's.  It has most of them in single or
 * double precision, and they are the same numbers here:
 *
 *	2/pi	0x3fe45f30 0x6dc9c883
 *	A3	0xbf2a01a0 0x1a018b59	B1	0xbf000000
 *	A4	0x3ec71de3 0xa5341531	B4	0x3efa01a0 0x1a01d423
 *	A5	0xbe5ae645 0x2a118ae4	B5	0xbe927e4f 0xb79d9fcf
 *	A6	0x3de61209 0x7aae8da1	B6	0x3e21eed9 0x0612c972
 *	A7	0xbd6aaa77 0xccc994f5	B7	0xbda9396f 0x9f45ac19
 *					B8	0x3d2ac4d0 0xd6011ee3
 *	P3	0xbef2baa5 0xa8924f04	Q3	0xbf346f59 0xb39ba65f
 *					Q4	0x3ea0b759 0xf50f8688
 */
static const struct fp_ext fp_trig_twobypi =
	FPT_EXT(0, 0x3ffe, 0xa2f9836e, 0x4e441800);

static const struct fp_ext fp_sin_a1 =
	FPT_EXT(1, 0x3ffc, 0xaaaaaaaa, 0xaaaaaa99);
static const struct fp_ext fp_sin_a2 =
	FPT_EXT(0, 0x3ff8, 0x88888888, 0x888859af);
static const struct fp_ext fp_sin_a3 =
	FPT_EXT(1, 0x3ff2, 0xd00d00d0, 0x0c5ac800);
static const struct fp_ext fp_sin_a4 =
	FPT_EXT(0, 0x3fec, 0xb8ef1d29, 0xa0a98800);
static const struct fp_ext fp_sin_a5 =
	FPT_EXT(1, 0x3fe5, 0xd7322950, 0x8c572000);
static const struct fp_ext fp_sin_a6 =
	FPT_EXT(0, 0x3fde, 0xb0904bd5, 0x746d0800);
static const struct fp_ext fp_sin_a7 =
	FPT_EXT(1, 0x3fd6, 0xd553be66, 0x4ca7a800);

static const struct fp_ext fp_cos_b1 =
	FPT_EXT(1, 0x3ffe, 0x80000000, 0x00000000);
static const struct fp_ext fp_cos_b2 =
	FPT_EXT(0, 0x3ffa, 0xaaaaaaaa, 0xaaaaab5e);
static const struct fp_ext fp_cos_b3 =
	FPT_EXT(1, 0x3ff5, 0xb60b60b6, 0x0b61d438);
static const struct fp_ext fp_cos_b4 =
	FPT_EXT(0, 0x3fef, 0xd00d00d0, 0x0ea11800);
static const struct fp_ext fp_cos_b5 =
	FPT_EXT(1, 0x3fe9, 0x93f27dbc, 0xecfe7800);
static const struct fp_ext fp_cos_b6 =
	FPT_EXT(0, 0x3fe2, 0x8f76c830, 0x964b9000);
static const struct fp_ext fp_cos_b7 =
	FPT_EXT(1, 0x3fda, 0xc9cb7cfa, 0x2d60c800);
static const struct fp_ext fp_cos_b8 =
	FPT_EXT(0, 0x3fd2, 0xd62686b0, 0x08f71800);

static const struct fp_ext fp_tan_p1 =
	FPT_EXT(1, 0x3ffc, 0x8895a6c5, 0xfb423bca);
static const struct fp_ext fp_tan_p2 =
	FPT_EXT(0, 0x3ff6, 0xe073d3fc, 0x199c4a00);
static const struct fp_ext fp_tan_p3 =
	FPT_EXT(1, 0x3fef, 0x95d52d44, 0x92782000);
static const struct fp_ext fp_tan_q1 =
	FPT_EXT(1, 0x3ffd, 0xeef57e0d, 0xa84bc8ce);
static const struct fp_ext fp_tan_q2 =
	FPT_EXT(0, 0x3ff9, 0xd23cd684, 0x15d95fa1);
static const struct fp_ext fp_tan_q3 =
	FPT_EXT(1, 0x3ff3, 0xa37acd9c, 0xdd32f800);
static const struct fp_ext fp_tan_q4 =
	FPT_EXT(0, 0x3fea, 0x85bacfa8, 0x7c344000);

/*
 * 2^-40 in the form of fpt_compact().  Below it ssin, stan, satan and
 * stanh of the package return their argument and scos returns 1, less
 * a tiny term, each as the program rounds it.
 */
#define FP_TRIG_TINY	0x3fd78000

/* More passes than fp_trig_reducex() can need: see there. */
#define FP_TRIG_PASSES	600

/*
 * The argument reduction for 15 pi and more (REDUCEX): r becomes
 * r - N pi/2, at most pi/4 in magnitude, and the low bits of N are
 * returned.
 *
 * Each pass takes a multiple of 2^L pi/2 off the remainder, which is
 * held as R + r with r below the last bit of R.  For an R of exponent K
 * a pass has L = K - 27 and leaves a remainder below 2^L pi/4, whose
 * exponent is K - 28 at most.  From 16383 the exponent is down to 28,
 * where the last pass has L = 0, after 585 passes at most.  The largest
 * number takes 560 passes, and 565 are the most that were found for a
 * number of the largest exponent (0x7ffe0000 0xd43fc644 0x9585fe0b).
 *
 * The end of the loop does not rest on that arithmetic alone: no task
 * may stay in the kernel for a slip in it.  After FP_TRIG_PASSES passes
 * that were not the last the reduction gives up and returns false: no
 * number of this remainder is to be taken for a result, and the
 * instruction returns the NaN of an operand error.
 */
static bool fp_trig_reducex(struct fp_ext *r, int *np)
{
	struct fp_ext rl = { .exp = 0 };
	struct fp_ext n, p, w, wl, c;
	unsigned int pass;
	int k, l;

	/*
	 * An argument this large could overflow in the first pass: take
	 * 2^16383 pi/2 off first, in two pieces of which the first is exact.
	 */
	if (fpt_compact(r) == 0x7ffeffff) {
		c = (struct fp_ext)FPT_EXT(!r->sign, 0x7ffe, 0xc90fdaa2, 0);
		w = (struct fp_ext)FPT_EXT(!r->sign, 0x7fdc, 0x85a308d3, 0);
		fpt_add(r, &c);
		rl = *r;
		fpt_add(r, &w);
		fpt_sub(&rl, r);
		fpt_add(&rl, &w);
	}

	for (pass = 1; ; pass++) {
		/* R is 2^K or more; the last pass has L = 0 */
		k = r->exp - 0x3fff;
		l = k <= 28 ? 0 : k - 27;

		/*
		 * N = R * 2^-L * 2/pi to the nearest integer: adding and
		 * subtracting 2^63 of the sign of R rounds it.
		 */
		c = (struct fp_ext)FPT_EXT(0, 0x3ffe - l, 0xa2f9836e,
					   0x4e44152a);
		n = *r;
		fpt_mul(&n, &c);
		c = (struct fp_ext)FPT_EXT(r->sign, 0x3fff + 63, 0x80000000, 0);
		fpt_add(&n, &c);
		fpt_sub(&n, &c);

		/* W = N * P1 and w = N * P2, with P1 + P2 = 2^L pi/2 */
		c = (struct fp_ext)FPT_EXT(0, 0x3fff + l, 0xc90fdaa2, 0);
		w = n;
		fpt_mul(&w, &c);
		c = (struct fp_ext)FPT_EXT(0, 0x3fdd + l, 0x85a308d3, 0);
		wl = n;
		fpt_mul(&wl, &c);

		/* P + p = W + w, with p below the last bit of P */
		p = w;
		fpt_add(&p, &wl);
		fpt_sub(&w, &p);

		/* A = R - P and a = r - p */
		fpt_sub(r, &p);
		fpt_add(&w, &wl);
		p = *r;
		fpt_sub(&rl, &w);

		/* the new R = A + a */
		fpt_add(r, &rl);
		if (!l)
			break;
		if (pass == FP_TRIG_PASSES)
			return false;

		/* and the new r = (A - R) + a */
		fpt_sub(&p, r);
		fpt_add(&rl, &p);

		/* let other tasks run: a pass is about 2600 instructions */
		if (!(pass & 7))
			cond_resched();
	}
	*np = fpt_to_int(&n);

	return true;
}

/*
 * Reduce the argument of FSIN, FCOS, FSINCOS and FTAN, which is 2^-40
 * or more in magnitude: r = x - N pi/2, at most pi/4 in magnitude, and
 * the low bits of N go to np.  False if the reduction failed, which no
 * argument makes it do: see fp_trig_reducex().
 */
static bool fp_trig_reduce(struct fp_ext *r, int *np,
			   const struct fp_ext *x)
{
	struct fp_ext n;
	int k;

	*r = *x;
	if (fpt_compact(x) >= 0x4004bc7e)
		return fp_trig_reducex(r, np);

	/* below 15 pi: N pi/2 from the table, as Y1 + Y2 */
	n = *x;
	fpt_mul(&n, &fp_trig_twobypi);
	k = fpt_to_int(&n);
	fpt_sub(r, &fpt_pitbl[k + 32][0]);
	fpt_sub(r, &fpt_pitbl[k + 32][1]);
	*np = k;

	return true;
}

/*
 * sin(r + k pi/2) for an r of at most pi/4 in magnitude: the sine or
 * the cosine of r by the polynomial of the package, with the sign of
 * the quadrant.  The last operation is fpt_last_add()'s, and res is
 * what it leaves.  r is not kept.
 */
static void fp_sin_poly(struct fp_ext *res, struct fp_ext *r, int k,
			struct fpt_env *env)
{
	struct fp_ext s, t, a, b;

	/* S = R * R and T = S * S */
	s = *r;
	fpt_mul(&s, r);
	t = s;
	fpt_mul(&t, &s);

	if (!(k & 1)) {
		/*
		 * With R' = R of the quadrant's sign, the sine is
		 * R' + R' * S * ([A1 + T * (A3 + T * (A5 + T * A7))] +
		 *		  [S * (A2 + T * (A4 + T * A6))])
		 */
		if (k & 2)
			r->sign = !r->sign;
		a = fp_sin_a7;
		b = fp_sin_a6;
		fpt_mul(&a, &t);
		fpt_mul(&b, &t);
		fpt_add(&a, &fp_sin_a5);
		fpt_add(&b, &fp_sin_a4);
		fpt_mul(&a, &t);
		fpt_mul(&b, &t);
		fpt_add(&a, &fp_sin_a3);
		fpt_add(&b, &fp_sin_a2);
		fpt_mul(&t, &a);
		fpt_mul(&b, &s);
		fpt_add(&t, &fp_sin_a1);
		fpt_mul(&s, r);
		fpt_add(&t, &b);
		fpt_mul(&s, &t);

		fpt_last_add(&s, r, env);
	} else {
		/*
		 * With S' = S of the quadrant's sign, the cosine is
		 * +-1 + S' * ([B1 + T * (B3 + T * (B5 + T * B7))] +
		 *	       [S * (B2 + T * (B4 + T * (B6 + T * B8)))])
		 */
		*r = s;
		r->sign = (k & 2) != 0;
		b = fp_cos_b8;
		a = fp_cos_b7;
		fpt_mul(&b, &t);
		fpt_mul(&a, &t);
		fpt_add(&b, &fp_cos_b6);
		fpt_add(&a, &fp_cos_b5);
		fpt_mul(&b, &t);
		fpt_mul(&a, &t);
		fpt_add(&b, &fp_cos_b4);
		fpt_add(&a, &fp_cos_b3);
		fpt_mul(&b, &t);
		fpt_mul(&t, &a);
		fpt_add(&b, &fp_cos_b2);
		fpt_add(&t, &fp_cos_b1);
		fpt_mul(&s, &b);
		fpt_add(&s, &t);
		fpt_mul(&s, r);

		fpt_pow2(&t, 0);
		t.sign = (k & 2) != 0;
		fpt_last_add(&s, &t, env);
	}
	*res = s;
}

struct fp_ext *fp_fsin(struct fp_ext *dest, struct fp_ext *src)
{
	struct fpt_env env;
	struct fp_ext r;
	int n;

	dprint(PINSTR, "fsin\n");

	if (fpt_special(FPT_FSIN, dest, src))
		return dest;

	/* below 2^-40 the sine is the argument, as the program rounds it */
	if (fpt_compact(src) < FP_TRIG_TINY)
		return fpt_computed(dest, src, NULL);

	fpt_enter(&env);
	if (!fp_trig_reduce(&r, &n, src))
		return fpt_operr(dest, &env);
	fp_sin_poly(&r, &r, n, &env);

	return fpt_computed(dest, &r, &env);
}

/*
 * The cosine of an argument below 2^-40: 1 less a tiny term, so that
 * the program's rounding mode decides.  This is what that subtraction
 * leaves, unrounded: the number below 1 and more than half a unit of
 * its last place.
 */
static const struct fp_ext fp_cos_small = {
	.lowmant = 0xff,
	.exp = 0x3ffe,
	.mant.m32 = { 0xffffffff, 0xffffffff },
};

struct fp_ext *fp_fcos(struct fp_ext *dest, struct fp_ext *src)
{
	struct fpt_env env;
	struct fp_ext r;
	int n;

	dprint(PINSTR, "fcos\n");

	if (fpt_special(FPT_FCOS, dest, src))
		return dest;

	if (fpt_compact(src) < FP_TRIG_TINY)
		return fpt_computed(dest, &fp_cos_small, NULL);

	/* the cosine is the sine a quarter turn on */
	fpt_enter(&env);
	if (!fp_trig_reduce(&r, &n, src))
		return fpt_operr(dest, &env);
	fp_sin_poly(&r, &r, n + 1, &env);

	return fpt_computed(dest, &r, &env);
}

struct fp_ext *fp_ftan(struct fp_ext *dest, struct fp_ext *src)
{
	struct fp_ext r, s, p, q;
	struct fpt_env env;
	int n;

	dprint(PINSTR, "ftan\n");

	if (fpt_special(FPT_FTAN, dest, src))
		return dest;

	/* below 2^-40 the tangent is the argument, as the program rounds it */
	if (fpt_compact(src) < FP_TRIG_TINY)
		return fpt_computed(dest, src, NULL);

	fpt_enter(&env);
	if (!fp_trig_reduce(&r, &n, src))
		return fpt_operr(dest, &env);

	/*
	 * With S = R * R, tan(R) is U / V with
	 * U = R + R * S * (P1 + S * (P2 + S * P3)) and
	 * V = 1 + S * (Q1 + S * (Q2 + S * (Q3 + S * Q4)))
	 */
	s = r;
	fpt_mul(&s, &r);
	q = fp_tan_q4;
	p = fp_tan_p3;
	fpt_mul(&q, &s);
	fpt_mul(&p, &s);
	fpt_add(&q, &fp_tan_q3);
	fpt_add(&p, &fp_tan_p2);
	fpt_mul(&q, &s);
	fpt_mul(&p, &s);
	fpt_add(&q, &fp_tan_q2);
	fpt_add(&p, &fp_tan_p1);
	fpt_mul(&q, &s);
	fpt_mul(&p, &s);
	fpt_add(&q, &fp_tan_q1);
	fpt_mul(&p, &r);
	fpt_mul(&q, &s);
	fpt_add(&p, &r);
	fpt_pow2(&s, 0);
	fpt_add(&q, &s);

	/* the last operation: tan(R), or -cot(R) = V / -U for an odd N */
	if (n & 1) {
		p.sign = !p.sign;
		fpt_last_div(&q, &p, &env);
		return fpt_computed(dest, &q, &env);
	}
	fpt_last_div(&p, &q, &env);

	return fpt_computed(dest, &p, &env);
}

/*
 * The inverse trigonometric instructions follow satan, sasin and sacos
 * of the package (satan.S, sasin.S and sacos.S) in the same way.  The
 * package states an error below 2 units in the last place for the arc
 * tangent and below 3 for the arc sine and the arc cosine, which it
 * computes by the arc tangent.  It does not keep these bounds for the
 * arc tangent and the arc sine: results more than 2 and more than 3
 * units off were found.
 *
 * The constants of satan.S, which has them in double precision:
 *
 *	A1	0xbfc2476f 0x4e1da28e	B1	0xbfd55555 0x55555555
 *	A2	0x4002ac69 0x34a26db3	B2	0x3fc99999 0x99998fa9
 *	A3	0xbff6687e 0x314987d8	B3	0xbfc24924 0x921872f9
 *	C1	0xbfd55555 0x55555536	B4	0x3fbc71c6 0x46940220
 *	C2	0x3fc99999 0x9996263e	B5	0xbfb744ee 0x7faf45db
 *	C3	0xbfc24924 0x827107b8	B6	0x3fb34444 0x7f876989
 *	C4	0x3fbc7187 0x962d1d7d
 *	C5	0xbfb70bf3 0x98539e6a
 */
static const struct fp_ext fp_atan_a1 =
	FPT_EXT(1, 0x3ffc, 0x923b7a70, 0xed147000);
static const struct fp_ext fp_atan_a2 =
	FPT_EXT(0, 0x4000, 0x956349a5, 0x136d9800);
static const struct fp_ext fp_atan_a3 =
	FPT_EXT(1, 0x3fff, 0xb343f18a, 0x4c3ec000);

static const struct fp_ext fp_atan_b1 =
	FPT_EXT(1, 0x3ffd, 0xaaaaaaaa, 0xaaaaa800);
static const struct fp_ext fp_atan_b2 =
	FPT_EXT(0, 0x3ffc, 0xcccccccc, 0xcc7d4800);
static const struct fp_ext fp_atan_b3 =
	FPT_EXT(1, 0x3ffc, 0x92492490, 0xc397c800);
static const struct fp_ext fp_atan_b4 =
	FPT_EXT(0, 0x3ffb, 0xe38e3234, 0xa0110000);
static const struct fp_ext fp_atan_b5 =
	FPT_EXT(1, 0x3ffb, 0xba2773fd, 0x7a2ed800);
static const struct fp_ext fp_atan_b6 =
	FPT_EXT(0, 0x3ffb, 0x9a2223fc, 0x3b4c4800);

static const struct fp_ext fp_atan_c1 =
	FPT_EXT(1, 0x3ffd, 0xaaaaaaaa, 0xaaa9b000);
static const struct fp_ext fp_atan_c2 =
	FPT_EXT(0, 0x3ffc, 0xcccccccc, 0xb131f000);
static const struct fp_ext fp_atan_c3 =
	FPT_EXT(1, 0x3ffc, 0x92492413, 0x883dc000);
static const struct fp_ext fp_atan_c4 =
	FPT_EXT(0, 0x3ffb, 0xe38c3cb1, 0x68ebe800);
static const struct fp_ext fp_atan_c5 =
	FPT_EXT(1, 0x3ffb, 0xb85f9cc2, 0x9cf35000);

/* pi/2 in 64 bits, and the smallest normalized number */
static const struct fp_ext fp_atan_piby2 =
	FPT_EXT(0, 0x3fff, 0xc90fdaa2, 0x2168c235);
static const struct fp_ext fp_atan_tiny =
	FPT_EXT(0, 0x0001, 0x80000000, 0x00000000);

/*
 * The arc tangent of a normalized number as the sum of two terms: the
 * sum is the last operation of satan, and it is left to the caller,
 * who knows how it is to be rounded.  x is neither res nor last.
 */
static void fp_atan_terms(struct fp_ext *res, struct fp_ext *last,
			  const struct fp_ext *x)
{
	unsigned int compact = fpt_compact(x);
	struct fp_ext u, y, z, p, q;

	if (compact >= 0x3ffb8000 && compact <= 0x4002ffff) {
		/*
		 * From 1/16 to 16: with F the first five bits of X and a
		 * sixth that is set, atan(X) = atan(F) + atan(U) for
		 * U = (X - F) / (1 + X * F), atan(F) from the table of
		 * eight binades of sixteen F each
		 */
		q = *x;
		q.mant.m32[0] = (q.mant.m32[0] & 0xf8000000) | 0x04000000;
		q.mant.m32[1] = 0;
		u = *x;
		y = *x;
		fpt_mul(&y, &q);
		fpt_sub(&u, &q);
		fpt_pow2(&p, 0);
		fpt_add(&y, &p);
		fpt_div(&u, &y);
		*last = fpt_atantbl[(x->exp - 0x3ffb) << 4 |
				    (x->mant.m32[0] >> 27 & 15)];
		last->sign = x->sign;

		/*
		 * With V = U * U, atan(U) is
		 * U + A1 * U * V * (A2 + V * (A3 + V))
		 */
		y = u;
		fpt_mul(&y, &u);
		p = fp_atan_a3;
		fpt_add(&p, &y);
		fpt_mul(&p, &y);
		fpt_mul(&y, &u);
		fpt_add(&p, &fp_atan_a2);
		fpt_mul(&y, &fp_atan_a1);
		fpt_mul(&y, &p);
		fpt_add(&u, &y);
		*res = u;
	} else if (compact < FP_TRIG_TINY) {
		/* below 2^-40 the arc tangent is the argument: add nothing */
		*res = *x;
		last->lowmant = 0;
		last->sign = x->sign;
		last->exp = 0;
		last->mant.m64 = 0;
	} else if (compact < 0x3ffb8000) {
		/*
		 * Below 1/16: with Y = X * X and Z = Y * Y, atan(X) is
		 * X + X * Y * ([B1 + Z * (B3 + Z * B5)] +
		 *		[Y * (B2 + Z * (B4 + Z * B6))])
		 */
		y = *x;
		fpt_mul(&y, x);
		z = y;
		fpt_mul(&z, &y);
		p = fp_atan_b6;
		q = fp_atan_b5;
		fpt_mul(&p, &z);
		fpt_mul(&q, &z);
		fpt_add(&p, &fp_atan_b4);
		fpt_add(&q, &fp_atan_b3);
		fpt_mul(&p, &z);
		fpt_mul(&z, &q);
		fpt_add(&p, &fp_atan_b2);
		fpt_add(&z, &fp_atan_b1);
		fpt_mul(&p, &y);
		fpt_mul(&y, x);
		fpt_add(&z, &p);
		fpt_mul(&y, &z);
		*res = y;
		*last = *x;
	} else if (compact <= 0x40638000) {
		/*
		 * From 16 on: atan(X) = sign(X) * pi/2 + atan(X') for
		 * X' = -1 / X.  With Y = X' * X' and Z = Y * Y, atan(X') is
		 * X' + X' * Y * ([C1 + Z * (C3 + Z * C5)] +
		 *		  [Y * (C2 + Z * C4)])
		 */
		fpt_pow2(&u, 0);
		u.sign = 1;
		fpt_div(&u, x);
		y = u;
		fpt_mul(&y, &u);
		z = y;
		fpt_mul(&z, &y);
		p = fp_atan_c5;
		q = fp_atan_c4;
		fpt_mul(&p, &z);
		fpt_mul(&q, &z);
		fpt_add(&p, &fp_atan_c3);
		fpt_add(&q, &fp_atan_c2);
		fpt_mul(&z, &p);
		fpt_mul(&q, &y);
		fpt_add(&z, &fp_atan_c1);
		fpt_mul(&y, &u);
		fpt_add(&z, &q);
		fpt_mul(&y, &z);
		fpt_add(&y, &u);
		*res = y;
		*last = fp_atan_piby2;
		last->sign = x->sign;
	} else {
		/* beyond 2^100: sign(X) * (pi/2 less a tiny term) */
		*res = fp_atan_piby2;
		res->sign = x->sign;
		*last = fp_atan_tiny;
		last->sign = !x->sign;
	}
}

struct fp_ext *fp_fasin(struct fp_ext *dest, struct fp_ext *src)
{
	struct fp_ext a, b, x;
	struct fpt_env env;

	dprint(PINSTR, "fasin\n");

	/* a normalized operand is below 1 in magnitude from here on */
	if (fpt_special(FPT_FASIN, dest, src))
		return dest;

	/* asin(x) = atan(x / sqrt((1 - x) * (1 + x))) */
	fpt_enter(&env);
	fpt_pow2(&a, 0);
	b = a;
	fpt_sub(&a, src);
	fpt_add(&b, src);
	fpt_mul(&a, &b);
	fpt_sqrt(&b, &a);
	x = *src;
	fpt_div(&x, &b);
	fp_atan_terms(&a, &b, &x);
	fpt_last_add(&a, &b, &env);

	return fpt_computed(dest, &a, &env);
}

struct fp_ext *fp_facos(struct fp_ext *dest, struct fp_ext *src)
{
	struct fp_ext a, b, x;
	struct fpt_env env;

	dprint(PINSTR, "facos\n");

	/* a normalized operand is below 1 in magnitude from here on */
	if (fpt_special(FPT_FACOS, dest, src))
		return dest;

	/* acos(x) = 2 * atan(sqrt((1 - x) / (1 + x))) */
	fpt_enter(&env);
	fpt_pow2(&b, 0);
	a = *src;
	a.sign = !a.sign;
	fpt_add(&a, &b);
	fpt_add(&b, src);
	fpt_div(&a, &b);
	fpt_sqrt(&x, &a);
	fp_atan_terms(&a, &b, &x);
	fpt_add(&a, &b);

	/* the last operation is the doubling */
	fpt_last_add(&a, &a, &env);

	return fpt_computed(dest, &a, &env);
}

struct fp_ext *fp_fatan(struct fp_ext *dest, struct fp_ext *src)
{
	struct fp_ext a, b;
	struct fpt_env env;

	dprint(PINSTR, "fatan\n");

	if (fpt_special(FPT_FATAN, dest, src))
		return dest;

	fpt_enter(&env);
	fp_atan_terms(&a, &b, src);
	fpt_last_add(&a, &b, &env);

	return fpt_computed(dest, &a, &env);
}

/*
 * 16381 log 2 in two parts: T1 and T2 of ssinh.S and scosh.S, which
 * have them in double precision as 0x40c62d38 0xd3d64634 and
 * 0x3d6f90ae 0xb1e75cc7.
 */
static const struct fp_ext fp_hyp_t1 =
	FPT_EXT(0, 0x400c, 0xb169c69e, 0xb231a000);
static const struct fp_ext fp_hyp_t2 =
	FPT_EXT(0, 0x3fd6, 0xfc85758f, 0x3ae63800);

/*
 * FSINH follows ssinh of Motorola's floating-point package for the
 * 68040 (arch/m68k/fpsp040/ssinh.S; its notice is at the head of this
 * file): with Z = e^|x| - 1,
 *
 *	sinh x = sign(x) * (Z + Z / (1 + Z)) / 2
 *
 * and from 16380 log 2 on, which is as far as the package takes e^|x|
 * itself (it would overflow from 16384 log 2 on, before the result
 * does), sign(x) * 2^16380 * e^(|x| - 16381 log 2).  The package states
 * an error below 3 units in the last place.
 */
struct fp_ext *fp_fsinh(struct fp_ext *dest, struct fp_ext *src)
{
	unsigned int compact;
	struct fpt_env env;
	struct fp_ext y, z, r;

	dprint(PINSTR, "fsinh\n");

	if (fpt_special(FPT_FSINH, dest, src))
		return dest;

	fpt_enter(&env);
	compact = fpt_compact(src);
	y = *src;
	y.sign = 0;
	if (compact > 0x400cb2b3) {
		fpt_overflow(&r, src->sign, &env);
	} else if (compact > 0x400cb167) {
		fpt_sub(&y, &fp_hyp_t1);
		fpt_sub(&y, &fp_hyp_t2);
		fp_etox(&r, &y, NULL);
		fpt_pow2(&z, 16380);
		z.sign = src->sign;
		fpt_last_mul(&r, &z, &env);
	} else {
		fp_etoxm1(&z, &y, NULL);
		y = z;
		fpt_pow2(&r, 0);
		fpt_add(&y, &r);
		r = z;
		fpt_div(&r, &y);
		fpt_add(&r, &z);
		fpt_pow2(&z, -1);
		z.sign = src->sign;
		fpt_last_mul(&r, &z, &env);
	}

	return fpt_computed(dest, &r, &env);
}

/*
 * FCOSH follows scosh of the package (arch/m68k/fpsp040/scosh.S): with
 * Z = e^|x|,
 *
 *	cosh x = Z / 2 + (1 / 4) / (Z / 2)
 *
 * and from 16380 log 2 on 2^16380 * e^(|x| - 16381 log 2).  The package
 * states an error below 3 units in the last place.
 */
struct fp_ext *fp_fcosh(struct fp_ext *dest, struct fp_ext *src)
{
	unsigned int compact;
	struct fpt_env env;
	struct fp_ext y, z, r;

	dprint(PINSTR, "fcosh\n");

	if (fpt_special(FPT_FCOSH, dest, src))
		return dest;

	fpt_enter(&env);
	compact = fpt_compact(src);
	y = *src;
	y.sign = 0;
	if (compact > 0x400cb2b3) {
		fpt_overflow(&r, 0, &env);
	} else if (compact > 0x400cb167) {
		fpt_sub(&y, &fp_hyp_t1);
		fpt_sub(&y, &fp_hyp_t2);
		fp_etox(&r, &y, NULL);
		fpt_pow2(&z, 16380);
		fpt_last_mul(&r, &z, &env);
	} else {
		fp_etox(&r, &y, NULL);
		fpt_pow2(&z, -1);
		fpt_mul(&r, &z);
		fpt_pow2(&z, -2);
		fpt_div(&z, &r);
		fpt_last_add(&r, &z, &env);
	}

	return fpt_computed(dest, &r, &env);
}

struct fp_ext *fp_ftanh(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("ftanh\n");

	fp_monadic_check(dest, src);

	return dest;
}

struct fp_ext *fp_fatanh(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("fatanh\n");

	fp_monadic_check(dest, src);

	return dest;
}

/*
 * FSINCOS: the sine is the result of the instruction and the cosine
 * goes to the register creg.  ssincos of the package computes both with
 * polynomials in S alone, where FSIN and FCOS split theirs in two, so
 * the last bit of a result may not be that of FSIN or FCOS.
 */
static struct fp_ext *fp_fsincos(struct fp_ext *dest, struct fp_ext *src,
				 unsigned int creg)
{
	struct fp_ext *cdest = &FPDATA->fpreg[creg];
	struct fp_ext r, s, a, b, t;
	struct fpt_env env;
	int n;

	dprint(PINSTR, "fsincos\n");

	/* the cosine has an entry in the table wherever the sine has one */
	t = *src;
	if (fpt_special(FPT_FSINCOS, &s, src)) {
		fpt_special(FPT_FSINCOS_COS, &b, &t);
		fpt_second(cdest, &b, NULL);
		*dest = s;
		return dest;
	}

	/* below 2^-40 the sine is the argument, as the program rounds it */
	if (fpt_compact(src) < FP_TRIG_TINY) {
		fpt_second(cdest, &fp_cos_small, NULL);
		return fpt_computed(dest, src, NULL);
	}

	fpt_enter(&env);
	if (!fp_trig_reduce(&r, &n, src)) {
		fpt_operr(dest, &env);
		fpt_second(cdest, dest, NULL);
		return dest;
	}

	/*
	 * With S = R * R, sin(R) - R is
	 * R * S * (A1 + S * (A2 + S * (A3 + S * (A4 + S * (A5 + S * (A6 +
	 * S * A7)))))) and cos(R) - 1 is
	 * S * (B1 + S * (B2 + S * (B3 + S * (B4 + S * (B5 + S * (B6 + S *
	 * (B7 + S * B8)))))))
	 */
	s = r;
	fpt_mul(&s, &r);
	a = fp_sin_a7;
	b = fp_cos_b8;
	fpt_mul(&a, &s);
	fpt_mul(&b, &s);
	fpt_add(&a, &fp_sin_a6);
	fpt_add(&b, &fp_cos_b7);
	fpt_mul(&a, &s);
	fpt_mul(&b, &s);
	fpt_add(&a, &fp_sin_a5);
	fpt_add(&b, &fp_cos_b6);
	fpt_mul(&a, &s);
	fpt_mul(&b, &s);
	fpt_add(&a, &fp_sin_a4);
	fpt_add(&b, &fp_cos_b5);
	fpt_mul(&a, &s);
	fpt_mul(&b, &s);
	fpt_add(&a, &fp_sin_a3);
	fpt_add(&b, &fp_cos_b4);
	fpt_mul(&a, &s);
	fpt_mul(&b, &s);
	fpt_add(&a, &fp_sin_a2);
	fpt_add(&b, &fp_cos_b3);
	fpt_mul(&a, &s);
	fpt_mul(&b, &s);
	fpt_add(&a, &fp_sin_a1);
	fpt_add(&b, &fp_cos_b2);
	fpt_mul(&a, &s);
	fpt_mul(&b, &s);

	/*
	 * The quadrant: for an even N the sine and the cosine of the
	 * argument are sin(R) and cos(R), for an odd N they are cos(R)
	 * and -sin(R), and bit 1 of N changes the sign of both.
	 */
	s.sign = (n & 2) != 0;
	r.sign ^= ((n >> 1) ^ n) & 1;
	fpt_mul(&a, &r);
	fpt_add(&b, &fp_cos_b1);
	fpt_mul(&b, &s);

	/* the last operation of each: as the program rounds */
	fpt_last_add(&a, &r, &env);
	fpt_pow2(&t, 0);
	t.sign = s.sign;
	fpt_last_add(&b, &t, &env);
	if (n & 1) {
		fpt_second(cdest, &a, &env);
		return fpt_computed(dest, &b, &env);
	}
	fpt_second(cdest, &b, &env);
	return fpt_computed(dest, &a, &env);
}

struct fp_ext *fp_fsincos0(struct fp_ext *dest, struct fp_ext *src)
{
	return fp_fsincos(dest, src, 0);
}

struct fp_ext *fp_fsincos1(struct fp_ext *dest, struct fp_ext *src)
{
	return fp_fsincos(dest, src, 1);
}

struct fp_ext *fp_fsincos2(struct fp_ext *dest, struct fp_ext *src)
{
	return fp_fsincos(dest, src, 2);
}

struct fp_ext *fp_fsincos3(struct fp_ext *dest, struct fp_ext *src)
{
	return fp_fsincos(dest, src, 3);
}

struct fp_ext *fp_fsincos4(struct fp_ext *dest, struct fp_ext *src)
{
	return fp_fsincos(dest, src, 4);
}

struct fp_ext *fp_fsincos5(struct fp_ext *dest, struct fp_ext *src)
{
	return fp_fsincos(dest, src, 5);
}

struct fp_ext *fp_fsincos6(struct fp_ext *dest, struct fp_ext *src)
{
	return fp_fsincos(dest, src, 6);
}

struct fp_ext *fp_fsincos7(struct fp_ext *dest, struct fp_ext *src)
{
	return fp_fsincos(dest, src, 7);
}
