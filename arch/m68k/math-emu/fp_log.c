/*

  fp_log.c: floating-point math routines for the Linux-m68k
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
 * floating-point package for the 68040 (arch/m68k/fpsp040): fp_fetox(),
 * fp_fetoxm1(), fp_ftwotox(), fp_ftentox() and fp_flogn(), with the
 * functions and the constants above them that they use, are the
 * package's algorithms written in C, and the constants are the package's
 * in another form.  The package comes with this notice
 * (arch/m68k/fpsp040/README):
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
 * fp_fsqrt(), fp_flognp1(), fp_flog10(), fp_flog2(), fp_fgetexp() and
 * fp_fgetman() are not taken from the package.
 */

#include "fp_emu.h"
#include "fp_log.h"
#include "fp_trans.h"
#include "multi_arith.h"

struct fp_ext *fp_fsqrt(struct fp_ext *dest, struct fp_ext *src)
{
	union fp_mant64 x, root, rem, trial;
	unsigned long low, remh, trialh;
	int i, exp;

	dprint(PINSTR, "fsqrt\n");

	fp_monadic_check(dest, src);

	if (IS_ZERO(dest))
		return dest;

	if (dest->sign) {
		fp_set_nan(dest);
		return dest;
	}
	if (IS_INF(dest))
		return dest;

	/* shift up the mantissa of a denormalized number */
	exp = dest->exp;
	if ((long)dest->mant.m32[0] >= 0)
		exp -= fp_overnormalize(dest);

	/*
	 * The mantissa m is an integer of 64 bits with the highest bit
	 * set, and the number is m * 2^(exp - 0x3fff - 63).  Let x be
	 * m * 2^63 for an odd exp and m * 2^64 for an even one.  Then
	 * 2^126 <= x < 2^128, and the square root of the number is
	 *
	 *	sqrt(x) * 2^((exp + 0x3fff) / 2 - 0x3fff - 63)
	 *
	 * with a division that rounds down.  The integer part of sqrt(x)
	 * has 64 bits with the highest bit set: the mantissa of the result.
	 */
	x = dest->mant;
	low = 0;
	if (exp & 1) {
		low = (x.m32[1] & 1) << 1;
		x.m64 >>= 1;
	}

	/*
	 * The root bit by bit, as in long division: each step shifts the
	 * next two bits of x, from x itself and then from low, into the
	 * remainder and subtracts 4 * root + 1 from it if it can.  The
	 * remainder is at most twice the root, so that it has up to 66
	 * bits with these two.
	 */
	root.m64 = 0;
	rem.m64 = 0;
	for (i = 0; i < 64; i++) {
		/*
		 * the two bits that the left shift takes off the top of the
		 * remainder, which stays below twice the root and so needs
		 * no more than these
		 */
		remh = rem.m32[0] >> 30;
		rem.m64 = rem.m64 << 2 | x.m32[0] >> 30;
		x.m64 = x.m64 << 2 | low;
		low = 0;

		trialh = root.m32[0] >> 30;
		trial.m64 = root.m64 << 2 | 1;
		root.m64 <<= 1;
		if (remh > trialh || (remh == trialh && rem.m64 >= trial.m64)) {
			remh -= trialh + (rem.m64 < trial.m64);
			rem.m64 -= trial.m64;
			root.m32[1] |= 1;
		}
	}

	dest->exp = (exp + 0x3fff) / 2;
	dest->mant = root;

	/*
	 * The low mantissa byte tells the rounding what follows these 64
	 * bits: nothing if there is no remainder, and more than half a
	 * unit if the remainder is more than the root.  It is never just
	 * a half.
	 */
	if (remh || rem.m64 > root.m64)
		dest->lowmant = 0x81;
	else
		dest->lowmant = rem.m64 != 0;

	return dest;
}

/*
 * FETOX follows setox of Motorola's floating-point package for the
 * 68040 (arch/m68k/fpsp040/setox.S; its notice is at the head of this
 * file) step by step, so that it computes the same digits: with
 * N = 64 x / log 2 to the nearest integer, N = 64 M + J and
 * R = x - N log 2 / 64,
 *
 *	e^x = 2^M * 2^(J/64) * e^R
 *
 * where 2^(J/64) comes from a table and e^R - 1 from a polynomial in R.
 * The package states an error below 0.85 units in the last place.
 *
 * The constants are the package's.  It has most of them in single or
 * double precision, and they are the same numbers here: 64/log 2 is
 * 0x42b8aa3b, L1 0xbc317218, A1 0x3f000000, A2 0x3fc55555 0x55554018,
 * A3 0x3fa55555 0x55554431, A4 0x3c088895 and A5 0x3ab60b70.
 */
static const struct fp_ext fp_exp_64byln2 =
	FPT_EXT(0, 0x4005, 0xb8aa3b00, 0x00000000);
/* L1 + L2 is -log 2 / 64 to 88 bits */
static const struct fp_ext fp_exp_l1 =
	FPT_EXT(1, 0x3ff8, 0xb1721800, 0x00000000);
static const struct fp_ext fp_exp_l2 =
	FPT_EXT(0, 0x3fdc, 0x82e30865, 0x4361c4c6);
static const struct fp_ext fp_exp_a1 =
	FPT_EXT(0, 0x3ffe, 0x80000000, 0x00000000);
static const struct fp_ext fp_exp_a2 =
	FPT_EXT(0, 0x3ffc, 0xaaaaaaaa, 0xaa00c000);
static const struct fp_ext fp_exp_a3 =
	FPT_EXT(0, 0x3ffa, 0xaaaaaaaa, 0xaa218800);
static const struct fp_ext fp_exp_a4 =
	FPT_EXT(0, 0x3ff8, 0x88889500, 0x00000000);
static const struct fp_ext fp_exp_a5 =
	FPT_EXT(0, 0x3ff5, 0xb60b7000, 0x00000000);

/*
 * e^x for a normalized number or zero.  With env the last operation is
 * rounded as the program asked, without env to nearest: see
 * fpt_last_add().  The hyperbolic instructions use this as well.
 */
void fp_etox(struct fp_ext *res, const struct fp_ext *x,
	     struct fpt_env *env)
{
	struct fp_ext n, r, s, p, q;
	unsigned int compact;
	int k, j, m, m1 = 0;

	/* below 2^-65: 1 + x, which the rounding mode decides */
	if (x->exp < 0x3fbe) {
		fpt_pow2(res, 0);
		fpt_last_add(res, x, env);
		return;
	}

	/* beyond 16480 log 2 the result is out of range */
	compact = fpt_compact(x);
	if (compact > 0x400cb27c) {
		if (x->sign)
			fpt_underflow(res, env);
		else
			fpt_overflow(res, 0, env);
		return;
	}

	/* N, as an integer and as a number */
	r = *x;
	fpt_mul(&r, &fp_exp_64byln2);
	k = fpt_to_int(&r);
	fpt_from_int(&n, k);
	j = k & 63;
	m = k >> 6;
	/*
	 * From 16380 log 2 on 2^M may be out of range where the result is
	 * not: scale by 2^M1 and 2^M with M1 + M as M was, M1 about half.
	 */
	if (compact >= 0x400cb167) {
		m1 = m >> 1;
		m -= m1;
	}

	/* R = (x + N * L1) + N * L2; the first sum is exact */
	r = n;
	fpt_mul(&r, &fp_exp_l1);
	q = n;
	fpt_mul(&q, &fp_exp_l2);
	fpt_add(&r, x);
	fpt_add(&r, &q);

	/*
	 * With S = R * R, e^R - 1 is
	 * [R + R * S * (A2 + S * A4)] + [S * (A1 + S * (A3 + S * A5))]
	 */
	s = r;
	fpt_mul(&s, &r);
	p = fp_exp_a5;
	fpt_mul(&p, &s);
	q = s;
	fpt_mul(&q, &fp_exp_a4);
	fpt_add(&p, &fp_exp_a3);
	fpt_add(&q, &fp_exp_a2);
	fpt_mul(&p, &s);
	fpt_mul(&q, &s);
	fpt_add(&p, &fp_exp_a1);
	fpt_mul(&q, &r);
	fpt_mul(&p, &s);
	fpt_add(&r, &q);
	fpt_add(&r, &p);

	/* 2^(J/64) * e^R = T + (T * (e^R - 1) + t) */
	fpt_mul(&r, &fpt_exptbl[j][0]);
	fpt_add(&r, &fpt_exptbl[j][1]);
	fpt_add(&r, &fpt_exptbl[j][0]);

	if (compact >= 0x400cb167) {
		fpt_pow2(&s, m1);
		fpt_mul(&r, &s);
	}

	/* the last operation: times 2^M */
	fpt_pow2(&s, m);
	fpt_last_mul(&r, &s, env);
	*res = r;
}

struct fp_ext *fp_fetox(struct fp_ext *dest, struct fp_ext *src)
{
	struct fpt_env env;
	struct fp_ext r;

	dprint(PINSTR, "fetox\n");

	if (fpt_special(FPT_FETOX, dest, src))
		return dest;

	fpt_enter(&env);
	fp_etox(&r, src, &env);

	return fpt_computed(dest, &r, &env);
}

/*
 * FETOXM1 follows setoxm1 of the same file.  From 1/4 to 70 log 2 it
 * reduces the operand as FETOX does and computes
 *
 *	e^x - 1 = 2^M * (2^(J/64) + 2^(J/64) * (e^R - 1) - 2^-M)
 *
 * in an order of additions that depends on M.  Below 1/4 it is a
 * polynomial in x, below 2^-65 x itself, and beyond 70 log 2 either e^x
 * or -1.  The package states an error below 0.85 units in the last
 * place.
 *
 * The constants of single and double precision in the package's form:
 *
 *	A2	0x3fc55555 0x55555555	B3	0x3fa55555 0x55555555
 *	A3	0x3fa55555 0x55554f5a	B4	0x3f811111 0x11111111
 *	A4	0x3f811111 0x11174385	B5	0x3f56c16c 0x16c170e2
 *	A5	0x3ab60b6a		B6	0x3f2a01a0 0x1a019df3
 *	A6	0x3950097b		B7	0x3efa01a0 0x19d7cb68
 *	B9	0x3493f281		B8	0x3ec71de3 0xa5774682
 *	B10	0x32d73220		B11	0x310f8290
 *	B12	0x2f30caa8
 *
 * A1 and B1 are 1/2, which FETOX has as its own A1.
 */
static const struct fp_ext fp_em1_a2 =
	FPT_EXT(0, 0x3ffc, 0xaaaaaaaa, 0xaaaaa800);
static const struct fp_ext fp_em1_a3 =
	FPT_EXT(0, 0x3ffa, 0xaaaaaaaa, 0xaa7ad000);
static const struct fp_ext fp_em1_a4 =
	FPT_EXT(0, 0x3ff8, 0x88888888, 0xba1c2800);
static const struct fp_ext fp_em1_a5 =
	FPT_EXT(0, 0x3ff5, 0xb60b6a00, 0x00000000);
static const struct fp_ext fp_em1_a6 =
	FPT_EXT(0, 0x3ff2, 0xd0097b00, 0x00000000);

static const struct fp_ext fp_em1_b2 =
	FPT_EXT(0, 0x3ffc, 0xaaaaaaaa, 0xaaaaaaab);
static const struct fp_ext fp_em1_b3 =
	FPT_EXT(0, 0x3ffa, 0xaaaaaaaa, 0xaaaaa800);
static const struct fp_ext fp_em1_b4 =
	FPT_EXT(0, 0x3ff8, 0x88888888, 0x88888800);
static const struct fp_ext fp_em1_b5 =
	FPT_EXT(0, 0x3ff5, 0xb60b60b6, 0x0b871000);
static const struct fp_ext fp_em1_b6 =
	FPT_EXT(0, 0x3ff2, 0xd00d00d0, 0x0cef9800);
static const struct fp_ext fp_em1_b7 =
	FPT_EXT(0, 0x3fef, 0xd00d00ce, 0xbe5b4000);
static const struct fp_ext fp_em1_b8 =
	FPT_EXT(0, 0x3fec, 0xb8ef1d2b, 0xba341000);
static const struct fp_ext fp_em1_b9 =
	FPT_EXT(0, 0x3fe9, 0x93f28100, 0x00000000);
static const struct fp_ext fp_em1_b10 =
	FPT_EXT(0, 0x3fe5, 0xd7322000, 0x00000000);
static const struct fp_ext fp_em1_b11 =
	FPT_EXT(0, 0x3fe2, 0x8f829000, 0x00000000);
static const struct fp_ext fp_em1_b12 =
	FPT_EXT(0, 0x3fde, 0xb0caa800, 0x00000000);

/* e^x - 1 for a normalized number below 1/4 in magnitude */
static void fp_etoxm1_small(struct fp_ext *res, const struct fp_ext *x,
			    struct fpt_env *env)
{
	struct fp_ext s, p, q;

	/*
	 * Below 2^-65: x less the tiny term 2^-16382, which the rounding
	 * mode decides.  Against the smallest operands that term is not
	 * tiny: for them the package scales x by 2^140 first, and the term
	 * is then lost in every rounding mode.
	 *
	 * Here the package's code is not followed.  Its description
	 * (setox.S, step 8.1) and the comment on its compare have these
	 * operands end at 2^-16312, but the constant of the compare,
	 * 0x00330000, is 2^-16332.  From there to 2^-16312 the code
	 * returns x - 2^-16382, which is up to 8192 units in the last
	 * place from e^x - 1.  The description's border is taken instead.
	 * Motorola's package for the 68060 and QEMU's FPU have the code's
	 * constant, so that this differs from all three in that range.
	 */
	if (x->exp < 0x3fbe) {
		*res = *x;
		fpt_pow2(&s, -16382);
		s.sign = 1;
		if (x->exp < 0x0047) {
			fpt_pow2(&p, 140);
			fpt_mul(res, &p);
			fpt_add(res, &s);
			fpt_pow2(&p, -140);
			fpt_last_mul(res, &p, env);
		} else {
			fpt_last_add(res, &s, env);
		}
		return;
	}

	/*
	 * With S = x * x, e^x - 1 is x + (S * B1 + Q) and Q is
	 * [x * S * (B2 + S * (B4 + S * (B6 + S * (B8 + S * (B10 +
	 *					       S * B12)))))] +
	 * [S * S * (B3 + S * (B5 + S * (B7 + S * (B9 + S * B11))))]
	 */
	s = *x;
	fpt_mul(&s, x);
	p = fp_em1_b12;
	fpt_mul(&p, &s);
	q = fp_em1_b11;
	fpt_add(&p, &fp_em1_b10);
	fpt_mul(&q, &s);
	fpt_mul(&p, &s);
	fpt_add(&q, &fp_em1_b9);
	fpt_add(&p, &fp_em1_b8);
	fpt_mul(&q, &s);
	fpt_mul(&p, &s);
	fpt_add(&q, &fp_em1_b7);
	fpt_add(&p, &fp_em1_b6);
	fpt_mul(&q, &s);
	fpt_mul(&p, &s);
	fpt_add(&q, &fp_em1_b5);
	fpt_add(&p, &fp_em1_b4);
	fpt_mul(&q, &s);
	fpt_mul(&p, &s);
	fpt_add(&q, &fp_em1_b3);
	fpt_add(&p, &fp_em1_b2);
	fpt_mul(&q, &s);
	fpt_mul(&p, &s);
	fpt_mul(&q, &s);
	fpt_mul(&p, x);
	fpt_mul(&s, &fp_exp_a1);
	fpt_add(&p, &q);
	fpt_add(&s, &p);

	fpt_last_add(&s, x, env);
	*res = s;
}

/*
 * e^x - 1 for a normalized number.  With env the last operation is
 * rounded as the program asked, without env to nearest: see
 * fpt_last_add().  The hyperbolic instructions use this as well.
 */
void fp_etoxm1(struct fp_ext *res, const struct fp_ext *x,
	       struct fpt_env *env)
{
	struct fp_ext n, r, s, p, q;
	const struct fp_ext *t;
	int k, m;

	if (x->exp < 0x3ffd) {
		fp_etoxm1_small(res, x, env);
		return;
	}

	/* beyond 70 log 2: e^x, or -1 with a tiny term */
	if (fpt_compact(x) > 0x4004c215) {
		if (!x->sign) {
			fp_etox(res, x, env);
			return;
		}
		fpt_pow2(res, 0);
		res->sign = 1;
		fpt_pow2(&s, -126);
		fpt_last_add(res, &s, env);
		return;
	}

	/* N, M, J and R as for e^x */
	r = *x;
	fpt_mul(&r, &fp_exp_64byln2);
	k = fpt_to_int(&r);
	fpt_from_int(&n, k);
	t = fpt_exptbl[k & 63];
	m = k >> 6;
	r = n;
	fpt_mul(&r, &fp_exp_l1);
	q = n;
	fpt_mul(&q, &fp_exp_l2);
	fpt_add(&r, x);
	fpt_add(&r, &q);

	/*
	 * With S = R * R, e^R - 1 is
	 * [R * S * (A2 + S * (A4 + S * A6))] +
	 * [R + S * (A1 + S * (A3 + S * A5))]
	 */
	s = r;
	fpt_mul(&s, &r);
	p = fp_em1_a6;
	fpt_mul(&p, &s);
	q = s;
	fpt_mul(&q, &fp_em1_a5);
	fpt_add(&p, &fp_em1_a4);
	fpt_add(&q, &fp_em1_a3);
	fpt_mul(&p, &s);
	fpt_mul(&q, &s);
	fpt_add(&p, &fp_em1_a2);
	fpt_add(&q, &fp_exp_a1);
	fpt_mul(&p, &s);
	fpt_mul(&s, &q);
	fpt_mul(&p, &r);
	fpt_add(&r, &s);
	fpt_add(&r, &p);

	/* p = T * (e^R - 1); then T + t + p - 2^-M, the small terms first */
	fpt_mul(&r, &t[0]);
	fpt_pow2(&s, -m);
	s.sign = 1;
	if (m > 63) {
		q = t[1];
		fpt_add(&q, &s);
		fpt_add(&r, &q);
		fpt_add(&r, &t[0]);
	} else if (m < -3) {
		fpt_add(&r, &t[1]);
		fpt_add(&r, &t[0]);
		fpt_add(&r, &s);
	} else {
		q = t[0];
		fpt_add(&r, &t[1]);
		fpt_add(&q, &s);
		fpt_add(&r, &q);
	}

	/* the last operation: times 2^M */
	fpt_pow2(&s, m);
	fpt_last_mul(&r, &s, env);
	*res = r;
}

struct fp_ext *fp_fetoxm1(struct fp_ext *dest, struct fp_ext *src)
{
	struct fpt_env env;
	struct fp_ext r;

	dprint(PINSTR, "fetoxm1\n");

	if (fpt_special(FPT_FETOXM1, dest, src))
		return dest;

	fpt_enter(&env);
	fp_etoxm1(&r, src, &env);

	return fpt_computed(dest, &r, &env);
}

/*
 * FTWOTOX and FTENTOX follow stwotox and stentox of the package
 * (arch/m68k/fpsp040/stwotox.S).  With N = 64 x, or 64 x log 10 / log 2,
 * to the nearest integer and N = 64 (M + M') + J, where M is about M',
 *
 *	2^x or 10^x = 2^M' * (2^M * 2^(J/64)) * e^R
 *
 * where 2^(J/64) comes from a table of this file's own and e^R - 1 from
 * a polynomial in R, the rest of x as a natural logarithm.  The package
 * states an error below 2 units in the last place.
 *
 * The constants of double precision in the package's form:
 *
 *	A1	0x3fe00000 0x00000000	64 log 10 / log 2
 *	A2	0x3fc55555 0x55554a54		0x406a934f 0x0979a371
 *	A3	0x3fa55555 0x55554cc1	the first part of log 2 / 64 log 10
 *	A4	0x3f811112 0x302c712c		0x3f734413 0x509f8000
 *	A5	0x3f56c16d 0x6f7bd0b2
 *
 * A1 is 1/2, which FETOX has as its own A1.
 */
static const struct fp_ext fp_log2 =
	FPT_EXT(0, 0x3ffe, 0xb17217f7, 0xd1cf79ac);
static const struct fp_ext fp_log10 =
	FPT_EXT(0, 0x4000, 0x935d8ddd, 0xaaa8ac17);
static const struct fp_ext fp_ten_l2ten64 =
	FPT_EXT(0, 0x4006, 0xd49a784b, 0xcd1b8800);
static const struct fp_ext fp_ten_l10two1 =
	FPT_EXT(0, 0x3ff7, 0x9a209a84, 0xfc000000);
static const struct fp_ext fp_ten_l10two2 =
	FPT_EXT(1, 0x3fcd, 0xc0219dc1, 0xda994fd2);
static const struct fp_ext fp_two_a2 =
	FPT_EXT(0, 0x3ffc, 0xaaaaaaaa, 0xaa52a000);
static const struct fp_ext fp_two_a3 =
	FPT_EXT(0, 0x3ffa, 0xaaaaaaaa, 0xaa660800);
static const struct fp_ext fp_two_a4 =
	FPT_EXT(0, 0x3ff8, 0x88889181, 0x63896000);
static const struct fp_ext fp_two_a5 =
	FPT_EXT(0, 0x3ff5, 0xb60b6b7b, 0xde859000);

/*
 * The operands that take no reduction: below 2^-70 the result is 1 + x,
 * which the rounding mode decides, and beyond the bound, which is 16480
 * as a power of two, it is out of range.
 */
static bool fp_twoten_edge(struct fp_ext *res, const struct fp_ext *x,
			   unsigned int bound, struct fpt_env *env)
{
	unsigned int compact = fpt_compact(x);

	if (compact < 0x3fb98000) {
		fpt_pow2(res, 0);
		fpt_last_add(res, x, env);
	} else if (compact <= bound) {
		return false;
	} else if (x->sign) {
		fpt_underflow(res, env);
	} else {
		fpt_overflow(res, 0, env);
	}

	return true;
}

/* The result from R, which r is, and N: expr of the package. */
static void fp_twoten_expr(struct fp_ext *r, int n, struct fpt_env *env)
{
	struct fp_ext s, p, q, fact1, fact2;
	int l = n >> 6, m = l >> 1;

	/* Fact1 + Fact2 = 2^M * 2^(J/64) */
	fact1 = fpt_exp2tbl[n & 63][0];
	fact1.exp += m;
	fact2 = fpt_exp2tbl[n & 63][1];
	fact2.exp += m;

	/*
	 * With S = R * R, e^R - 1 is
	 * [R + R * S * (A2 + S * A4)] + [S * (A1 + S * (A3 + S * A5))]
	 */
	s = *r;
	fpt_mul(&s, r);
	p = fp_two_a5;
	q = fp_two_a4;
	fpt_mul(&p, &s);
	fpt_mul(&q, &s);
	fpt_add(&p, &fp_two_a3);
	fpt_add(&q, &fp_two_a2);
	fpt_mul(&p, &s);
	fpt_mul(&q, &s);
	fpt_add(&p, &fp_exp_a1);
	fpt_mul(&q, r);
	fpt_mul(&p, &s);
	fpt_add(r, &q);
	fpt_add(r, &p);

	/* Fact1 + (Fact1 * (e^R - 1) + Fact2) */
	fpt_mul(r, &fact1);
	fpt_add(r, &fact2);
	fpt_add(r, &fact1);

	/* the last operation: times 2^M' */
	fpt_pow2(&s, l - m);
	fpt_last_mul(r, &s, env);
}

struct fp_ext *fp_ftwotox(struct fp_ext *dest, struct fp_ext *src)
{
	struct fpt_env env;
	struct fp_ext n, r;
	int k;

	dprint(PINSTR, "ftwotox\n");

	if (fpt_special(FPT_FTWOTOX, dest, src))
		return dest;

	fpt_enter(&env);
	if (!fp_twoten_edge(&r, src, 0x400d80c0, &env)) {
		/* R = (x - N / 64) * log 2 */
		n = *src;
		fpt_pow2(&r, 6);
		fpt_mul(&n, &r);
		k = fpt_to_int(&n);
		fpt_from_int(&n, k);
		fpt_pow2(&r, -6);
		fpt_mul(&n, &r);
		r = *src;
		fpt_sub(&r, &n);
		fpt_mul(&r, &fp_log2);
		fp_twoten_expr(&r, k, &env);
	}

	return fpt_computed(dest, &r, &env);
}

struct fp_ext *fp_ftentox(struct fp_ext *dest, struct fp_ext *src)
{
	struct fpt_env env;
	struct fp_ext n, r, q;
	int k;

	dprint(PINSTR, "ftentox\n");

	if (fpt_special(FPT_FTENTOX, dest, src))
		return dest;

	fpt_enter(&env);
	if (!fp_twoten_edge(&r, src, 0x400b9b07, &env)) {
		/*
		 * R = ((x - N * L1) - N * L2) * log 10, with
		 * L1 + L2 = log 2 / 64 log 10
		 */
		n = *src;
		fpt_mul(&n, &fp_ten_l2ten64);
		k = fpt_to_int(&n);
		fpt_from_int(&n, k);
		q = n;
		fpt_mul(&n, &fp_ten_l10two1);
		fpt_mul(&q, &fp_ten_l10two2);
		r = *src;
		fpt_sub(&r, &n);
		fpt_sub(&r, &q);
		fpt_mul(&r, &fp_log10);
		fp_twoten_expr(&r, k, &env);
	}

	return fpt_computed(dest, &r, &env);
}

/*
 * The logarithms follow slogn and slog2 of the package
 * (arch/m68k/fpsp040/slogn.S and slog2.S).  For x = 2^K * Y with Y from
 * 1 to 2, and F the first seven bits of Y with a one after them,
 *
 *	log x = K log 2 + log F + log(1 + U),	U = (Y - F) / F
 *
 * where 1/F and log F come from a table and log(1 + U) from a polynomial
 * in U.  Next to 1 it is a polynomial in U = 2 (x - 1) / (x + 1)
 * instead.  The package states an error below 2 units in the last place
 * for the natural logarithms and below 1.7 for those to base 10 and 2.
 * It does not keep these bounds: the sum from the table is right to
 * about a unit in the last place of 1/2 however small it is, and just
 * outside 15/16 to 17/16, where the other formula ends, the logarithm
 * is as small as 0.06.  Results more than ten units off were found
 * there, and more than twenty for the logarithm to base 10.
 *
 * The constants of double precision in the package's form:
 *
 *	A1	0xbfe00000 0x00000008	B1	0x3fb55555 0x55555555
 *	A2	0x3fd55555 0x555555a4	B2	0x3f899999 0x999995ec
 *	A3	0xbfcfffff 0xff6f7e97	B3	0x3f624924 0x928bccff
 *	A4	0x3fc99999 0x987d8730	B4	0x3f3c71c2 0xfe80c7e0
 *	A5	0xbfc555b5 0x848cb7db	B5	0x3f175496 0xadd7dad6
 *	A6	0x3fc2499a 0xb5e4040b
 */
static const struct fp_ext fp_log_a1 =
	FPT_EXT(1, 0x3ffe, 0x80000000, 0x00004000);
static const struct fp_ext fp_log_a2 =
	FPT_EXT(0, 0x3ffd, 0xaaaaaaaa, 0xaaad2000);
static const struct fp_ext fp_log_a3 =
	FPT_EXT(1, 0x3ffc, 0xfffffffb, 0x7bf4b800);
static const struct fp_ext fp_log_a4 =
	FPT_EXT(0, 0x3ffc, 0xccccccc3, 0xec398000);
static const struct fp_ext fp_log_a5 =
	FPT_EXT(1, 0x3ffc, 0xaaadac24, 0x65bed800);
static const struct fp_ext fp_log_a6 =
	FPT_EXT(0, 0x3ffc, 0x924cd5af, 0x20205800);

static const struct fp_ext fp_log_b1 =
	FPT_EXT(0, 0x3ffb, 0xaaaaaaaa, 0xaaaaa800);
static const struct fp_ext fp_log_b2 =
	FPT_EXT(0, 0x3ff8, 0xcccccccc, 0xccaf6000);
static const struct fp_ext fp_log_b3 =
	FPT_EXT(0, 0x3ff6, 0x92492494, 0x5e67f800);
static const struct fp_ext fp_log_b4 =
	FPT_EXT(0, 0x3ff3, 0xe38e17f4, 0x063f0000);
static const struct fp_ext fp_log_b5 =
	FPT_EXT(0, 0x3ff1, 0xbaa4b56e, 0xbed6b000);

/*
 * K log 2 + log F + log(1 + U) from Y - F, which u is, and the place j
 * of F in the table.  With V = U * U, log(1 + U) is
 * [U + V * (A1 + V * (A3 + V * A5))] + [U * V * (A2 + V * (A4 + V * A6))]
 */
static void fp_log_table(struct fp_ext *u, int k, unsigned int j,
			 struct fpt_env *env)
{
	struct fp_ext v, p, q, klog2;

	fpt_mul(u, &fpt_logtbl[j][0]);
	fpt_from_int(&klog2, k);
	fpt_mul(&klog2, &fp_log2);
	v = *u;
	fpt_mul(&v, u);
	p = v;
	q = v;
	fpt_mul(&p, &fp_log_a6);
	fpt_mul(&q, &fp_log_a5);
	fpt_add(&p, &fp_log_a4);
	fpt_add(&q, &fp_log_a3);
	fpt_mul(&p, &v);
	fpt_mul(&q, &v);
	fpt_add(&p, &fp_log_a2);
	fpt_add(&q, &fp_log_a1);
	fpt_mul(&p, &v);
	fpt_mul(&q, &v);
	fpt_mul(&p, u);
	fpt_add(u, &q);
	fpt_add(&p, &fpt_logtbl[j][1]);
	fpt_add(u, &p);

	fpt_last_add(u, &klog2, env);
}

/*
 * log x for a normalized positive x, which y is, away from 1: K is its
 * exponent, or what the caller knows it to be.
 */
static void fp_log_main(struct fp_ext *y, int k, struct fpt_env *env)
{
	struct fp_ext f;

	y->exp = 0x3fff;
	f = (struct fp_ext)FPT_EXT(0, 0x3fff,
		(y->mant.m32[0] & 0xfe000000) | 0x01000000, 0);
	fpt_sub(y, &f);
	fp_log_table(y, k, f.mant.m32[0] >> 25 & 63, env);
}

/*
 * log((2 + U) / (2 - U)) for U = num / den, into num.  With V = U * U
 * and W = V * V it is
 * U + U * V * ([B1 + W * (B3 + W * B5)] + [V * (B2 + W * B4)])
 */
static void fp_log_near1(struct fp_ext *num, const struct fp_ext *den,
			 struct fpt_env *env)
{
	struct fp_ext u, v, w, a, b;

	fpt_div(num, den);
	u = *num;
	v = u;
	fpt_mul(&v, &u);
	w = v;
	fpt_mul(&w, &v);
	a = fp_log_b5;
	b = fp_log_b4;
	fpt_mul(&a, &w);
	fpt_mul(&b, &w);
	fpt_add(&a, &fp_log_b3);
	fpt_add(&b, &fp_log_b2);
	fpt_mul(&w, &a);
	fpt_mul(&b, &v);
	fpt_add(&w, &fp_log_b1);
	fpt_mul(&v, &u);
	fpt_add(&w, &b);
	fpt_mul(&v, &w);

	fpt_last_add(&v, &u, env);
	*num = v;
}

/*
 * log x for a positive number, normalized or denormalized, that is not
 * 1.  With env the last operation is rounded as the program asked,
 * without env to nearest: see fpt_last_add().
 */
static void fp_logn(struct fp_ext *res, const struct fp_ext *x,
		    struct fpt_env *env)
{
	struct fp_ext one, den;
	unsigned int compact;
	int k;

	*res = *x;
	k = x->exp - 0x3fff;
	/* a denormalized number has a lower exponent than its field says */
	if ((long)res->mant.m32[0] >= 0)
		k -= fp_overnormalize(res);

	/* from about 15/16 to 17/16: U = 2 (x - 1) / (x + 1) */
	compact = fpt_compact(res);
	if (compact >= 0x3ffef07d && compact <= 0x3fff8841) {
		fpt_pow2(&one, 0);
		fpt_sub(res, &one);
		den = *x;
		fpt_add(&den, &one);
		fpt_add(res, res);
		fp_log_near1(res, &den, env);
		return;
	}

	fp_log_main(res, k, env);
}

struct fp_ext *fp_flogn(struct fp_ext *dest, struct fp_ext *src)
{
	struct fpt_env env;
	struct fp_ext r;

	dprint(PINSTR, "flogn\n");

	if (fpt_special(FPT_FLOGN, dest, src))
		return dest;

	fpt_enter(&env);
	fp_logn(&r, src, &env);

	return fpt_computed(dest, &r, &env);
}

struct fp_ext *fp_flognp1(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("flognp1\n");

	fp_monadic_check(dest, src);

	return dest;
}

struct fp_ext *fp_flog10(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("flog10\n");

	fp_monadic_check(dest, src);

	return dest;
}

struct fp_ext *fp_flog2(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("flog2\n");

	fp_monadic_check(dest, src);

	return dest;
}

struct fp_ext *fp_fgetexp(struct fp_ext *dest, struct fp_ext *src)
{
	int exp;

	dprint(PINSTR, "fgetexp\n");

	fp_monadic_check(dest, src);

	if (IS_INF(dest)) {
		fp_set_nan(dest);
		return dest;
	}
	if (IS_ZERO(dest))
		return dest;

	/* the exponent that a denormalized number has when normalized */
	exp = dest->exp;
	if ((long)dest->mant.m32[0] >= 0)
		exp -= fp_overnormalize(dest);

	fp_conv_long2ext(dest, exp - 0x3FFF);

	fp_normalize_ext(dest);

	return dest;
}

struct fp_ext *fp_fgetman(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "fgetman\n");

	fp_monadic_check(dest, src);

	if (IS_ZERO(dest))
		return dest;

	if (IS_INF(dest)) {
		fp_set_nan(dest);
		return dest;
	}

	/* the mantissa of a denormalized number lacks its highest bits */
	if ((long)dest->mant.m32[0] >= 0)
		fp_overnormalize(dest);

	dest->exp = 0x3FFF;

	return dest;
}

