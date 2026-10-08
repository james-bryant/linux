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
 * with the functions and the constants above them that they use, are the
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
 * fp_fsqrt(), fp_fetoxm1(), fp_ftwotox(), fp_ftentox(), fp_flogn(),
 * fp_flognp1(), fp_flog10(), fp_flog2(), fp_fgetexp() and fp_fgetman()
 * are not taken from the package.
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

struct fp_ext *fp_fetoxm1(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("fetoxm1\n");

	fp_monadic_check(dest, src);

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

struct fp_ext *fp_ftwotox(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("ftwotox\n");

	fp_monadic_check(dest, src);

	return dest;
}

struct fp_ext *fp_ftentox(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("ftentox\n");

	fp_monadic_check(dest, src);

	return dest;
}

struct fp_ext *fp_flogn(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("flogn\n");

	fp_monadic_check(dest, src);

	return dest;
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

