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

#include "fp_emu.h"
#include "fp_log.h"
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

struct fp_ext *fp_fetox(struct fp_ext *dest, struct fp_ext *src)
{
	uprint("fetox\n");

	fp_monadic_check(dest, src);

	return dest;
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

