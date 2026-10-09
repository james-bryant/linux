// SPDX-License-Identifier: GPL-2.0-or-later
/*

   fp_arith.c: floating-point math routines for the Linux-m68k
   floating point emulator.

   Copyright (c) 1998-1999 David Huggins-Daines.

   Somewhat based on the AlphaLinux floating point emulator, by David
   Mosberger-Tang.

 */

#include "fp_emu.h"
#include "multi_arith.h"
#include "fp_arith.h"

const struct fp_ext fp_QNaN =
{
	.exp = 0x7fff,
	.mant = { .m64 = ~0 }
};

const struct fp_ext fp_Inf =
{
	.exp = 0x7fff,
};

/* let's start with the easy ones */

struct fp_ext *fp_fabs(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "fabs\n");

	fp_monadic_check(dest, src);

	dest->sign = 0;

	return dest;
}

struct fp_ext *fp_fneg(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "fneg\n");

	fp_monadic_check(dest, src);

	dest->sign = !dest->sign;

	return dest;
}

/* Now, the slightly harder ones */

/*
 * fp_fadd: Implements the kernel of the FADD, FSADD, FDADD, FSUB and
 * FDSUB instructions.
 */

struct fp_ext *fp_fadd(struct fp_ext *dest, struct fp_ext *src)
{
	int diff;

	dprint(PINSTR, "fadd\n");

	fp_dyadic_check(dest, src);

	if (IS_INF(dest)) {
		/* infinity - infinity == NaN */
		if (IS_INF(src) && (src->sign != dest->sign))
			fp_set_nan(dest);
		return dest;
	}
	if (IS_INF(src)) {
		fp_copy_ext(dest, src);
		return dest;
	}

	if (IS_ZERO(dest)) {
		if (IS_ZERO(src)) {
			if (src->sign != dest->sign) {
				if (FPDATA->rnd == FPCR_ROUND_RM)
					dest->sign = 1;
				else
					dest->sign = 0;
			}
		} else
			fp_copy_ext(dest, src);
		return dest;
	}

	dest->lowmant = src->lowmant = 0;

	if ((diff = dest->exp - src->exp) > 0)
		__fp_denormalize(src, diff);
	else if ((diff = -diff) > 0)
		__fp_denormalize(dest, diff);

	if (dest->sign == src->sign) {
		if (fp_addmant(dest, src))
			if (!fp_addcarry(dest))
				return dest;
	} else {
		if (dest->mant.m64 < src->mant.m64) {
			fp_submant(dest, src, dest);
			dest->sign = !dest->sign;
		} else
			fp_submant(dest, dest, src);
		/* numbers that cancel: the zero of the rounding mode */
		if (!dest->mant.m64 && !dest->lowmant)
			dest->sign = FPDATA->rnd == FPCR_ROUND_RM;
	}

	return dest;
}

/* fp_fsub: Implements the kernel of the FSUB, FSSUB, and FDSUB
   instructions.

   Remember that the arguments are in assembler-syntax order! */

struct fp_ext *fp_fsub(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "fsub ");

	/* a NaN is returned as it is, also one in the source */
	fp_dyadic_check(dest, src);

	src->sign = !src->sign;
	return fp_fadd(dest, src);
}


/*
 * fp_fcmp: Implements the FCMP instruction, which sets the condition
 * codes as for dest - src without computing the difference: equal
 * infinities are equal, and their difference is not a number.
 *
 * The caller derives the condition codes from the value returned: a
 * NaN operand without its sign, a zero for equal operands, which is
 * negative only for a negative zero or infinity in dest, or else a
 * number with the sign of the difference.  It is never an infinity.
 */

struct fp_ext *fp_fcmp(struct fp_ext *dest, struct fp_ext *src)
{
	struct fp_ext *res = &FPDATA->temp[1];
	int cmp;

	dprint(PINSTR, "fcmp ");

	fp_copy_ext(res, dest);

	/* unordered: FCMP does not give N the sign of a NaN */
	if (!fp_normalize_ext(res)) {
		fp_check_snan(src);
		res->sign = 0;
		return res;
	}
	if (!fp_normalize_ext(src)) {
		fp_copy_ext(res, src);
		res->sign = 0;
		return res;
	}

	/* +0 and -0 are equal */
	if (IS_ZERO(src))
		src->sign = res->sign;

	/* both are normalized, and the mantissa of an infinity is ignored */
	if (res->sign != src->sign)
		cmp = 1;
	else if (res->exp != src->exp)
		cmp = res->exp > src->exp ? 1 : -1;
	else if (IS_INF(res) || res->mant.m64 == src->mant.m64)
		cmp = 0;
	else
		cmp = res->mant.m64 > src->mant.m64 ? 1 : -1;

	if (cmp) {
		if (cmp < 0)
			res->sign = !res->sign;
		res->exp = 0x3fff;
		res->mant.m64 = 1ULL << 63;
	} else {
		if (!IS_INF(res) && !IS_ZERO(res))
			res->sign = 0;
		res->exp = 0;
		res->mant.m64 = 0;
	}
	res->lowmant = 0;

	return res;
}

struct fp_ext *fp_ftst(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "ftst\n");

	(void)dest;

	return src;
}

struct fp_ext *fp_fmul(struct fp_ext *dest, struct fp_ext *src)
{
	union fp_mant128 temp;
	int exp;

	dprint(PINSTR, "fmul\n");

	fp_dyadic_check(dest, src);

	/* calculate the correct sign now, as it's necessary for infinities */
	dest->sign = src->sign ^ dest->sign;

	/* Handle infinities */
	if (IS_INF(dest)) {
		if (IS_ZERO(src))
			fp_set_nan(dest);
		return dest;
	}
	if (IS_INF(src)) {
		if (IS_ZERO(dest)) {
			fp_set_nan(dest);
		} else {
			/* the infinity gets the sign of the product */
			src->sign = dest->sign;
			fp_copy_ext(dest, src);
		}
		return dest;
	}

	/* Of course, as we all know, zero * anything = zero.  You may
	   not have known that it might be a positive or negative
	   zero... */
	if (IS_ZERO(dest) || IS_ZERO(src)) {
		dest->exp = 0;
		dest->mant.m64 = 0;
		dest->lowmant = 0;

		return dest;
	}

	exp = dest->exp + src->exp - 0x3ffe;

	/* shift up the mantissa for denormalized numbers,
	   so that the highest bit is set, this makes the
	   shift of the result below easier */
	if ((long)dest->mant.m32[0] >= 0)
		exp -= fp_overnormalize(dest);
	if ((long)src->mant.m32[0] >= 0)
		exp -= fp_overnormalize(src);

	/* now, do a 64-bit multiply with expansion */
	fp_multiplymant(&temp, dest, src);

	/* normalize it back to 64 bits and stuff it back into the
	   destination struct */
	if ((long)temp.m32[0] > 0) {
		exp--;
		fp_putmant128(dest, &temp, 1);
	} else
		fp_putmant128(dest, &temp, 0);

	if (exp >= 0x7fff) {
		fp_set_ovrflw(dest);
		return dest;
	}
	dest->exp = exp;
	if (exp < 0) {
		fp_set_sr(FPSR_EXC_UNFL);
		fp_denormalize(dest, -exp);
	}

	return dest;
}

/* fp_fdiv: Implements the "kernel" of the FDIV, FSDIV, FDDIV and
   FSGLDIV instructions.

   Note that the order of the operands is counter-intuitive: instead
   of src / dest, the result is actually dest / src. */

struct fp_ext *fp_fdiv(struct fp_ext *dest, struct fp_ext *src)
{
	union fp_mant128 temp;
	int exp;

	dprint(PINSTR, "fdiv\n");

	fp_dyadic_check(dest, src);

	/* calculate the correct sign now, as it's necessary for infinities */
	dest->sign = src->sign ^ dest->sign;

	/* Handle infinities */
	if (IS_INF(dest)) {
		/* infinity / infinity = NaN (quiet, as always) */
		if (IS_INF(src))
			fp_set_nan(dest);
		/* infinity / anything else = infinity (with appropriate sign) */
		return dest;
	}
	if (IS_INF(src)) {
		/* anything / infinity = zero (with appropriate sign) */
		dest->exp = 0;
		dest->mant.m64 = 0;
		dest->lowmant = 0;

		return dest;
	}

	/* zeroes */
	if (IS_ZERO(dest)) {
		/* zero / zero = NaN */
		if (IS_ZERO(src))
			fp_set_nan(dest);
		/* zero / anything else = zero */
		return dest;
	}
	if (IS_ZERO(src)) {
		/* anything / zero = infinity (with appropriate sign) */
		fp_set_sr(FPSR_EXC_DZ);
		dest->exp = 0x7fff;
		dest->mant.m64 = 0;

		return dest;
	}

	exp = dest->exp - src->exp + 0x3fff;

	/* shift up the mantissa for denormalized numbers,
	   so that the highest bit is set, this makes lots
	   of things below easier */
	if ((long)dest->mant.m32[0] >= 0)
		exp -= fp_overnormalize(dest);
	if ((long)src->mant.m32[0] >= 0)
		exp += fp_overnormalize(src);

	/* now, do the 64-bit divide */
	fp_dividemant(&temp, dest, src);

	/* normalize it back to 64 bits and stuff it back into the
	   destination struct */
	if (!temp.m32[0]) {
		exp--;
		fp_putmant128(dest, &temp, 32);
	} else
		fp_putmant128(dest, &temp, 31);

	if (exp >= 0x7fff) {
		fp_set_ovrflw(dest);
		return dest;
	}
	dest->exp = exp;
	if (exp < 0) {
		fp_set_sr(FPSR_EXC_UNFL);
		fp_denormalize(dest, -exp);
	}

	return dest;
}

struct fp_ext *fp_fsglmul(struct fp_ext *dest, struct fp_ext *src)
{
	int exp;

	dprint(PINSTR, "fsglmul\n");

	fp_dyadic_check(dest, src);

	/* what an earlier instruction left there is not for fp_set_ovrflw() */
	dest->lowmant = 0;

	/* calculate the correct sign now, as it's necessary for infinities */
	dest->sign = src->sign ^ dest->sign;

	/* Handle infinities */
	if (IS_INF(dest)) {
		if (IS_ZERO(src))
			fp_set_nan(dest);
		return dest;
	}
	if (IS_INF(src)) {
		if (IS_ZERO(dest)) {
			fp_set_nan(dest);
		} else {
			/* the infinity gets the sign of the product */
			src->sign = dest->sign;
			fp_copy_ext(dest, src);
		}
		return dest;
	}

	/* Of course, as we all know, zero * anything = zero.  You may
	   not have known that it might be a positive or negative
	   zero... */
	if (IS_ZERO(dest) || IS_ZERO(src)) {
		dest->exp = 0;
		dest->mant.m64 = 0;
		dest->lowmant = 0;

		return dest;
	}

	exp = dest->exp + src->exp - 0x3ffe;

	/*
	 * normalize the operands, so that a denormalized one keeps its
	 * significant bits when it is truncated to single precision
	 */
	if ((long)dest->mant.m32[0] >= 0)
		exp -= fp_overnormalize(dest);
	if ((long)src->mant.m32[0] >= 0)
		exp -= fp_overnormalize(src);

	/* do a 32-bit multiply */
	fp_mul64(dest->mant.m32[0], dest->mant.m32[1],
		 dest->mant.m32[0] & 0xffffff00,
		 src->mant.m32[0] & 0xffffff00);

	/*
	 * the product of two normalized numbers is in [1, 4): normalize it
	 * to [2, 4) before the exponent is tested for overflow
	 */
	if ((long)dest->mant.m32[0] >= 0) {
		exp--;
		dest->mant.m64 <<= 1;
	}

	if (exp >= 0x7fff) {
		/* the bits that single precision drops make it inexact */
		if ((dest->mant.m32[0] & 0xff) | dest->mant.m32[1])
			dest->lowmant = 1;
		fp_set_ovrflw(dest);
		return dest;
	}
	dest->exp = exp;
	if (exp < 0) {
		fp_set_sr(FPSR_EXC_UNFL);
		fp_denormalize(dest, -exp);
	}

	return dest;
}

struct fp_ext *fp_fsgldiv(struct fp_ext *dest, struct fp_ext *src)
{
	int exp;
	unsigned long quot, rem;

	dprint(PINSTR, "fsgldiv\n");

	fp_dyadic_check(dest, src);

	/* what an earlier instruction left there is not for fp_set_ovrflw() */
	dest->lowmant = 0;

	/* calculate the correct sign now, as it's necessary for infinities */
	dest->sign = src->sign ^ dest->sign;

	/* Handle infinities */
	if (IS_INF(dest)) {
		/* infinity / infinity = NaN (quiet, as always) */
		if (IS_INF(src))
			fp_set_nan(dest);
		/* infinity / anything else = infinity (with approprate sign) */
		return dest;
	}
	if (IS_INF(src)) {
		/* anything / infinity = zero (with appropriate sign) */
		dest->exp = 0;
		dest->mant.m64 = 0;
		dest->lowmant = 0;

		return dest;
	}

	/* zeroes */
	if (IS_ZERO(dest)) {
		/* zero / zero = NaN */
		if (IS_ZERO(src))
			fp_set_nan(dest);
		/* zero / anything else = zero */
		return dest;
	}
	if (IS_ZERO(src)) {
		/* anything / zero = infinity (with appropriate sign) */
		fp_set_sr(FPSR_EXC_DZ);
		dest->exp = 0x7fff;
		dest->mant.m64 = 0;

		return dest;
	}

	exp = dest->exp - src->exp + 0x3fff;

	/*
	 * normalize the operands first: a denormalized divisor whose top
	 * 24 mantissa bits are zero would otherwise be truncated to zero
	 * below, and fp_div64() would divide by zero and trap in the kernel
	 */
	if ((long)dest->mant.m32[0] >= 0)
		exp -= fp_overnormalize(dest);
	if ((long)src->mant.m32[0] >= 0)
		exp += fp_overnormalize(src);

	/*
	 * truncate both mantissas to single precision: the top 24 bits of
	 * the high longword, nothing of the low one.  The low longword must
	 * go too, not only be left out of the divide: fp_sub64() below
	 * subtracts the whole 64-bit mantissas, and a non-zero low longword
	 * would borrow into the 24-bit high longwords and give a quotient
	 * that is wrong, not merely inaccurate.  PRM 5-112: "the extraneous
	 * mantissa bits are truncated prior to the division".
	 */
	dest->mant.m32[0] &= 0xffffff00;
	dest->mant.m32[1] = 0;
	src->mant.m32[0] &= 0xffffff00;
	src->mant.m32[1] = 0;

	/* do the 32-bit divide */
	if (dest->mant.m32[0] >= src->mant.m32[0]) {
		fp_sub64(dest->mant, src->mant);
		fp_div64(quot, rem, dest->mant.m32[0], 0, src->mant.m32[0]);
		dest->mant.m32[0] = 0x80000000 | (quot >> 1);
		dest->mant.m32[1] = (quot & 1) | rem;	/* only for rounding */
	} else {
		fp_div64(quot, rem, dest->mant.m32[0], 0, src->mant.m32[0]);
		dest->mant.m32[0] = quot;
		dest->mant.m32[1] = rem;		/* only for rounding */
		exp--;
	}

	if (exp >= 0x7fff) {
		/* the bits that single precision drops make it inexact */
		if ((dest->mant.m32[0] & 0xff) | dest->mant.m32[1])
			dest->lowmant = 1;
		fp_set_ovrflw(dest);
		return dest;
	}
	dest->exp = exp;
	if (exp < 0) {
		fp_set_sr(FPSR_EXC_UNFL);
		fp_denormalize(dest, -exp);
	}

	return dest;
}

/* fp_roundint: Internal rounding function for use by several of these
   emulated instructions.

   This one rounds off the fractional part using the rounding mode
   specified. */

static void fp_roundint(struct fp_ext *dest, int mode)
{
	union fp_mant64 oldmant;
	unsigned long mask;

	if (!fp_normalize_ext(dest))
		return;

	/* infinities and zeroes */
	if (IS_INF(dest) || IS_ZERO(dest))
		return;

	/* first truncate the lower bits */
	oldmant = dest->mant;
	switch (dest->exp) {
	case 0 ... 0x3ffe:
		dest->mant.m64 = 0;
		break;
	case 0x3fff ... 0x401e:
		dest->mant.m32[0] &= 0xffffffffU << (0x401e - dest->exp);
		dest->mant.m32[1] = 0;
		if (oldmant.m64 == dest->mant.m64)
			return;
		break;
	case 0x401f ... 0x403e:
		dest->mant.m32[1] &= 0xffffffffU << (0x403e - dest->exp);
		if (oldmant.m32[1] == dest->mant.m32[1])
			return;
		break;
	default:
		return;
	}
	fp_set_sr(FPSR_EXC_INEX2);

	/* We might want to normalize upwards here... however, since
	   we know that this is only called on the output of fp_fdiv,
	   or with the input to fp_fint or fp_fintrz, and the inputs
	   to all these functions are either normal or denormalized
	   (no subnormals allowed!), there's really no need.

	   In the case of fp_fdiv, observe that 0x80000000 / 0xffff =
	   0xffff8000, and the same holds for 128-bit / 64-bit. (i.e. the
	   smallest possible normal dividend and the largest possible normal
	   divisor will still produce a normal quotient, therefore, (normal
	   << 64) / normal is normal in all cases) */

	switch (mode) {
	case FPCR_ROUND_RN:
		switch (dest->exp) {
		case 0 ... 0x3ffd:
			return;
		case 0x3ffe:
			/* As noted above, the input is always normal, so the
			   guard bit (bit 63) is always set.  therefore, the
			   only case in which we will NOT round to 1.0 is when
			   the input is exactly 0.5. */
			if (oldmant.m64 == (1ULL << 63))
				return;
			break;
		case 0x3fff ... 0x401d:
			mask = 1 << (0x401d - dest->exp);
			if (!(oldmant.m32[0] & mask))
				return;
			if (oldmant.m32[0] & (mask << 1))
				break;
			/* in two steps: the count reaches 32 at 0x401d */
			if (!(oldmant.m32[0] << (dest->exp - 0x3ffe) << 1) &&
			    !oldmant.m32[1])
				return;
			break;
		case 0x401e:
			if (!(oldmant.m32[1] & 0x80000000))
				return;
			if (oldmant.m32[0] & 1)
				break;
			if (!(oldmant.m32[1] << 1))
				return;
			break;
		case 0x401f ... 0x403d:
			mask = 1 << (0x403d - dest->exp);
			if (!(oldmant.m32[1] & mask))
				return;
			if (oldmant.m32[1] & (mask << 1))
				break;
			/* and here at 0x403d */
			if (!(oldmant.m32[1] << (dest->exp - 0x401e) << 1))
				return;
			break;
		default:
			return;
		}
		break;
	case FPCR_ROUND_RZ:
		return;
	default:
		if (dest->sign ^ (mode - FPCR_ROUND_RM))
			break;
		return;
	}

	switch (dest->exp) {
	case 0 ... 0x3ffe:
		dest->exp = 0x3fff;
		dest->mant.m64 = 1ULL << 63;
		break;
	case 0x3fff ... 0x401e:
		mask = 1UL << (0x401e - dest->exp);
		if (dest->mant.m32[0] += mask)
			break;
		dest->mant.m32[0] = 0x80000000;
		dest->exp++;
		break;
	case 0x401f ... 0x403e:
		mask = 1UL << (0x403e - dest->exp);
		if (dest->mant.m32[1] += mask)
			break;
		if (dest->mant.m32[0] += 1)
                        break;
		dest->mant.m32[0] = 0x80000000;
                dest->exp++;
		break;
	}
}

/* modrem_kernel: Implementation of the FREM and FMOD instructions
   (which are exactly the same, except for the rounding used on the
   intermediate value) */

static struct fp_ext *modrem_kernel(struct fp_ext *dest, struct fp_ext *src,
				    int mode)
{
	unsigned long long rem, div;
	unsigned int quot = 0;
	int exp, last, sign, carry = 0;

	fp_dyadic_check(dest, src);

	/* Infinities and zeros */
	if (IS_INF(dest) || IS_ZERO(src)) {
		fp_set_nan(dest);
		return dest;
	}
	sign = dest->sign ^ src->sign;
	if (IS_ZERO(dest) || IS_INF(src)) {
		/* the quotient is zero */
		fp_set_quotient(sign << 7);
		return dest;
	}

	/* shift up the mantissa of denormalized numbers */
	exp = dest->exp;
	if ((long)dest->mant.m32[0] >= 0)
		exp -= fp_overnormalize(dest);
	last = src->exp;
	if ((long)src->mant.m32[0] >= 0)
		last -= fp_overnormalize(src);
	rem = dest->mant.m64;
	div = src->mant.m64;

	/*
	 * Divide the mantissas as integers, bit by bit: one bit of the
	 * quotient for every exponent from that of the dividend down to
	 * that of the divisor.  The remainder is exact, and of the
	 * quotient only the low bits are kept.  FREM divides one bit
	 * further: the bit that is worth one half.
	 */
	if (mode == FPCR_ROUND_RN)
		last--;
	if (exp >= last) {
		for (;;) {
			quot <<= 1;
			if (carry || rem >= div) {
				rem -= div;
				quot |= 1;
			}
			if (exp == last)
				break;
			exp--;
			carry = rem >> 63;
			rem <<= 1;
			/*
			 * the loop runs up to 32831 times for far-apart
			 * operands; let other tasks run meanwhile
			 */
			if (!(exp & 0x3ff))
				cond_resched();
		}
		if (mode == FPCR_ROUND_RN) {
			/*
			 * With that bit the next multiple of the divisor
			 * is as near or nearer.  It is the one to take
			 * unless it is as near and odd, and the remainder
			 * is what is missing to it.
			 */
			carry = quot & 1;
			quot >>= 1;
			if (carry) {
				if (rem || (quot & 1)) {
					quot++;
					dest->sign = !dest->sign;
				}
				rem = div - rem;
			}
		}
	}

	dest->mant.m64 = rem;
	dest->exp = exp;
	/* a remainder of denormalized numbers; no bit is lost */
	if (exp < 0)
		fp_denormalize(dest, -exp);

	/* set the quotient byte */
	fp_set_quotient((quot & 0x7f) | (sign << 7));
	return dest;
}

/* fp_fmod: Implements the kernel of the FMOD instruction.

   Again, the argument order is backwards.  The result, as defined in
   the Motorola manuals, is:

   fmod(src,dest) = (dest - (src * floor(dest / src))) */

struct fp_ext *fp_fmod(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "fmod\n");
	return modrem_kernel(dest, src, FPCR_ROUND_RZ);
}

/* fp_frem: Implements the kernel of the FREM instruction.

   frem(src,dest) = (dest - (src * round(dest / src)))
 */

struct fp_ext *fp_frem(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "frem\n");
	return modrem_kernel(dest, src, FPCR_ROUND_RN);
}

struct fp_ext *fp_fint(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "fint\n");

	fp_copy_ext(dest, src);

	fp_roundint(dest, FPDATA->rnd);

	return dest;
}

struct fp_ext *fp_fintrz(struct fp_ext *dest, struct fp_ext *src)
{
	dprint(PINSTR, "fintrz\n");

	fp_copy_ext(dest, src);

	fp_roundint(dest, FPCR_ROUND_RZ);

	return dest;
}

struct fp_ext *fp_fscale(struct fp_ext *dest, struct fp_ext *src)
{
	int scale, oldround;
	unsigned int oldsr;

	dprint(PINSTR, "fscale\n");

	fp_dyadic_check(dest, src);

	/* what an earlier instruction left there is not for fp_set_ovrflw() */
	dest->lowmant = 0;

	/* Infinities */
	if (IS_INF(src)) {
		fp_set_nan(dest);
		return dest;
	}
	if (IS_INF(dest))
		return dest;

	/* zeroes */
	if (IS_ZERO(src) || IS_ZERO(dest))
		return dest;

	if (src->exp >= 0x400d) {
		/* 2^14 or more: always an overflow or an underflow */
		scale = src->sign ? -0x10000 : 0x10000;
	} else {
		/*
		 * src must be rounded with round to zero,
		 * and that is not inexact.
		 */
		oldround = FPDATA->rnd;
		oldsr = FPDATA->fpsr;
		FPDATA->rnd = FPCR_ROUND_RZ;
		scale = fp_conv_ext2long(src);
		FPDATA->rnd = oldround;
		FPDATA->fpsr = oldsr;
	}

	/* new exponent */
	scale += dest->exp;

	if (scale >= 0x7fff) {
		fp_set_ovrflw(dest);
		return dest;
	}
	dest->exp = scale;
	if (scale < 0) {
		fp_set_sr(FPSR_EXC_UNFL);
		fp_denormalize(dest, -scale);
	}

	return dest;
}

