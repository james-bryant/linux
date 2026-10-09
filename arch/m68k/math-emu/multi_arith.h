/* SPDX-License-Identifier: GPL-2.0-or-later */
/* multi_arith.h: multi-precision integer arithmetic functions, needed
   to do extended-precision floating point.

   (c) 1998 David Huggins-Daines.

   Somewhat based on arch/alpha/math-emu/ieee-math.c, which is (c)
   David Mosberger-Tang.

 */

/* Note:

   These are not general multi-precision math routines.  Rather, they
   implement the subset of integer arithmetic that we need in order to
   multiply, divide, and normalize 128-bit unsigned mantissae.  */

#ifndef _MULTI_ARITH_H
#define _MULTI_ARITH_H

#include "fp_emu.h"

static inline void fp_denormalize(struct fp_ext *reg, unsigned int cnt)
{
	/*
	 * bits already shifted out, e.g. by a multiply, are below the new
	 * ones and must survive as the sticky bit for the rounding
	 */
	unsigned char sticky = reg->lowmant ? 1 : 0;

	reg->exp += cnt;

	switch (cnt) {
	case 0 ... 8:
		reg->lowmant = reg->mant.m32[1] << (8 - cnt);
		reg->mant.m32[1] = (reg->mant.m32[1] >> cnt) |
				   (reg->mant.m32[0] << (32 - cnt));
		reg->mant.m32[0] = reg->mant.m32[0] >> cnt;
		break;
	case 9 ... 32:
		reg->lowmant = reg->mant.m32[1] >> (cnt - 8);
		if (reg->mant.m32[1] << (40 - cnt))
			reg->lowmant |= 1;
		reg->mant.m32[1] = (reg->mant.m32[1] >> cnt) |
				   (reg->mant.m32[0] << (32 - cnt));
		reg->mant.m32[0] = reg->mant.m32[0] >> cnt;
		break;
	case 33 ... 39:
		asm volatile ("bfextu %1{%2,#8},%0" : "=d" (reg->lowmant)
			: "m" (reg->mant.m32[0]), "d" (64 - cnt));
		if (reg->mant.m32[1] << (40 - cnt))
			reg->lowmant |= 1;
		reg->mant.m32[1] = reg->mant.m32[0] >> (cnt - 32);
		reg->mant.m32[0] = 0;
		break;
	case 40 ... 71:
		reg->lowmant = reg->mant.m32[0] >> (cnt - 40);
		if ((reg->mant.m32[0] << (72 - cnt)) || reg->mant.m32[1])
			reg->lowmant |= 1;
		reg->mant.m32[1] = reg->mant.m32[0] >> (cnt - 32);
		reg->mant.m32[0] = 0;
		break;
	default:
		reg->lowmant = reg->mant.m32[0] || reg->mant.m32[1];
		reg->mant.m32[0] = 0;
		reg->mant.m32[1] = 0;
		break;
	}
	reg->lowmant |= sticky;
}

static inline int fp_overnormalize(struct fp_ext *reg)
{
	int shift;

	if (reg->mant.m32[0]) {
		asm ("bfffo %1{#0,#32},%0" : "=d" (shift) : "dm" (reg->mant.m32[0]));
		reg->mant.m32[0] = (reg->mant.m32[0] << shift) | (reg->mant.m32[1] >> (32 - shift));
		reg->mant.m32[1] = (reg->mant.m32[1] << shift);
	} else {
		asm ("bfffo %1{#0,#32},%0" : "=d" (shift) : "dm" (reg->mant.m32[1]));
		reg->mant.m32[0] = (reg->mant.m32[1] << shift);
		reg->mant.m32[1] = 0;
		shift += 32;
	}

	return shift;
}

/*
 * The multi-word add, subtract and shift below hand the carry from one
 * instruction to the next in the X condition-code bit.  Each is written
 * as a single asm statement: the compiler may neither reorder the
 * instructions nor place one of its own that writes X between them, both
 * of which it is free to do with the separate statements this once was.
 */
static inline int fp_addmant(struct fp_ext *dest, struct fp_ext *src)
{
	int carry;

	asm volatile ("add.b %4,%0\n\t"
		      "addx.l %5,%1\n\t"
		      "addx.l %6,%2\n\t"
		      "addx.l %3,%3"
		      : "=&d" (dest->lowmant), "=&d" (dest->mant.m32[1]),
			"=&d" (dest->mant.m32[0]), "=&d" (carry)
		      : "d" (src->lowmant), "d" (src->mant.m32[1]),
			"d" (src->mant.m32[0]), "0" (dest->lowmant),
			"1" (dest->mant.m32[1]), "2" (dest->mant.m32[0]),
			"3" (0));

	return carry;
}

static inline int fp_addcarry(struct fp_ext *reg)
{
	/* shift the carry in first: fp_set_ovrflw() reads the mantissa */
	reg->lowmant = (reg->mant.m32[1] << 7) | (reg->lowmant ? 1 : 0);
	reg->mant.m32[1] = (reg->mant.m32[1] >> 1) |
			   (reg->mant.m32[0] << 31);
	reg->mant.m32[0] = (reg->mant.m32[0] >> 1) | 0x80000000;
	if (++reg->exp == 0x7fff) {
		fp_set_ovrflw(reg);
		return 0;
	}

	return 1;
}

static inline void fp_submant(struct fp_ext *dest, struct fp_ext *src1,
			      struct fp_ext *src2)
{
	asm volatile ("sub.b %3,%0\n\t"
		      "subx.l %4,%1\n\t"
		      "subx.l %5,%2"
		      : "=&d" (dest->lowmant), "=&d" (dest->mant.m32[1]),
			"=&d" (dest->mant.m32[0])
		      : "d" (src2->lowmant), "d" (src2->mant.m32[1]),
			"d" (src2->mant.m32[0]), "0" (src1->lowmant),
			"1" (src1->mant.m32[1]), "2" (src1->mant.m32[0]));
}

#define fp_mul64(desth, destl, src1, src2) ({				\
	asm ("mulu.l %2,%1:%0" : "=d" (destl), "=d" (desth)		\
		: "dm" (src1), "0" (src2));				\
})
#define fp_div64(quot, rem, srch, srcl, div)				\
	asm ("divu.l %2,%1:%0" : "=d" (quot), "=d" (rem)		\
		: "dm" (div), "1" (srch), "0" (srcl))
#define fp_add64(dest1, dest2, src1, src2) ({				\
	asm ("add.l %3,%0\n\t"						\
	     "addx.l %2,%1"						\
		: "=&d" (dest2), "=&d" (dest1)				\
		: "d" (src1), "dm" (src2), "0" (dest2), "1" (dest1));	\
})
#define fp_addx96(dest, src) ({						\
	asm volatile ("add.l %3,%0\n\t"					\
		      "addx.l %4,%1\n\t"				\
		      "addx.l %5,%2"					\
		: "=&d" ((dest)->m32[2]), "=&d" ((dest)->m32[1]),	\
		  "=&d" ((dest)->m32[0])				\
		: "dm" ((src).m32[1]), "d" ((src).m32[0]), "d" (0UL),	\
		  "0" ((dest)->m32[2]), "1" ((dest)->m32[1]),		\
		  "2" ((dest)->m32[0]));				\
})
#define fp_sub64(dest, src) ({						\
	asm ("sub.l %3,%0\n\t"						\
	     "subx.l %2,%1"						\
		: "=&d" ((dest).m32[1]), "=&d" ((dest).m32[0])		\
		: "d" ((src).m32[0]), "dm" ((src).m32[1]),		\
		  "0" ((dest).m32[1]), "1" ((dest).m32[0]));		\
})
#define fp_sub96c(dest, srch, srcm, srcl) ({				\
	char carry;							\
	asm ("sub.l %4,%0\n\t"						\
	     "subx.l %5,%1\n\t"						\
	     "subx.l %6,%2\n\t"						\
	     "scs %3"							\
		: "=&d" ((dest).m32[2]), "=&d" ((dest).m32[1]),		\
		  "=&d" ((dest).m32[0]), "=d" (carry)			\
		: "dm" (srcl), "d" (srcm), "d" (srch),			\
		  "0" ((dest).m32[2]), "1" ((dest).m32[1]),		\
		  "2" ((dest).m32[0]));					\
	carry;								\
})

static inline void fp_multiplymant(union fp_mant128 *dest, struct fp_ext *src1,
				   struct fp_ext *src2)
{
	union fp_mant64 temp;

	fp_mul64(dest->m32[0], dest->m32[1], src1->mant.m32[0], src2->mant.m32[0]);
	fp_mul64(dest->m32[2], dest->m32[3], src1->mant.m32[1], src2->mant.m32[1]);

	fp_mul64(temp.m32[0], temp.m32[1], src1->mant.m32[0], src2->mant.m32[1]);
	fp_addx96(dest, temp);

	fp_mul64(temp.m32[0], temp.m32[1], src1->mant.m32[1], src2->mant.m32[0]);
	fp_addx96(dest, temp);
}

static inline void fp_dividemant(union fp_mant128 *dest, struct fp_ext *src,
				 struct fp_ext *div)
{
	union fp_mant128 tmp;
	union fp_mant64 tmp64;
	unsigned long *mantp = dest->m32;
	unsigned long fix, rem, first, dummy;
	int i;

	/* the algorithm below requires dest to be smaller than div,
	   but both have the high bit set */
	if (src->mant.m64 >= div->mant.m64) {
		fp_sub64(src->mant, div->mant);
		*mantp = 1;
	} else
		*mantp = 0;
	mantp++;

	/* basic idea behind this algorithm: we can't divide two 64bit numbers
	   (AB/CD) directly, but we can calculate AB/C0, but this means this
	   quotient is off by C0/CD, so we have to multiply the first result
	   to fix the result, after that we have nearly the correct result
	   and only a few corrections are needed. */

	/* C0/CD can be precalculated, but it's an 64bit division again, but
	   we can make it a bit easier, by dividing first through C so we get
	   10/1D and now only a single shift and the value fits into 32bit. */
	fix = 0x80000000;
	dummy = div->mant.m32[1] / div->mant.m32[0] + 1;
	dummy = (dummy >> 1) | fix;
	fp_div64(fix, dummy, fix, 0, dummy);
	fix--;

	for (i = 0; i < 3; i++, mantp++) {
		if (src->mant.m32[0] == div->mant.m32[0]) {
			fp_div64(first, rem, 0, src->mant.m32[1], div->mant.m32[0]);

			fp_mul64(*mantp, dummy, first, fix);
			*mantp += fix;
		} else {
			fp_div64(first, rem, src->mant.m32[0], src->mant.m32[1], div->mant.m32[0]);

			fp_mul64(*mantp, dummy, first, fix);
		}

		fp_mul64(tmp.m32[0], tmp.m32[1], div->mant.m32[0], first - *mantp);
		fp_add64(tmp.m32[0], tmp.m32[1], 0, rem);
		tmp.m32[2] = 0;

		fp_mul64(tmp64.m32[0], tmp64.m32[1], *mantp, div->mant.m32[1]);
		fp_sub96c(tmp, 0, tmp64.m32[0], tmp64.m32[1]);

		src->mant.m32[0] = tmp.m32[1];
		src->mant.m32[1] = tmp.m32[2];

		while (!fp_sub96c(tmp, 0, div->mant.m32[0], div->mant.m32[1])) {
			src->mant.m32[0] = tmp.m32[1];
			src->mant.m32[1] = tmp.m32[2];
			*mantp += 1;
		}
	}

	/*
	 * src now holds the remainder.  If it is not zero, the quotient
	 * goes on below its last bit here: set that bit, which
	 * fp_putmant128() takes into the sticky bit and nowhere else, so
	 * that the rounding sees an inexact quotient.
	 */
	if (src->mant.m64)
		dest->m32[3] |= 1;
}

static inline void fp_putmant128(struct fp_ext *dest, union fp_mant128 *src,
				 int shift)
{
	unsigned long tmp, dummy;

	switch (shift) {
	case 0:
		dest->mant.m64 = src->m64[0];
		dest->lowmant = src->m32[2] >> 24;
		if (src->m32[3] || (src->m32[2] << 8))
			dest->lowmant |= 1;
		break;
	case 1:
		asm volatile ("lsl.l #1,%0\n\t"
			      "roxl.l #1,%1\n\t"
			      "roxl.l #1,%2"
			: "=d" (tmp), "=d" (dest->mant.m32[1]),
			  "=d" (dest->mant.m32[0])
			: "0" (src->m32[2]), "1" (src->m32[1]),
			  "2" (src->m32[0]));
		dest->lowmant = tmp >> 24;
		if (src->m32[3] || (tmp << 8))
			dest->lowmant |= 1;
		break;
	case 31:
		asm volatile ("lsr.l #1,%3\n\t"
			      "roxr.l #1,%0\n\t"
			      "roxr.l #1,%1\n\t"
			      "roxr.l #1,%2"
			: "=d" (dest->mant.m32[0]), "=d" (dest->mant.m32[1]),
			  "=d" (tmp), "=d" (dummy)
			: "0" (src->m32[1]), "1" (src->m32[2]),
			  "2" (src->m32[3]), "3" (src->m32[0]));
		dest->lowmant = tmp >> 24;
		if (src->m32[3] << 7)
			dest->lowmant |= 1;
		break;
	case 32:
		dest->mant.m32[0] = src->m32[1];
		dest->mant.m32[1] = src->m32[2];
		dest->lowmant = src->m32[3] >> 24;
		if (src->m32[3] << 8)
			dest->lowmant |= 1;
		break;
	}
}

#endif	/* _MULTI_ARITH_H */
