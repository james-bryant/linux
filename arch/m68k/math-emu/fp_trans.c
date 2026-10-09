// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * fp_trans.c: what the transcendental instructions of the Linux-m68k
 * floating point emulator share.
 *
 * The instructions follow the algorithms of Motorola's floating-point
 * package for the 68040 (arch/m68k/fpsp040).  The package computes in
 * extended precision and rounds to nearest whatever the program's FPCR
 * says, and only the last operation of a function is rounded as the
 * program asked.  An instruction here does the same:
 *
 *	if (fpt_special(FPT_FXXX, dest, src))
 *		return dest;
 *	fpt_enter(&env);
 *	... fpt_add(), fpt_mul() and the others, on local variables ...
 *	fpt_last_add(&res, &last, &env);
 *	return fpt_computed(dest, &res, &env);
 *
 * The last operation leaves the nearest rounding and is not rounded
 * here.  Some algorithms of the package use others: the hyperbolic
 * sine takes an exponential.  The one that is used gets no fpt_env,
 * and its last operation is then an operation like the others.
 *
 * Only the layer puts the program's rounding mode, precision and FPSR
 * back: in the last operation, and in what ends an algorithm for one
 * that comes there without a last operation.  That is fpt_computed(),
 * or fpt_operr() for an algorithm that cannot go on.  No path of an
 * algorithm can leave them changed for the final rounding or for the
 * next instruction.
 *
 * What an instruction returns for a zero, an infinity or a denormalized
 * number, and which exception bits a computed result gets, is not part
 * of an algorithm: the first is the table fpt_table[] and the second is
 * fpt_flags().  The manuals of the MC68881 and MC68882 leave some of it
 * open, and the package decided those cases here.
 */

/*
 * This file contains a modified version of parts of Motorola's
 * floating-point package for the 68040 (arch/m68k/fpsp040): the table
 * fpt_table[] and the results fpt_results[] follow what the package's
 * routines for special operands return (do_func.S, tbldo.S and the
 * function's own file) wherever the manuals leave a result open, and
 * fpt_flags() sets INEX2 by the package's rule.  The package comes with
 * this notice (arch/m68k/fpsp040/README):
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
 * The rest of the file, the operations and the state that fpt_enter()
 * keeps for the last one, is not taken from the package.
 */

#include <linux/bitops.h>

#include "fp_emu.h"
#include "fp_arith.h"
#include "fp_log.h"
#include "fp_trans.h"

/* The classes of operands that an instruction can have a rule for. */
enum {
	FPT_ZERO,		/* a zero */
	FPT_DENORM,		/* a denormalized number */
	FPT_LT1,		/* a normalized number below 1 in magnitude */
	FPT_EQ1,		/* 1 or -1 */
	FPT_GT1,		/* a number above 1 in magnitude */
	FPT_INF,		/* an infinity */
	FPT_CLASSES
};

/* a class and the sign of the operand: the columns of fpt_table[] */
#define FPT_P(class)	(2 * (class))
#define FPT_N(class)	(2 * (class) + 1)

/* What the instruction returns for an operand of a class. */
enum {
	FPT_COMPUTE,		/* what its algorithm computes */
	FPT_OPERAND,		/* the operand */
	FPT_R_ZERO,		/* fpt_results[] from here on */
	FPT_R_ONE,
	FPT_R_ONE_UP,
	FPT_R_ONE_DOWN,
	FPT_R_INF,
	FPT_R_NAN,
	FPT_R_NEG = 0x80	/* flag: the result with a minus sign */
};

/*
 * A result may have bits below its mantissa, like the result of an
 * operation: it is then rounded to the program's precision in the
 * program's rounding mode.
 */
static const struct fp_ext fpt_results[] = {
	[FPT_R_ZERO] = {
		.exp = 0,
	},
	[FPT_R_ONE] = {
		.exp = 0x3fff,
		.mant.m32 = { 0x80000000, 0 },
	},
	/* 1 and a tiny term */
	[FPT_R_ONE_UP] = {
		.lowmant = 0x01,
		.exp = 0x3fff,
		.mant.m32 = { 0x80000000, 0 },
	},
	/* 1 less a tiny term */
	[FPT_R_ONE_DOWN] = {
		.lowmant = 0xff,
		.exp = 0x3ffe,
		.mant.m32 = { 0xffffffff, 0xffffffff },
	},
	[FPT_R_INF] = {
		.exp = 0x7fff,
	},
	/* the NaN of an operand error */
	[FPT_R_NAN] = {
		.exp = 0x7fff,
		.mant.m32 = { 0xffffffff, 0xffffffff },
	},
};

/* the exception status byte of FPSR */
#define FPT_INEX2	(1 << (FPSR_EXC_INEX2 - 8))
#define FPT_DZ		(1 << (FPSR_EXC_DZ - 8))
#define FPT_OPERR	(1 << (FPSR_EXC_OPERR - 8))

/*
 * The special operands: for every instruction and every class of
 * operand, of either sign, the result and the exception bits that come
 * with it.  The condition codes and the accrued exception bits follow
 * from these as for every other instruction.  A class without an entry
 * is computed.  A NaN is no class: every instruction returns a NaN
 * operand.
 *
 * Where the M68000 Family Programmer's Reference Manual and the
 * MC68881/MC68882 User's Manual give the result, it is theirs.  They
 * say nothing of denormalized operands; those are as in Motorola's
 * package for the 68040.  No entry was checked on a coprocessor.
 *
 * UNFL is in no entry.  A denormalized operand that is returned has it
 * with INEX2 in the package, and here it comes from fp_finalrounding(),
 * which sets UNFL for every result that is denormalized or underflows
 * in the program's precision.  An entry could not take it away.
 */
static const struct {
	unsigned char result;
	unsigned char exc;
} fpt_table[FPT_INSNS][2 * FPT_CLASSES] = {
	/*
	 * The package (setoxd) adds a tiny term of the operand's sign to 1
	 * for a denormalized operand, so that the rounding mode decides.
	 */
	[FPT_FETOX] = {
		[FPT_P(FPT_ZERO)]	= { FPT_R_ONE },
		[FPT_N(FPT_ZERO)]	= { FPT_R_ONE },
		[FPT_P(FPT_DENORM)]	= { FPT_R_ONE_UP, FPT_INEX2 },
		[FPT_N(FPT_DENORM)]	= { FPT_R_ONE_DOWN, FPT_INEX2 },
		[FPT_P(FPT_INF)]	= { FPT_R_INF },
		[FPT_N(FPT_INF)]	= { FPT_R_ZERO },
	},
	/*
	 * e^x - 1 of a denormalized number is that number as far as the
	 * precision shows, and the package (t_extdnrm) reports it as an
	 * underflow.
	 */
	[FPT_FETOXM1] = {
		[FPT_P(FPT_ZERO)]	= { FPT_OPERAND },
		[FPT_N(FPT_ZERO)]	= { FPT_OPERAND },
		[FPT_P(FPT_DENORM)]	= { FPT_OPERAND, FPT_INEX2 },
		[FPT_N(FPT_DENORM)]	= { FPT_OPERAND, FPT_INEX2 },
		[FPT_P(FPT_INF)]	= { FPT_R_INF },
		[FPT_N(FPT_INF)]	= { FPT_R_ONE | FPT_R_NEG },
	},
	/* the logarithms to base 10 and 2 as FLOGN */
	[FPT_FLOG10] = {
		[FPT_P(FPT_ZERO)]	= { FPT_R_INF | FPT_R_NEG, FPT_DZ },
		[FPT_N(FPT_ZERO)]	= { FPT_R_INF | FPT_R_NEG, FPT_DZ },
		[FPT_N(FPT_DENORM)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_N(FPT_LT1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_P(FPT_EQ1)]	= { FPT_R_ZERO },
		[FPT_N(FPT_EQ1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_N(FPT_GT1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_P(FPT_INF)]	= { FPT_R_INF },
		[FPT_N(FPT_INF)]	= { FPT_R_NAN, FPT_OPERR },
	},
	[FPT_FLOG2] = {
		[FPT_P(FPT_ZERO)]	= { FPT_R_INF | FPT_R_NEG, FPT_DZ },
		[FPT_N(FPT_ZERO)]	= { FPT_R_INF | FPT_R_NEG, FPT_DZ },
		[FPT_N(FPT_DENORM)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_N(FPT_LT1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_P(FPT_EQ1)]	= { FPT_R_ZERO },
		[FPT_N(FPT_EQ1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_N(FPT_GT1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_P(FPT_INF)]	= { FPT_R_INF },
		[FPT_N(FPT_INF)]	= { FPT_R_NAN, FPT_OPERR },
	},
	/*
	 * The logarithm of a zero of either sign is minus infinity with DZ,
	 * and that of a negative number is an operand error.  The logarithm
	 * of 1 is +0 and is exact: the package (sslogn, sslog10, sslog2)
	 * does not compute it.
	 */
	[FPT_FLOGN] = {
		[FPT_P(FPT_ZERO)]	= { FPT_R_INF | FPT_R_NEG, FPT_DZ },
		[FPT_N(FPT_ZERO)]	= { FPT_R_INF | FPT_R_NEG, FPT_DZ },
		[FPT_N(FPT_DENORM)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_N(FPT_LT1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_P(FPT_EQ1)]	= { FPT_R_ZERO },
		[FPT_N(FPT_EQ1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_N(FPT_GT1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_P(FPT_INF)]	= { FPT_R_INF },
		[FPT_N(FPT_INF)]	= { FPT_R_NAN, FPT_OPERR },
	},
	/*
	 * A denormalized operand as for FETOXM1.  For -1 the manuals have
	 * DZ with a NaN where they describe the instruction and with minus
	 * infinity where they describe the exception: the package (sslognp1)
	 * returns minus infinity.
	 */
	[FPT_FLOGNP1] = {
		[FPT_P(FPT_ZERO)]	= { FPT_OPERAND },
		[FPT_N(FPT_ZERO)]	= { FPT_OPERAND },
		[FPT_P(FPT_DENORM)]	= { FPT_OPERAND, FPT_INEX2 },
		[FPT_N(FPT_DENORM)]	= { FPT_OPERAND, FPT_INEX2 },
		[FPT_N(FPT_EQ1)]	= { FPT_R_INF | FPT_R_NEG, FPT_DZ },
		[FPT_N(FPT_GT1)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_P(FPT_INF)]	= { FPT_R_INF },
		[FPT_N(FPT_INF)]	= { FPT_R_NAN, FPT_OPERR },
	},
	/* a denormalized operand as for FETOXM1 */
	[FPT_FSIN] = {
		[FPT_P(FPT_ZERO)]	= { FPT_OPERAND },
		[FPT_N(FPT_ZERO)]	= { FPT_OPERAND },
		[FPT_P(FPT_DENORM)]	= { FPT_OPERAND, FPT_INEX2 },
		[FPT_N(FPT_DENORM)]	= { FPT_OPERAND, FPT_INEX2 },
		[FPT_P(FPT_INF)]	= { FPT_R_NAN, FPT_OPERR },
		[FPT_N(FPT_INF)]	= { FPT_R_NAN, FPT_OPERR },
	},
	/* a denormalized operand as for FETOX (stentoxd, stwotoxd) */
	[FPT_FTENTOX] = {
		[FPT_P(FPT_ZERO)]	= { FPT_R_ONE },
		[FPT_N(FPT_ZERO)]	= { FPT_R_ONE },
		[FPT_P(FPT_DENORM)]	= { FPT_R_ONE_UP, FPT_INEX2 },
		[FPT_N(FPT_DENORM)]	= { FPT_R_ONE_DOWN, FPT_INEX2 },
		[FPT_P(FPT_INF)]	= { FPT_R_INF },
		[FPT_N(FPT_INF)]	= { FPT_R_ZERO },
	},
	[FPT_FTWOTOX] = {
		[FPT_P(FPT_ZERO)]	= { FPT_R_ONE },
		[FPT_N(FPT_ZERO)]	= { FPT_R_ONE },
		[FPT_P(FPT_DENORM)]	= { FPT_R_ONE_UP, FPT_INEX2 },
		[FPT_N(FPT_DENORM)]	= { FPT_R_ONE_DOWN, FPT_INEX2 },
		[FPT_P(FPT_INF)]	= { FPT_R_INF },
		[FPT_N(FPT_INF)]	= { FPT_R_ZERO },
	},
};

/*
 * The exception bits of a computed result, beyond those that its last
 * operation and its rounding set.  The package sets INEX2 for every
 * computed result, also for one that happens to be exact.  The manuals
 * only say that INEX2 may be set then.
 */
static void fpt_flags(void)
{
	fp_set_sr(FPSR_EXC_INEX2);
}

/*
 * Normalize the result of an operation.  The bits below the mantissa
 * move up with it; nothing is rounded.  A number too small to normalize
 * stays denormalized, with an exponent of zero.
 */
static void fpt_normalize(struct fp_ext *reg)
{
	unsigned long hi = reg->mant.m32[0], lo = reg->mant.m32[1];
	unsigned long low = reg->lowmant;
	unsigned int shift;

	if (IS_INF(reg) || (long)hi < 0)
		return;
	if (!(hi | lo | low)) {
		reg->exp = 0;
		return;
	}

	if (hi)
		shift = 32 - fls(hi);
	else if (lo)
		shift = 64 - fls(lo);
	else
		shift = 72 - fls(low);
	if (shift > reg->exp)
		shift = reg->exp;
	reg->exp -= shift;

	for (; shift >= 8; shift -= 8) {
		hi = hi << 8 | lo >> 24;
		lo = lo << 8 | low;
		low = 0;
	}
	if (shift) {
		hi = hi << shift | lo >> (32 - shift);
		lo = lo << shift | low >> (8 - shift);
		low <<= shift;
	}
	reg->mant.m32[0] = hi;
	reg->mant.m32[1] = lo;
	reg->lowmant = low;
}

/*
 * Round the result of an operation to the 64 bits of the mantissa, to
 * the nearest and in a tie to the even one.  The emulator's own
 * rounding of an operand does this only with
 * CONFIG_M68KFPU_EMU_EXTRAPREC, and in the program's rounding mode.
 */
static void fpt_round(struct fp_ext *reg)
{
	unsigned char low;

	fpt_normalize(reg);
	low = reg->lowmant;
	reg->lowmant = 0;
	if (IS_INF(reg) || !(low & 0x80))
		return;
	if (!(low & 0x7f) && !(reg->mant.m32[1] & 1))
		return;

	if (++reg->mant.m32[1] || ++reg->mant.m32[0])
		return;
	reg->mant.m32[0] = 0x80000000;
	if (++reg->exp == 0x7fff)
		reg->mant.m32[0] = 0;
}

/*
 * The operations of an algorithm: dest = dest <op> src in extended
 * precision, rounded to nearest, after fpt_enter().
 * The emulator's functions modify both operands, hence the copy of the
 * source, which may then be a constant.
 */
void fpt_add(struct fp_ext *dest, const struct fp_ext *src)
{
	struct fp_ext tmp = *src;

	fp_fadd(dest, &tmp);
	fpt_round(dest);
}

void fpt_sub(struct fp_ext *dest, const struct fp_ext *src)
{
	struct fp_ext tmp = *src;

	fp_fsub(dest, &tmp);
	fpt_round(dest);
}

void fpt_mul(struct fp_ext *dest, const struct fp_ext *src)
{
	struct fp_ext tmp = *src;

	fp_fmul(dest, &tmp);
	fpt_round(dest);
}

void fpt_div(struct fp_ext *dest, const struct fp_ext *src)
{
	struct fp_ext tmp = *src;

	fp_fdiv(dest, &tmp);
	fpt_round(dest);
}

/* dest = sqrt(src) */
void fpt_sqrt(struct fp_ext *dest, const struct fp_ext *src)
{
	struct fp_ext tmp = *src;

	fp_fsqrt(dest, &tmp);
	fpt_round(dest);
}

/*
 * The integer nearest to a number, the even one in a tie, after
 * fpt_enter().  The number must be less than 2^31 in magnitude.
 */
int fpt_to_int(const struct fp_ext *src)
{
	struct fp_ext tmp = *src;

	return fp_conv_ext2long(&tmp);
}

void fpt_from_int(struct fp_ext *dest, int val)
{
	fp_conv_long2ext(dest, val);
	fp_normalize_ext(dest);
}

/*
 * Compute in extended precision and round to nearest from here on, as
 * the package does.  The emulator's functions set exception bits in
 * FPSR as they go and read the rounding mode and precision there, so
 * the program's are kept aside until the last operation.
 */
void fpt_enter(struct fpt_env *env)
{
	env->fpsr = FPDATA->fpsr;
	env->prec = FPDATA->prec;
	env->rnd = FPDATA->rnd;
	env->entered = true;
	FPDATA->prec = FPCR_PRECISION_X;
	FPDATA->rnd = FPCR_ROUND_RN;
}

/*
 * Back to the program's rounding mode and precision.  The exception
 * bits that the operations since fpt_enter() have set are dropped: only
 * the last operation's count.
 *
 * An algorithm does not call this.  Its last operation does,
 * fpt_last_add() or one of the others, and fpt_restore() does for an
 * algorithm that ended without one.
 */
static void fpt_leave(struct fpt_env *env)
{
	FPDATA->fpsr = env->fpsr;
	FPDATA->prec = env->prec;
	FPDATA->rnd = env->rnd;
	env->entered = false;
}

/*
 * This is where the layer makes sure of the program's state: every
 * function that ends an algorithm or rounds a result of it comes here
 * first.  env is what the algorithm entered with, or NULL if it never
 * did.  If the algorithm has not left, because the path that it took
 * has no last operation or because it forgot one, the program's
 * rounding mode, precision and FPSR are put back here, before
 * fp_finalrounding() and the next instruction read them.
 */
static void fpt_restore(struct fpt_env *env)
{
	if (env && env->entered)
		fpt_leave(env);
}

/*
 * The last operation of an algorithm, as the program rounds: it puts
 * the program's rounding mode and precision back first, and leaves its
 * result for fpt_computed() as the emulator's functions leave it.
 *
 * Without env the algorithm is used by another one that goes on
 * computing with the result: then this is an operation like the
 * others.
 */
void fpt_last_add(struct fp_ext *dest, const struct fp_ext *src,
		  struct fpt_env *env)
{
	struct fp_ext tmp = *src;

	if (env)
		fpt_leave(env);
	fp_fadd(dest, &tmp);
	if (!env)
		fpt_round(dest);
}

void fpt_last_mul(struct fp_ext *dest, const struct fp_ext *src,
		  struct fpt_env *env)
{
	struct fp_ext tmp = *src;

	if (env)
		fpt_leave(env);
	fp_fmul(dest, &tmp);
	if (!env)
		fpt_round(dest);
}

void fpt_last_div(struct fp_ext *dest, const struct fp_ext *src,
		  struct fpt_env *env)
{
	struct fp_ext tmp = *src;

	if (env)
		fpt_leave(env);
	fp_fdiv(dest, &tmp);
	if (!env)
		fpt_round(dest);
}

/*
 * A result beyond the largest number, of this sign: the overflow of
 * the program's rounding mode and precision, which the square of a
 * large number is.  It is the last operation of an algorithm.
 */
void fpt_overflow(struct fp_ext *res, int sign, struct fpt_env *env)
{
	struct fp_ext tmp;

	fpt_pow2(res, 16383);
	tmp = *res;
	tmp.sign = sign;
	fpt_last_mul(res, &tmp, env);
}

/* A positive result below the smallest number: the underflow, likewise. */
void fpt_underflow(struct fp_ext *res, struct fpt_env *env)
{
	struct fp_ext tmp;

	fpt_pow2(res, -16382);
	tmp = *res;
	fpt_last_mul(res, &tmp, env);
}

/*
 * Hand a result to fp_finalrounding(), which rounds it to the program's
 * precision in the program's mode and sets INEX2, UNFL, OVFL, the
 * condition codes and the accrued bits.  The result comes as the last
 * operation left it, with the bits below the mantissa in the low
 * mantissa byte.
 *
 * It is normalized here, because without
 * CONFIG_M68KFPU_EMU_EXTRAPREC fp_finalrounding() would lose those bits
 * when it does that itself.
 *
 * For single and double precision fp_finalrounding() rounds twice,
 * first to extended precision.  To make the second rounding the only
 * one, the first is taken from it: the bits below the mantissa only
 * tell the second rounding that there is something below its guard
 * bit, and the last bit of the mantissa can tell it as well.
 */
static void fpt_store(struct fp_ext *dest, const struct fp_ext *res)
{
	*dest = *res;
	fpt_normalize(dest);
	if (FPDATA->prec != FPCR_PRECISION_X && dest->lowmant &&
	    !IS_INF(dest)) {
		dest->mant.m32[1] |= 1;
		dest->lowmant = 0;
	}
}

/*
 * The result of an algorithm: the value of its last operation, which
 * was done in the program's rounding mode and has not been rounded.
 * env is what the algorithm entered with, or NULL if it never did.
 */
struct fp_ext *fpt_computed(struct fp_ext *dest, const struct fp_ext *res,
			    struct fpt_env *env)
{
	fpt_restore(env);
	fpt_store(dest, res);
	fpt_flags();

	return dest;
}

/*
 * The result of an algorithm that cannot go on: that of an operand
 * error, the NaN with OPERR that an instruction returns for an operand
 * outside its domain.
 */
struct fp_ext *fpt_operr(struct fp_ext *dest, struct fpt_env *env)
{
	fpt_restore(env);
	fpt_store(dest, &fpt_results[FPT_R_NAN]);
	fp_set_sr(FPSR_EXC_OPERR);

	return dest;
}

/*
 * Normalize the operand, and return its result if the instruction has
 * a rule for it: true if dest is the result.  If not, the operand is a
 * normalized number or, for an instruction that computes these, a
 * denormalized one.
 */
bool fpt_special(enum fpt_insn insn, struct fp_ext *dest, struct fp_ext *src)
{
	unsigned int class, result;

	/* a NaN is returned; fp_finalrounding() tests it for SNAN */
	if (!fp_normalize_ext(src)) {
		*dest = *src;
		return true;
	}
	/* an operand has nothing below its mantissa */
	src->lowmant = 0;

	if (IS_INF(src))
		class = FPT_INF;
	else if (!src->mant.m64)
		class = FPT_ZERO;
	else if ((long)src->mant.m32[0] >= 0)
		class = FPT_DENORM;
	else if (src->exp < 0x3fff)
		class = FPT_LT1;
	else if (src->exp == 0x3fff && src->mant.m32[0] == 0x80000000 &&
		 !src->mant.m32[1])
		class = FPT_EQ1;
	else
		class = FPT_GT1;
	class = src->sign ? FPT_N(class) : FPT_P(class);

	result = fpt_table[insn][class].result;
	switch (result & ~FPT_R_NEG) {
	case FPT_COMPUTE:
		return false;
	case FPT_OPERAND:
		fpt_store(dest, src);
		break;
	default:
		fpt_store(dest, &fpt_results[result & ~FPT_R_NEG]);
		break;
	}
	if (result & FPT_R_NEG)
		dest->sign = 1;
	FPDATA->fpsr |= fpt_table[insn][class].exc << 8;

	return true;
}
