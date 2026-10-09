/* SPDX-License-Identifier: GPL-2.0-or-later */
/*
 * fp_trans.h: what the transcendental instructions of the Linux-m68k
 * floating point emulator share.  See fp_trans.c.
 */

#ifndef _FP_TRANS_H
#define _FP_TRANS_H

#include "fp_emu.h"

/* The transcendental instructions: the rows of the special operands. */
enum fpt_insn {
	FPT_FACOS,
	FPT_FASIN,
	FPT_FATAN,
	FPT_FATANH,
	FPT_FCOS,
	FPT_FCOSH,
	FPT_FETOX,
	FPT_FETOXM1,
	FPT_FLOG10,
	FPT_FLOG2,
	FPT_FLOGN,
	FPT_FLOGNP1,
	FPT_FSIN,
	FPT_FSINCOS,
	FPT_FSINH,
	FPT_FTAN,
	FPT_FTANH,
	FPT_FTENTOX,
	FPT_FTWOTOX,
	FPT_INSNS
};

/*
 * A constant: its sign, its exponent and the halves of its mantissa.
 * The constants of single and double precision of Motorola's package
 * are written as the same numbers in extended precision, which is what
 * an operand of those formats is in an operation of the coprocessor.
 */
#define FPT_EXT(s, e, hi, lo)						\
	{ .sign = (s), .exp = (e), .mant.m32 = { (hi), (lo) } }

/* fp_tables.c: the tables of the package, two numbers an entry */
extern const struct fp_ext fpt_exptbl[64][2];
extern const struct fp_ext fpt_exp2tbl[64][2];
extern const struct fp_ext fpt_logtbl[64][2];
extern const struct fp_ext fpt_pitbl[65][2];

/*
 * What fpt_enter() takes from the program, until the last operation of
 * the algorithm or what ends it gives it back.
 */
struct fpt_env {
	unsigned int fpsr;
	unsigned short prec;
	unsigned short rnd;
	bool entered;
};

bool fpt_special(enum fpt_insn insn, struct fp_ext *dest, struct fp_ext *src);
void fpt_enter(struct fpt_env *env);
struct fp_ext *fpt_computed(struct fp_ext *dest, const struct fp_ext *res,
			    struct fpt_env *env);
struct fp_ext *fpt_operr(struct fp_ext *dest, struct fpt_env *env);

void fpt_add(struct fp_ext *dest, const struct fp_ext *src);
void fpt_sub(struct fp_ext *dest, const struct fp_ext *src);
void fpt_mul(struct fp_ext *dest, const struct fp_ext *src);
void fpt_div(struct fp_ext *dest, const struct fp_ext *src);
void fpt_sqrt(struct fp_ext *dest, const struct fp_ext *src);

void fpt_last_add(struct fp_ext *dest, const struct fp_ext *src,
		  struct fpt_env *env);
void fpt_last_mul(struct fp_ext *dest, const struct fp_ext *src,
		  struct fpt_env *env);
void fpt_last_div(struct fp_ext *dest, const struct fp_ext *src,
		  struct fpt_env *env);
void fpt_overflow(struct fp_ext *res, int sign, struct fpt_env *env);
void fpt_underflow(struct fp_ext *res, struct fpt_env *env);

/* fp_log.c: the algorithms that other instructions use */
void fp_etox(struct fp_ext *res, const struct fp_ext *x,
	     struct fpt_env *env);
void fp_etoxm1(struct fp_ext *res, const struct fp_ext *x,
	       struct fpt_env *env);
void fp_lognp1(struct fp_ext *res, const struct fp_ext *x,
	       struct fpt_env *env);

int fpt_to_int(const struct fp_ext *src);
void fpt_from_int(struct fp_ext *dest, int val);

/* 2^exp */
static inline void fpt_pow2(struct fp_ext *dest, int exp)
{
	dest->lowmant = 0;
	dest->sign = 0;
	dest->exp = 0x3fff + exp;
	dest->mant.m32[0] = 0x80000000;
	dest->mant.m32[1] = 0;
}

/*
 * The magnitude of a normalized number in the form that Motorola's
 * package compares with its thresholds: the exponent and the first 16
 * bits of the mantissa.
 */
static inline unsigned int fpt_compact(const struct fp_ext *reg)
{
	return reg->exp << 16 | reg->mant.m32[0] >> 16;
}

#endif /* _FP_TRANS_H */
