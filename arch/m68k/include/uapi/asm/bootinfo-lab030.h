/* SPDX-License-Identifier: GPL-2.0 WITH Linux-syscall-note */
/*
 * asm/bootinfo-lab030.h -- m68030-lab specific boot information definitions
 */

#ifndef _UAPI_ASM_M68K_BOOTINFO_LAB030_H
#define _UAPI_ASM_M68K_BOOTINFO_LAB030_H

/*
 * Machine type (BI_MACHTYPE).  This number is private and provisional: it
 * has not been allocated in <asm/bootinfo.h>.  It spells "LAB3", far away
 * from the allocated machine types, so that it cannot collide with one of
 * them, and it is to be replaced by an allocated number before the platform
 * is submitted.
 */
#define MACH_LAB030		0x4c414233

#endif /* _UAPI_ASM_M68K_BOOTINFO_LAB030_H */
