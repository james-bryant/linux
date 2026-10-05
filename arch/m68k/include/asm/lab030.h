/* SPDX-License-Identifier: GPL-2.0 */
/*
 * m68030-lab board definitions
 *
 * The address map and the CPLD registers are provisional.
 */

#ifndef _ASM_M68K_LAB030_H
#define _ASM_M68K_LAB030_H

/* CPLD registers: big-endian longwords */
#define LAB030_CPLD_BASE		0xff010000
#define LAB030_CPLD_SIZE		0x1000

#define LAB030_CPLD_OVERLAY		0x00
#define LAB030_CPLD_BOARD_ID		0x04
#define LAB030_CPLD_TICK_CTRL		0x08
#define LAB030_CPLD_TICK_STATUS		0x0c
#define LAB030_CPLD_TICK_COUNT		0x10
#define LAB030_CPLD_RESET		0x14

#define LAB030_TICK_CTRL_ENABLE		0x1
#define LAB030_TICK_STATUS_PENDING	0x1	/* write 1 to clear */
#define LAB030_RESET_BOARD		0x1

/*
 * TICK_COUNT runs at 1.8432 MHz and wraps every 18432 counts.  Each wrap
 * is a tick: 100 a second.
 */
#define LAB030_TICK_CLOCK_FREQ		1843200
#define LAB030_TICK_CYCLES		18432

/* 16550A console UART: one register every four bytes */
#define LAB030_UART_BASE		0xff030000
#define LAB030_UART_REGSHIFT		2
#define LAB030_UART_CLOCK_FREQ		1843200

/*
 * CompactFlash socket in true-IDE mode: one register every two bytes.  The
 * control block has a single register, device control and alternate status.
 */
#define LAB030_CF_CMD_BASE		0xff040000
#define LAB030_CF_CMD_SIZE		0x10
#define LAB030_CF_CTL_BASE		0xff050000
#define LAB030_CF_CTL_DEVCTL		0x0c
#define LAB030_CF_CTL_DEVCTL_SIZE	0x2
#define LAB030_CF_REGSHIFT		1

#endif /* _ASM_M68K_LAB030_H */
