// SPDX-License-Identifier: GPL-2.0
/*
 * m68030-lab board support
 */

#include <linux/clockchips.h>
#include <linux/clocksource.h>
#include <linux/err.h>
#include <linux/init.h>
#include <linux/interrupt.h>
#include <linux/io.h>
#include <linux/platform_device.h>
#include <linux/serial_8250.h>
#include <linux/serial_core.h>
#include <linux/timex.h>

#include <asm/config.h>
#include <asm/irq.h>
#include <asm/lab030.h>
#include <asm/machdep.h>
#include <asm/setup.h>

static void __iomem *lab030_cpld;

static u32 lab030_cpld_read(unsigned int reg)
{
	return in_be32(lab030_cpld + reg);
}

static void lab030_cpld_write(u32 value, unsigned int reg)
{
	out_be32(lab030_cpld + reg, value);
}

static void lab030_get_model(char *model)
{
	sprintf(model, "m68030-lab");
}

static void lab030_reset(void)
{
	local_irq_disable();
	lab030_cpld_write(LAB030_RESET_BOARD, LAB030_CPLD_RESET);
	while (1)
		;
}

static void __init lab030_init_IRQ(void)
{
}

static int lab030_tick_shutdown(struct clock_event_device *evt)
{
	lab030_cpld_write(0, LAB030_CPLD_TICK_CTRL);

	return 0;
}

static int lab030_tick_set_periodic(struct clock_event_device *evt)
{
	lab030_cpld_write(LAB030_TICK_CTRL_ENABLE, LAB030_CPLD_TICK_CTRL);

	return 0;
}

static struct clock_event_device lab030_tick = {
	.name			= "lab030-tick",
	.features		= CLOCK_EVT_FEAT_PERIODIC,
	.rating			= 200,
	.set_state_shutdown	= lab030_tick_shutdown,
	.set_state_periodic	= lab030_tick_set_periodic,
};

static u32 clk_total;

static irqreturn_t lab030_tick_interrupt(int irq, void *dev_id)
{
	lab030_cpld_write(LAB030_TICK_STATUS_PENDING, LAB030_CPLD_TICK_STATUS);
	clk_total += LAB030_TICK_CYCLES;
	lab030_tick.event_handler(&lab030_tick);

	return IRQ_HANDLED;
}

static bool lab030_tick_pending(void)
{
	return lab030_cpld_read(LAB030_CPLD_TICK_STATUS) &
	       LAB030_TICK_STATUS_PENDING;
}

static u64 lab030_read_clk(struct clocksource *cs)
{
	unsigned long flags;
	bool pending, tmp;
	u32 count, ticks;

	local_irq_save(flags);
	/*
	 * PENDING is set when the counter wraps, and the interrupt that adds
	 * the wrap to clk_total clears it.  Read it on both sides of the
	 * counter, in case the counter wraps in between.
	 */
	tmp = lab030_tick_pending();
	count = lab030_cpld_read(LAB030_CPLD_TICK_COUNT);
	pending = lab030_tick_pending();
	if (pending != tmp)
		count = lab030_cpld_read(LAB030_CPLD_TICK_COUNT);
	if (pending)
		count += LAB030_TICK_CYCLES;
	ticks = clk_total + count;
	local_irq_restore(flags);

	return ticks;
}

static struct clocksource lab030_clk = {
	.name	= "lab030-count",
	.rating	= 250,
	.read	= lab030_read_clk,
	.mask	= CLOCKSOURCE_MASK(32),
	.flags	= CLOCK_SOURCE_IS_CONTINUOUS,
};

static void __init lab030_sched_init(void)
{
	if (request_irq(IRQ_AUTO_6, lab030_tick_interrupt, IRQF_TIMER, "timer",
			NULL))
		pr_err("Couldn't register timer interrupt\n");

	clocksource_register_hz(&lab030_clk, LAB030_TICK_CLOCK_FREQ);

	lab030_tick.cpumask = cpumask_of(0);
	clockevents_register_device(&lab030_tick);
}

static unsigned long lab030_random_get_entropy(void)
{
	return lab030_cpld_read(LAB030_CPLD_TICK_COUNT);
}

#ifdef CONFIG_SERIAL_8250_CONSOLE
static int __init lab030_earlycon_setup(struct earlycon_device *device,
					const char *options)
{
	device->port.regshift = LAB030_UART_REGSHIFT;

	return early_serial8250_setup(device, options);
}

EARLYCON_DECLARE(lab030uart, lab030_earlycon_setup);
#endif

void __init config_lab030(void)
{
	char earlycon[32];

	lab030_cpld = ioremap(LAB030_CPLD_BASE, LAB030_CPLD_SIZE);

	snprintf(earlycon, sizeof(earlycon), "lab030uart,mmio,0x%08x",
		 LAB030_UART_BASE);
	setup_earlycon(earlycon);

	mach_sched_init = lab030_sched_init;
	mach_init_IRQ = lab030_init_IRQ;
	mach_get_model = lab030_get_model;
	mach_reset = lab030_reset;
	mach_random_get_entropy = lab030_random_get_entropy;
}

static const struct plat_serial8250_port lab030_uart_data[] __initconst = {
	{
		.mapbase	= LAB030_UART_BASE,
		.irq		= IRQ_AUTO_4,
		.uartclk	= LAB030_UART_CLOCK_FREQ,
		.regshift	= LAB030_UART_REGSHIFT,
		.iotype		= UPIO_MEM,
		.flags		= UPF_BOOT_AUTOCONF | UPF_IOREMAP,
	},
	{ }
};

static int __init lab030_platform_init(void)
{
	struct platform_device *pdev;

	if (!MACH_IS_LAB030)
		return -ENODEV;

	pdev = platform_device_register_data(NULL, "serial8250",
					     PLAT8250_DEV_PLATFORM,
					     lab030_uart_data,
					     sizeof(lab030_uart_data));

	return PTR_ERR_OR_ZERO(pdev);
}
arch_initcall(lab030_platform_init);
