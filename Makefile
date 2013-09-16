all:
	gcc main.c calc_bonds.c calc_hueristics.c get_input.c init_fp.c SmoothPeak.c -lm -o Hueristics
