# Ordinary differential equation solver, Runge-Kutta-England technique.
# Copyright © 1988 Free Software Foundation, Inc.
# François Pinard <pinard@iro.umontreal.ca>, 1988.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 1, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

CC = gcc
CFLAGS	= -g
LDFLAGS	=

example pregithub: example.o rke.o
	$(CC) $(LDFLAGS) example.o rke.o -o example -lm

example.o: rke.h example.c
rke.o: rke.h rke.c

clean:
	rm -f *~ *.o example
