sources = bam_aux.c  bam.c  bam_import.c  bgzf.c  faidx.c  kstring.c razf.c  sam.c  sam_header.c

src_dir = .
objs_dir = objs
objs += $(patsubst %.c, $(objs_dir)/%.c.o, $(sources))
CCFLAGS = -O3 -Wall -funroll-loops
CC = gcc
AR = ar
target = libbamreader.a

all: dir $(objs)
	$(AR) rcs $(target) $(objs)

dir:
	mkdir -p $(objs_dir)

clean:
	-rm -rf $(objs_dir) $(target)

$(objs_dir)/%.c.o: $(src_dir)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $<

