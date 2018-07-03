c_sources = bntseq.c bwt.c bwt_gen.c bwtindex.c QSufSort.c utils.c bwtmisc.c bwtio.c is.c

objs_dir = objs
src_dir = .
objs += $(patsubst %.c, $(objs_dir)/%.c.o, $(c_sources))

AR = ar
CC = gcc
CCFLAGS = -O3 -funroll-loops -I $(src_dir) 

target=libgenomeindexer.a

all: dir $(objs)
	$(AR) rcs $(target) $(objs)

dir:
	mkdir -p $(objs_dir)

clean:
	rm -rf $(objs_dir)
	rm -f $(target)

$(objs_dir)/%.c.o: $(src_dir)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $<


