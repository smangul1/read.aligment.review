sources=Aligner.cpp BWT.cpp CigarAlign.cpp Genome.cpp main.cpp \
	MemEngine.cpp Options.cpp PairedEnd.cpp SAM.cpp Seed.cpp \
	SeqFileParser.cpp Sequence.cpp SingleEnd.cpp SuffixArray.cpp \
	Utils.cpp CSOptions.cpp BSOptions.cpp

root_dir = .
core_dir = core
objs_dir = objs
objs += $(patsubst %.cpp, $(objs_dir)/%.cpp.o, $(sources))
mylibs = -lpthread -lm -lz
EXEC = cushaw3

#compressed input
MACROS = -DCOMPRESSED_INPUT

#have the support for Intel SSSE3?
have_ssse3 = 1
ifneq ($(have_ssse3), 0)
	MACROS += -DHAVE_SSSE3
	SSE_ARCH = -msse4
else
	MACROS +=
	SSE_ARCH = -msse2
endif

#compile optionsK
CXXFLAGS += -static -O3 -Wall -fopenmp -funroll-loops $(MACROS) $(SSE_ARCH)
CXX = g++


#bam reader
BAMREADER = $(root_dir)/bamreader
DNAINDEXER = $(root_dir)/genomeindexer
mylibs += -L $(BAMREADER) -lbamreader -L$(DNAINDEXER) -lgenomeindexer
CXXFLAGS += -I $(BAMREADER) -I $(DNAINDEXER)

all: dir bamr dnai $(objs)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(objs) $(mylibs)
	strip $(EXEC)

dir:
	mkdir -p $(objs_dir)

bamr:
	make -C $(BAMREADER)

dnai:
	make -C $(DNAINDEXER)

clean:
	-rm -rf $(objs_dir) $(EXEC)
	make -C $(BAMREADER) clean
	make -C $(DNAINDEXER) clean

$(objs_dir)/%.cpp.o: $(core_dir)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

