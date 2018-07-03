################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../core/Aligner.cpp \
../core/BSOptions.cpp \
../core/BWT.cpp \
../core/CSOptions.cpp \
../core/CigarAlign.cpp \
../core/Genome.cpp \
../core/MemEngine.cpp \
../core/Options.cpp \
../core/PairedEnd.cpp \
../core/SAM.cpp \
../core/Seed.cpp \
../core/SeqFileParser.cpp \
../core/Sequence.cpp \
../core/SingleEnd.cpp \
../core/SuffixArray.cpp \
../core/Utils.cpp \
../core/main.cpp 

OBJS += \
./core/Aligner.o \
./core/BSOptions.o \
./core/BWT.o \
./core/CSOptions.o \
./core/CigarAlign.o \
./core/Genome.o \
./core/MemEngine.o \
./core/Options.o \
./core/PairedEnd.o \
./core/SAM.o \
./core/Seed.o \
./core/SeqFileParser.o \
./core/Sequence.o \
./core/SingleEnd.o \
./core/SuffixArray.o \
./core/Utils.o \
./core/main.o 

CPP_DEPS += \
./core/Aligner.d \
./core/BSOptions.d \
./core/BWT.d \
./core/CSOptions.d \
./core/CigarAlign.d \
./core/Genome.d \
./core/MemEngine.d \
./core/Options.d \
./core/PairedEnd.d \
./core/SAM.d \
./core/Seed.d \
./core/SeqFileParser.d \
./core/Sequence.d \
./core/SingleEnd.d \
./core/SuffixArray.d \
./core/Utils.d \
./core/main.d 


# Each subdirectory must supply rules for building sources it contributes
core/%.o: ../core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" -I /home/yongchao/cuda-workspace/cushaw3-v3.0.2/bamreader -msse4 -fopenmp -I /home/yongchao/cuda-workspace/cushaw3-v3.0.2/genomeindexer
	@echo 'Finished building: $<'
	@echo ' '


