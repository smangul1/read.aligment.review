################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../genomeindexer/QSufSort.c \
../genomeindexer/bntseq.c \
../genomeindexer/bwt.c \
../genomeindexer/bwt_gen.c \
../genomeindexer/bwtindex.c \
../genomeindexer/bwtio.c \
../genomeindexer/bwtmisc.c \
../genomeindexer/is.c \
../genomeindexer/utils.c 

OBJS += \
./genomeindexer/QSufSort.o \
./genomeindexer/bntseq.o \
./genomeindexer/bwt.o \
./genomeindexer/bwt_gen.o \
./genomeindexer/bwtindex.o \
./genomeindexer/bwtio.o \
./genomeindexer/bwtmisc.o \
./genomeindexer/is.o \
./genomeindexer/utils.o 

C_DEPS += \
./genomeindexer/QSufSort.d \
./genomeindexer/bntseq.d \
./genomeindexer/bwt.d \
./genomeindexer/bwt_gen.d \
./genomeindexer/bwtindex.d \
./genomeindexer/bwtio.d \
./genomeindexer/bwtmisc.d \
./genomeindexer/is.d \
./genomeindexer/utils.d 


# Each subdirectory must supply rules for building sources it contributes
genomeindexer/%.o: ../genomeindexer/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


