################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../dnaindexer/QSufSort.c \
../dnaindexer/bntseq.c \
../dnaindexer/bwt.c \
../dnaindexer/bwt_gen.c \
../dnaindexer/bwtindex.c \
../dnaindexer/bwtio.c \
../dnaindexer/bwtmisc.c \
../dnaindexer/is.c \
../dnaindexer/utils.c 

OBJS += \
./dnaindexer/QSufSort.o \
./dnaindexer/bntseq.o \
./dnaindexer/bwt.o \
./dnaindexer/bwt_gen.o \
./dnaindexer/bwtindex.o \
./dnaindexer/bwtio.o \
./dnaindexer/bwtmisc.o \
./dnaindexer/is.o \
./dnaindexer/utils.o 

C_DEPS += \
./dnaindexer/QSufSort.d \
./dnaindexer/bntseq.d \
./dnaindexer/bwt.d \
./dnaindexer/bwt_gen.d \
./dnaindexer/bwtindex.d \
./dnaindexer/bwtio.d \
./dnaindexer/bwtmisc.d \
./dnaindexer/is.d \
./dnaindexer/utils.d 


# Each subdirectory must supply rules for building sources it contributes
dnaindexer/%.o: ../dnaindexer/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


