 OBJ_SRCS := 
ASM_SRCS := 
C_SRCS := 
O_SRCS := 
S_UPPER_SRCS := 
EXECUTABLES := 
OBJS := 
C_DEPS := 
SUBDIRS := \
. \
 
 C_SRCS += \
../community-development.c \
../community-computation-weighted.c \
../community-computation-weighted-sequential.c \
../community-exchange.c \
../dynamic-weighted-graph.c \
../execution-handler.c \
../input-handler.c \
../main.c \
../sorted-linked-list.c \
../utilities.c

OBJS += \
./community-development.o \
./community-computation-weighted.o \
./community-exchange.o \
./dynamic-weighted-graph.o \
./execution-handler.o \
./input-handler.o \
./main.o \
./sorted-linked-list.o \
./utilities.o

C_DEPS += \
./community-development.d \
./community-computation-weighted-sequential.d \
./community-computation-weighted.d \
./community-exchange.d \
./dynamic-weighted-graph.d \
./execution-handler.d \
./input-handler.d \
./main.d \
./sorted-linked-list.d \
./utilities.d

# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -fopenmp -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


