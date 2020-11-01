RM := rm

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: find-communities-OpenMP

%.o: %.c
	gcc -fopenmp -O2 -g -w -c -fopenmp -o $@ $<

# Tool invocations
find-communities-OpenMP: $(OBJS) $(USER_OBJS)
	gcc -fopenmp -w -o "find-communities-OpenMP" $(OBJS) $(USER_OBJS)

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS) find-communities-OpenMP

.PHONY: all clean dependents
.SECONDARY:
