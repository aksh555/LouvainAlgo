RM := rm

USER_OBJS :=
LIBS :=

-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

all: find-communities-OpenMP

%.o: %.c
	${Q} gcc -fopenmp -O2 -g -w -c -fopenmp -o $@ $<

find-communities-OpenMP: $(OBJS) $(USER_OBJS)
	${Q} gcc -fopenmp -w -o "find-communities-OpenMP" $(OBJS) $(USER_OBJS)

clean:
	-$(RM) $(EXECUTABLES)$(OBJS) find-communities-OpenMP

.PHONY: all clean dependents
.SECONDARY:
