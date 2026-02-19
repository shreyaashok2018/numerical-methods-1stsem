# ============================================================================
# Makefile for EE1103 Static Library
# Author: EE25B126
# ============================================================================

# --- Compiler and Flags ---
CC = gcc
# CFLAGS: -Wall (all warnings), -Wextra (more warnings), -O2 (optimize), 
#         -g (debug symbols), -std=c11 (use C11 standard)
CFLAGS = -Wall -Wextra -O2 -g -std=c11

# --- Archiver and Flags ---
AR = ar
# ARFLAGS: r (replace old files), c (create archive), s (create index)
ARFLAGS = rcs

# --- File Definitions ---

# Target Library Name
TARGET = ee25b126_ee1103.a

# Header File (Dependency)
HEADER = ee25b126_ee1103.h

# Automatically find all .c files starting with 'ee25b126_'
SRCS = $(wildcard ee25b126_*.c)

# Create a list of object files (.o) from the .c file list
OBJS = $(SRCS:.c=.o)

# --- Rules ---

# Phony targets are not files (e.g., 'all', 'clean')
.PHONY: all clean

# Default rule: 'make' or 'make all' will build the target library
all: $(TARGET)

# Rule to build the static library
# Depends on all object files
$(TARGET): $(OBJS)
	@echo "AR  $@"
	$(AR) $(ARFLAGS) $@ $^
	@echo "Static library $(TARGET) created successfully."

# Pattern rule to compile any .c file into a .o file
# Depends on the .c file and the main .h file
# If you change the .h file, all .o files will be re-compiled
%.o: %.c $(HEADER)
	@echo "CC  $<"
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to clean up build files
clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET)
	@echo "Cleanup complete."

# Also create a lib-prefixed archive for standard -l linking
lib: $(TARGET)
	@echo "Creating libee25b126_ee1103.a for -l linking"
	@cp $(TARGET) libee25b126_ee1103.a

.PHONY: lib
