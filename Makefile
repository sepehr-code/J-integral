# J-Integral Calculation Makefile
# PhD Research Implementation (2010)

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g
LDFLAGS = -lm

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj
BINDIR = bin
DATADIR = data
EXAMPLEDIR = examples
TESTDIR = tests

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
MAIN_OBJ = $(OBJDIR)/main.o
LIB_OBJECTS = $(filter-out $(MAIN_OBJ), $(OBJECTS))

# Target executable
TARGET = $(BINDIR)/j_integral

# Library name
LIBNAME = libj_integral.a
LIBRARY = $(BINDIR)/$(LIBNAME)

# Example targets
EXAMPLE_SOURCES = $(wildcard $(EXAMPLEDIR)/*.c)
EXAMPLE_TARGETS = $(EXAMPLE_SOURCES:$(EXAMPLEDIR)/%.c=$(BINDIR)/%)

# Test targets
TEST_SOURCES = $(wildcard $(TESTDIR)/*.c)
TEST_TARGETS = $(TEST_SOURCES:$(TESTDIR)/%.c=$(BINDIR)/%)

# Default target
.PHONY: all
all: directories $(TARGET) $(LIBRARY)

# Create necessary directories
.PHONY: directories
directories:
	@mkdir -p $(OBJDIR) $(BINDIR) $(DATADIR) $(EXAMPLEDIR) $(TESTDIR)

# Main executable
$(TARGET): $(OBJECTS)
	@echo "Linking $@..."
	@$(CC) $(OBJECTS) -o $@ $(LDFLAGS)
	@echo "Built $@ successfully"

# Static library (without main.o)
$(LIBRARY): $(LIB_OBJECTS)
	@echo "Creating library $@..."
	@ar rcs $@ $(LIB_OBJECTS)
	@echo "Library $@ created successfully"

# Object files
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

# Examples
.PHONY: examples
examples: directories $(LIBRARY) $(EXAMPLE_TARGETS)

$(BINDIR)/%: $(EXAMPLEDIR)/%.c $(LIBRARY)
	@echo "Building example $@..."
	@$(CC) $(CFLAGS) -I$(INCDIR) $< -L$(BINDIR) -lj_integral -o $@ $(LDFLAGS)

# Tests
.PHONY: tests
tests: directories $(LIBRARY) $(TEST_TARGETS)

$(BINDIR)/%: $(TESTDIR)/%.c $(LIBRARY)
	@echo "Building test $@..."
	@$(CC) $(CFLAGS) -I$(INCDIR) $< -L$(BINDIR) -lj_integral -o $@ $(LDFLAGS)

# Run tests
.PHONY: check
check: tests
	@echo "Running tests..."
	@for test in $(TEST_TARGETS); do \
		echo "Running $$test..."; \
		$$test || exit 1; \
	done
	@echo "All tests passed!"

# Sample data files
.PHONY: sample-data
sample-data: directories
	@echo "Generating sample data files..."
	@echo "# Sample mesh file" > $(DATADIR)/sample_mesh.dat
	@echo "21 21" >> $(DATADIR)/sample_mesh.dat
	@echo "2.0 2.0" >> $(DATADIR)/sample_mesh.dat
	@echo "0.5 1.0 1.0 1.0" >> $(DATADIR)/sample_mesh.dat
	@echo "# Sample field data file" > $(DATADIR)/sample_field.dat
	@echo "# node_id x y ux uy sigma_xx sigma_yy sigma_xy" >> $(DATADIR)/sample_field.dat
	@echo "Sample data files created in $(DATADIR)/"

# Documentation
.PHONY: docs
docs:
	@echo "Generating documentation..."
	@if command -v doxygen >/dev/null 2>&1; then \
		doxygen Doxyfile 2>/dev/null || echo "Doxygen configuration not found"; \
	else \
		echo "Doxygen not installed. Install with: sudo apt-get install doxygen"; \
	fi

# Installation
.PHONY: install
install: all
	@echo "Installing J-integral calculator..."
	@sudo cp $(TARGET) /usr/local/bin/
	@sudo cp $(LIBRARY) /usr/local/lib/
	@sudo mkdir -p /usr/local/include/j_integral
	@sudo cp $(INCDIR)/*.h /usr/local/include/j_integral/
	@echo "Installation completed. Run 'j_integral --help' for usage."

# Uninstall
.PHONY: uninstall
uninstall:
	@echo "Uninstalling J-integral calculator..."
	@sudo rm -f /usr/local/bin/j_integral
	@sudo rm -f /usr/local/lib/$(LIBNAME)
	@sudo rm -rf /usr/local/include/j_integral
	@echo "Uninstallation completed."

# Clean targets
.PHONY: clean
clean:
	@echo "Cleaning build files..."
	@rm -rf $(OBJDIR) $(BINDIR)
	@echo "Clean completed."

.PHONY: clean-data
clean-data:
	@echo "Cleaning data files..."
	@rm -f $(DATADIR)/*.dat $(DATADIR)/*.txt $(DATADIR)/*.csv
	@echo "Data files cleaned."

.PHONY: clean-all
clean-all: clean clean-data
	@echo "Full clean completed."

# Development targets
.PHONY: debug
debug: CFLAGS += -DDEBUG -g3 -O0
debug: clean all

.PHONY: release
release: CFLAGS += -DNDEBUG -O3
release: clean all

.PHONY: profile
profile: CFLAGS += -pg
profile: LDFLAGS += -pg
profile: clean all

# Static analysis
.PHONY: lint
lint:
	@echo "Running static analysis..."
	@if command -v cppcheck >/dev/null 2>&1; then \
		cppcheck --enable=all --std=c99 $(SRCDIR)/ $(INCDIR)/; \
	else \
		echo "cppcheck not installed. Install with: sudo apt-get install cppcheck"; \
	fi

# Memory check
.PHONY: memcheck
memcheck: debug
	@echo "Running memory check..."
	@if command -v valgrind >/dev/null 2>&1; then \
		valgrind --leak-check=full --show-leak-kinds=all $(TARGET); \
	else \
		echo "valgrind not installed. Install with: sudo apt-get install valgrind"; \
	fi

# Performance profiling
.PHONY: perf
perf: profile
	@echo "Running performance analysis..."
	@if command -v gprof >/dev/null 2>&1; then \
		$(TARGET) && gprof $(TARGET) gmon.out > performance_report.txt; \
		echo "Performance report saved to performance_report.txt"; \
	else \
		echo "gprof not available. Install with: sudo apt-get install binutils"; \
	fi

# Code formatting
.PHONY: format
format:
	@echo "Formatting code..."
	@if command -v clang-format >/dev/null 2>&1; then \
		find $(SRCDIR) $(INCDIR) -name "*.c" -o -name "*.h" | xargs clang-format -i; \
		echo "Code formatting completed."; \
	else \
		echo "clang-format not installed. Install with: sudo apt-get install clang-format"; \
	fi

# Help target
.PHONY: help
help:
	@echo "J-Integral Calculation Makefile"
	@echo "==============================="
	@echo ""
	@echo "Available targets:"
	@echo "  all          - Build main executable and library (default)"
	@echo "  clean        - Remove build files"
	@echo "  clean-all    - Remove all generated files"
	@echo "  examples     - Build example programs"
	@echo "  tests        - Build test programs"
	@echo "  check        - Run all tests"
	@echo "  sample-data  - Generate sample input files"
	@echo "  install      - Install to system directories"
	@echo "  uninstall    - Remove from system directories"
	@echo ""
	@echo "Development targets:"
	@echo "  debug        - Build with debug symbols"
	@echo "  release      - Build optimized release version"
	@echo "  profile      - Build with profiling support"
	@echo "  lint         - Run static analysis"
	@echo "  memcheck     - Run memory leak detection"
	@echo "  perf         - Run performance profiling"
	@echo "  format       - Format source code"
	@echo "  docs         - Generate documentation"
	@echo ""
	@echo "Usage examples:"
	@echo "  make                    # Build everything"
	@echo "  make debug              # Debug build"
	@echo "  make tests check        # Build and run tests"
	@echo "  make install            # Install system-wide"
	@echo "  make clean-all          # Clean everything"

# Dependencies
-include $(OBJECTS:.o=.d)

# Generate dependency files
$(OBJDIR)/%.d: $(SRCDIR)/%.c
	@mkdir -p $(OBJDIR)
	@$(CC) $(CFLAGS) -I$(INCDIR) -MM -MT $(@:.d=.o) $< > $@

.PHONY: deps
deps: $(OBJECTS:.o=.d)

# Print variables (for debugging Makefile)
.PHONY: print-vars
print-vars:
	@echo "CC = $(CC)"
	@echo "CFLAGS = $(CFLAGS)"
	@echo "LDFLAGS = $(LDFLAGS)"
	@echo "SOURCES = $(SOURCES)"
	@echo "OBJECTS = $(OBJECTS)"
	@echo "TARGET = $(TARGET)"
	@echo "LIBRARY = $(LIBRARY)" 