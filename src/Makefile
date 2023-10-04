# configuration

TEST_CASE = ray-tracing-optical-flow
BUILD_DIR = ../build


# default benchmark parameters

NUMLVL      = 10
SUPERSAMPLE = 2
NUMMG       = 10

PARAMETERS = ${NUMLVL} ${SUPERSAMPLE} ${NUMMG}


# all

targets = \
    $(BUILD_DIR)/$(TEST_CASE)-base \
    $(BUILD_DIR)/$(TEST_CASE)-omp-host \


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR)


# build rules

$(BUILD_DIR)/$(TEST_CASE)-base: main.cpp ray-tracing.h ray-tracing-util.h intersect.h vec3.h optical-flow.h
	g++ -O3 -march=native -std=c++17 -o $(BUILD_DIR)/$(TEST_CASE)-base main.cpp  -lm -lpthread -lX11 -ljpeg -lpng

$(BUILD_DIR)/$(TEST_CASE)-omp-host: main.cpp ray-tracing.h ray-tracing-util.h intersect.h vec3.h optical-flow.h
	g++ -O3 -march=native -std=c++17 -fopenmp -o $(BUILD_DIR)/$(TEST_CASE)-omp-host main.cpp  -lm -lpthread -lX11 -ljpeg -lpng


# aliases without build directory

.PHONY: $(TEST_CASE)-base
$(TEST_CASE)-base: $(BUILD_DIR)/$(TEST_CASE)-base

.PHONY: $(TEST_CASE)-omp-host
$(TEST_CASE)-omp-host: $(BUILD_DIR)/$(TEST_CASE)-omp-host


# automated benchmark target

.PHONY: bench
bench: all
	@echo "Base:"
	$(BUILD_DIR)/$(TEST_CASE)-base $(PARAMETERS)
	#pnmtopng result.pnm > result-base.png
	#rm result.pnm
	#python3 convertToImage.py --input result.raw --output result-base-color.png --width $(NX) --height $(NY)
	@echo ""
	@echo "OpenMP Host:"
	$(BUILD_DIR)/$(TEST_CASE)-omp-host $(PARAMETERS)
	#pnmtopng result.pnm > result-base.png
	#rm result.pnm
	#python3 convertToImage.py --input result.raw --output result-base-color.png --width $(NX) --height $(NY)
	@echo ""


# clean target

.PHONY: clean
clean:
	rm $(targets)
