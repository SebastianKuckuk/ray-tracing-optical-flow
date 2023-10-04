# configuration

TEST_CASE = ray-tracing-optical-flow
BUILD_DIR = ../build
IMAGE_DIR = ../images


# default benchmark parameters

NUM_MGLVL    = 10
NUM_MG       = 10
SUPERSAMPLE  = 2
DT           = 0.001
NUM_TIMESTEP = 3
NUM_MPI      = 4

NX			= 512
NY			= $(NX)

PARAMETERS = ${NUM_MGLVL} ${SUPERSAMPLE} ${NUM_MG} ${DT} ${NUM_TIMESTEP}


# all

targets = \
    $(BUILD_DIR)/$(TEST_CASE)-base \
    $(BUILD_DIR)/$(TEST_CASE)-omp-host \


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(IMAGE_DIR)


# build rules

$(BUILD_DIR)/$(TEST_CASE)-base: main.cpp ray-tracing.h ray-tracing-util.h intersect.h vec3.h optical-flow.h
	mpicxx -O3 -march=native -std=c++17 -o $(BUILD_DIR)/$(TEST_CASE)-base main.cpp

$(BUILD_DIR)/$(TEST_CASE)-omp-host: main.cpp ray-tracing.h ray-tracing-util.h intersect.h vec3.h optical-flow.h
	mpicxx -O3 -march=native -std=c++17 -fopenmp -o $(BUILD_DIR)/$(TEST_CASE)-omp-host main.cpp


# aliases without build directory

.PHONY: $(TEST_CASE)-base
$(TEST_CASE)-base: $(BUILD_DIR)/$(TEST_CASE)-base

.PHONY: $(TEST_CASE)-omp-host
$(TEST_CASE)-omp-host: $(BUILD_DIR)/$(TEST_CASE)-omp-host


# automated benchmark target

.PHONY: bench
bench: all
#	@echo "Base:"
#	mpirun -n $(NUMMPI) $(BUILD_DIR)/$(TEST_CASE)-base $(PARAMETERS)
#	@echo ""
	@echo "OpenMP Host:"
	mpirun -n $(NUM_MPI) $(BUILD_DIR)/$(TEST_CASE)-omp-host $(PARAMETERS)
	@echo ""
	for img in $(IMAGE_DIR)/*.raw ; do \
		python3 convertToImage.py --input $${img} --output $${img}.png --width $(NX) --height $(NY) ; \
	done


# clean target

.PHONY: clean
clean:
	rm $(targets)
