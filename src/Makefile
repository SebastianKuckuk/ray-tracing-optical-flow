# configuration

TEST_CASE = ray-tracing-optical-flow
BUILD_DIR = ../build
IMAGE_DIR = ../images


# default benchmark parameters

NUM_MGLVL    = 10
NUM_MG       = 10
SUPERSAMPLE  = 2
DT           = 0.01
NUM_TIMESTEP = 4
NUM_MPI      = 8
NUM_REP      = 1
PRINT_IMG    = 0

NX			= 512   # has to be 2 ** NUM_MGLVL
NY			= $(NX)

PARAMETERS = ${NUM_MGLVL} ${SUPERSAMPLE} ${NUM_MG} ${DT} ${NUM_TIMESTEP} ${NUM_REP} ${PRINT_IMG}


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

$(BUILD_DIR)/$(TEST_CASE)-base: main.cpp ray-tracing.h ray-tracing-util.h intersect.h image.h vec3.h optical-flow.h timer.h
	mpicxx -O3 -march=native -std=c++17 -o $(BUILD_DIR)/$(TEST_CASE)-base main.cpp

$(BUILD_DIR)/$(TEST_CASE)-omp-host: main.cpp ray-tracing.h ray-tracing-util.h intersect.h image.h vec3.h optical-flow.h timer.h
	mpicxx -O3 -march=native -std=c++17 -fopenmp -o $(BUILD_DIR)/$(TEST_CASE)-omp-host main.cpp


# aliases without build directory

.PHONY: $(TEST_CASE)-base
$(TEST_CASE)-base: $(BUILD_DIR)/$(TEST_CASE)-base

.PHONY: $(TEST_CASE)-omp-host
$(TEST_CASE)-omp-host: $(BUILD_DIR)/$(TEST_CASE)-omp-host


# automated benchmark target

.PHONY: bench
bench: all
	@echo "Base:"
	mpirun -n $(NUM_MPI) $(BUILD_DIR)/$(TEST_CASE)-base $(PARAMETERS)
	@echo ""
#	@echo "OpenMP Host:"
#	mpirun -n $(NUM_MPI) $(BUILD_DIR)/$(TEST_CASE)-omp-host $(PARAMETERS)
#	@echo ""
	if [ 0 != ${PRINT_IMG} ] ; then \
		for img in $(IMAGE_DIR)/*.raw ; do \
			python3 convertToImage.py --input $${img} --output $${img}.png --width $(NX) --height $(NY) ; \
		done ; \
		ffmpeg -framerate 30 -pattern_type glob -i '../images/optical-flow-*.png' -c:v libx264 -pix_fmt yuv420p ../images/optical-flow.mp4 ; \
		ffmpeg -framerate 30 -pattern_type glob -i '../images/ray-tracing-*.png' -c:v libx264 -pix_fmt yuv420p ../images/ray-tracing.mp4 ; \
	fi
# clean target

.PHONY: clean
clean:
	rm $(targets)
