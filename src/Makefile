# configuration

TEST_CASE = ray-tracing-optical-flow
BUILD_DIR = ../build
DATA_DIR  = ../imageData
IMAGE_DIR = ../images
VIDEO_DIR = ../videos


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
    $(BUILD_DIR)/$(TEST_CASE)-base-prof \
    $(BUILD_DIR)/$(TEST_CASE)-omp-host-prof


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR) $(DATA_DIR) $(IMAGE_DIR) $(VIDEO_DIR)


# build rules

$(BUILD_DIR)/$(TEST_CASE)-base: main.cpp ray-tracing.h ray-tracing-util.h intersect.h image.h vec3.h optical-flow.h timer.h
	mpicxx -O3 -march=native -std=c++17 -o $(BUILD_DIR)/$(TEST_CASE)-base main.cpp

$(BUILD_DIR)/$(TEST_CASE)-omp-host: main.cpp ray-tracing.h ray-tracing-util.h intersect.h image.h vec3.h optical-flow.h timer.h
	mpicxx -O3 -march=native -std=c++17 -fopenmp -o $(BUILD_DIR)/$(TEST_CASE)-omp-host main.cpp

$(BUILD_DIR)/$(TEST_CASE)-base-prof: main.cpp ray-tracing.h ray-tracing-util.h intersect.h image.h vec3.h optical-flow.h timer.h
	mpicxx -O3 -march=native -std=c++17 -DLIKWID_PERFMON -I$(LIKWID_INCDIR) -L$(LIKWID_LIBDIR) -o $(BUILD_DIR)/$(TEST_CASE)-base-prof main.cpp -llikwid

$(BUILD_DIR)/$(TEST_CASE)-omp-host-prof: main.cpp ray-tracing.h ray-tracing-util.h intersect.h image.h vec3.h optical-flow.h timer.h
	mpicxx -O3 -march=native -std=c++17 -DLIKWID_PERFMON -I$(LIKWID_INCDIR) -L$(LIKWID_LIBDIR) -fopenmp -o $(BUILD_DIR)/$(TEST_CASE)-omp-host-prof main.cpp -llikwid


# aliases without build directory

.PHONY: $(TEST_CASE)-base
$(TEST_CASE)-base: $(BUILD_DIR)/$(TEST_CASE)-base

.PHONY: $(TEST_CASE)-omp-host
$(TEST_CASE)-omp-host: $(BUILD_DIR)/$(TEST_CASE)-omp-host

.PHONY: $(TEST_CASE)-base-prof
$(TEST_CASE)-base-prof: $(BUILD_DIR)/$(TEST_CASE)-base-prof

.PHONY: $(TEST_CASE)-omp-host-prof
$(TEST_CASE)-omp-host-prof: $(BUILD_DIR)/$(TEST_CASE)-omp-host-prof


# automated benchmark target

.PHONY: bench
bench: all
	@echo "Base:"
	mpirun -n $(NUM_MPI) $(BUILD_DIR)/$(TEST_CASE)-base $(PARAMETERS)
	@echo ""
#	@echo "OpenMP Host:"
#	mpirun -n $(NUM_MPI) $(BUILD_DIR)/$(TEST_CASE)-omp-host $(PARAMETERS)
#	@echo ""
#	@echo "Serial w/ Profiling:"
#	likwid-perfctr -g MEM_DP -C S0:0 -m $(BUILD_DIR)/$(TEST_CASE)-base-prof $(PARAMETERS)
#	@echo ""
#	@echo "Base w/ Profiling:"
#	likwid-mpirun -n $(NUM_MPI) -g MEM_DP -m $(BUILD_DIR)/$(TEST_CASE)-base-prof $(PARAMETERS)
#	@echo ""

	if [ 0 != ${PRINT_IMG} ] ; then \
		for img in $(DATA_DIR)/ray-tracing-*.raw ; do \
			python3 convertToImage.py --input $${img} --output $(IMAGE_DIR)/$$(basename $$img .raw).png --width $(NX) --height $(NY) --scale off ; \
		done ; \
		for img in $(DATA_DIR)/optical-flow-*.raw ; do \
			python3 convertToImage.py --input $${img} --output $(IMAGE_DIR)/$$(basename $$img .raw).png --width $(NX) --height $(NY) --scale absAuto ; \
		done ; \
		ffmpeg -y -framerate 30 -pattern_type glob -i '$(IMAGE_DIR)/ray-tracing-*.png'  -c:v libx264 -pix_fmt yuv420p $(VIDEO_DIR)/ray-tracing.mp4 ; \
		ffmpeg -y -framerate 30 -pattern_type glob -i '$(IMAGE_DIR)/optical-flow-*.png' -c:v libx264 -pix_fmt yuv420p $(VIDEO_DIR)/optical-flow.mp4 ; \
	fi
# clean target

.PHONY: clean
clean:
	rm $(targets)
