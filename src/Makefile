# configuration

TEST_CASE = ray-tracing
BUILD_DIR = ../../build/$(TEST_CASE)


# default benchmark parameters

NX    = 514
NY    = ${NX}
NWARM = 0
NIT   = 1

PARAMETERS = ${NX} ${NY} ${NWARM} ${NIT}


# all

targets = \
    $(BUILD_DIR)/$(TEST_CASE)-base \
#    $(BUILD_DIR)/$(TEST_CASE)-omp-host \
#    $(BUILD_DIR)/$(TEST_CASE)-cuda \
#    $(BUILD_DIR)/$(TEST_CASE)-openacc \
#    $(BUILD_DIR)/$(TEST_CASE)-omp-target \
#    $(BUILD_DIR)/$(TEST_CASE)-cpp-17 \
#    $(BUILD_DIR)/$(TEST_CASE)-kokkos-serial \
#    $(BUILD_DIR)/$(TEST_CASE)-kokkos-omp-host \
#    $(BUILD_DIR)/$(TEST_CASE)-kokkos-cuda \
#    $(BUILD_DIR)/$(TEST_CASE)-sycl-implicit \
#    $(BUILD_DIR)/$(TEST_CASE)-sycl-explicit \
#    $(BUILD_DIR)/$(TEST_CASE)-hip


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR)


# build rules

$(BUILD_DIR)/$(TEST_CASE)-base: $(TEST_CASE)-base.cpp $(TEST_CASE)-util.h vec3.h intersect.h conjugate-gradient-base.h multigrid-base.h
	g++ -O3 -march=native -std=c++17 -o $(BUILD_DIR)/$(TEST_CASE)-base $(TEST_CASE)-base.cpp  -lm -lpthread -lX11 -ljpeg -lpng

$(BUILD_DIR)/$(TEST_CASE)-omp-host: $(TEST_CASE)-omp-host.cpp $(TEST_CASE)-util.h ../util.h
	g++ -O3 -march=native -std=c++17 -fopenmp -o $(BUILD_DIR)/$(TEST_CASE)-omp-host $(TEST_CASE)-omp-host.cpp

$(BUILD_DIR)/$(TEST_CASE)-cuda: $(TEST_CASE)-cuda.cu $(TEST_CASE)-util.h ../util.h ../cuda-util.h
	nvc++ -O3 -std=c++17 -o $(BUILD_DIR)/$(TEST_CASE)-cuda $(TEST_CASE)-cuda.cu

$(BUILD_DIR)/$(TEST_CASE)-openacc: $(TEST_CASE)-openacc.cpp $(TEST_CASE)-util.h ../util.h
	nvc++ -fast -acc -target=gpu -o $(BUILD_DIR)/$(TEST_CASE)-openacc $(TEST_CASE)-openacc.cpp

$(BUILD_DIR)/$(TEST_CASE)-omp-target: $(TEST_CASE)-omp-target.cpp $(TEST_CASE)-util.h ../util.h
	nvc++ -O3 -mp=gpu -o $(BUILD_DIR)/$(TEST_CASE)-omp-target $(TEST_CASE)-omp-target.cpp # -Minfo=mp

$(BUILD_DIR)/$(TEST_CASE)-cpp-17: $(TEST_CASE)-cpp-17.cpp $(TEST_CASE)-util.h ../util.h
	nvc++ -O3 -march=native -std=c++17 -stdpar=gpu -gpu=cc75 -o $(BUILD_DIR)/$(TEST_CASE)-cpp-17 $(TEST_CASE)-cpp-17.cpp

$(BUILD_DIR)/$(TEST_CASE)-kokkos-serial: $(TEST_CASE)-kokkos.cpp $(TEST_CASE)-util.h ../util.h
	g++ -O3 -march=native -std=c++17 -I${HOME}/kokkos/install-serial/include -c $(TEST_CASE)-kokkos.cpp -o $(BUILD_DIR)/$(TEST_CASE)-kokkos-serial.o
	g++ $(BUILD_DIR)/$(TEST_CASE)-kokkos-serial.o -L${HOME}/kokkos/install-serial/lib -lkokkoscore -ldl -o $(BUILD_DIR)/$(TEST_CASE)-kokkos-serial
	rm $(BUILD_DIR)/$(TEST_CASE)-kokkos-serial.o

$(BUILD_DIR)/$(TEST_CASE)-kokkos-omp-host: $(TEST_CASE)-kokkos.cpp $(TEST_CASE)-util.h ../util.h
	g++ -O3 -march=native -std=c++17 -fopenmp -I${HOME}/kokkos/install-omp/include -c $(TEST_CASE)-kokkos.cpp -o $(BUILD_DIR)/$(TEST_CASE)-kokkos-omp-host.o
	g++ -fopenmp $(BUILD_DIR)/$(TEST_CASE)-kokkos-omp-host.o -L${HOME}/kokkos/install-omp/lib -lkokkoscore -ldl -o $(BUILD_DIR)/$(TEST_CASE)-kokkos-omp-host
	rm $(BUILD_DIR)/$(TEST_CASE)-kokkos-omp-host.o

$(BUILD_DIR)/$(TEST_CASE)-kokkos-cuda: $(TEST_CASE)-kokkos.cpp $(TEST_CASE)-util.h ../util.h
	${HOME}/kokkos/install-cuda/bin/nvcc_wrapper -O3 -std=c++17 -arch=sm_75 --expt-extended-lambda -I${HOME}/kokkos/install-cuda/include -c $(TEST_CASE)-kokkos.cpp -o $(BUILD_DIR)/$(TEST_CASE)-kokkos-cuda.o
	${HOME}/kokkos/install-cuda/bin/nvcc_wrapper -arch=sm_75 $(BUILD_DIR)/$(TEST_CASE)-kokkos-cuda.o -L${HOME}/kokkos/install-cuda/lib -lkokkoscore -ldl -o $(BUILD_DIR)/$(TEST_CASE)-kokkos-cuda
	rm $(BUILD_DIR)/$(TEST_CASE)-kokkos-cuda.o

$(BUILD_DIR)/$(TEST_CASE)-sycl-implicit: $(TEST_CASE)-sycl-implicit.cpp $(TEST_CASE)-util.h ../util.h ../sycl-util.h
	icpx -O3 -march=native -std=c++17 -fsycl -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend --cuda-gpu-arch=sm_75 -o $(BUILD_DIR)/$(TEST_CASE)-sycl-implicit $(TEST_CASE)-sycl-implicit.cpp

$(BUILD_DIR)/$(TEST_CASE)-sycl-explicit: $(TEST_CASE)-sycl-explicit.cpp $(TEST_CASE)-util.h ../util.h ../sycl-util.h
	icpx -O3 -march=native -std=c++17 -fsycl -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend --cuda-gpu-arch=sm_75 -o $(BUILD_DIR)/$(TEST_CASE)-sycl-explicit $(TEST_CASE)-sycl-explicit.cpp

$(BUILD_DIR)/$(TEST_CASE)-hip: $(TEST_CASE)-hip.hip $(TEST_CASE)-util.h ../util.h ../hip-util.h
	hipcc -O3 -std=c++17 -o $(BUILD_DIR)/$(TEST_CASE)-hip $(TEST_CASE)-hip.hip


# aliases without build directory

.PHONY: $(TEST_CASE)-base
$(TEST_CASE)-base: $(BUILD_DIR)/$(TEST_CASE)-base

.PHONY: $(TEST_CASE)-omp-host
$(TEST_CASE)-omp-host: $(BUILD_DIR)/$(TEST_CASE)-omp-host

.PHONY: $(TEST_CASE)-cuda
$(TEST_CASE)-cuda: $(BUILD_DIR)/$(TEST_CASE)-cuda

.PHONY: $(TEST_CASE)-openacc
$(TEST_CASE)-openacc: $(BUILD_DIR)/$(TEST_CASE)-openacc

.PHONY: $(TEST_CASE)-omp-target
$(TEST_CASE)-omp-target: $(BUILD_DIR)/$(TEST_CASE)-omp-target

.PHONY: $(TEST_CASE)-cpp-17
$(TEST_CASE)-cpp-17: $(BUILD_DIR)/$(TEST_CASE)-cpp-17

.PHONY: $(TEST_CASE)-kokkos-serial
$(TEST_CASE)-kokkos-serial: $(BUILD_DIR)/$(TEST_CASE)-kokkos-serial

.PHONY: $(TEST_CASE)-kokkos-omp-host
$(TEST_CASE)-kokkos-omp-host: $(BUILD_DIR)/$(TEST_CASE)-kokkos-omp-host

.PHONY: $(TEST_CASE)-kokkos-cuda
$(TEST_CASE)-kokkos-cuda: $(BUILD_DIR)/$(TEST_CASE)-kokkos-cuda

.PHONY: $(TEST_CASE)-sycl-implicit
$(TEST_CASE)-sycl-implicit: $(BUILD_DIR)/$(TEST_CASE)-sycl-implicit

.PHONY: $(TEST_CASE)-sycl-explicit
$(TEST_CASE)-sycl-explicit: $(BUILD_DIR)/$(TEST_CASE)-sycl-explicit

.PHONY: $(TEST_CASE)-hip
$(TEST_CASE)-hip: $(BUILD_DIR)/$(TEST_CASE)-hip


# automated benchmark target

.PHONY: bench
bench: all
	@echo "Base:"
	$(BUILD_DIR)/$(TEST_CASE)-base $(PARAMETERS)
	#pnmtopng result.pnm > result-base.png
	#rm result.pnm
	#python3 convertToImage.py --input result.raw --output result-base-color.png --width $(NX) --height $(NY)
	@echo ""
#	@echo "OpenMP Host:"
#	$(BUILD_DIR)/$(TEST_CASE)-omp-host $(PARAMETERS)
#	@echo ""
#	@echo "CUDA:"
#	$(BUILD_DIR)/$(TEST_CASE)-cuda $(PARAMETERS)
#	@echo ""
#	@echo "OpenACC:"
#	$(BUILD_DIR)/$(TEST_CASE)-openacc $(PARAMETERS)
#	@echo ""
##	@echo "OpenMP Target:"
##	$(BUILD_DIR)/$(TEST_CASE)-omp-target $(PARAMETERS)
##	@echo ""
#	@echo "CPP-17:"
#	$(BUILD_DIR)/$(TEST_CASE)-cpp-17 $(PARAMETERS)
#	@echo ""
#	@echo "Kokkos serial:"
#	$(BUILD_DIR)/$(TEST_CASE)-kokkos-serial $(PARAMETERS)
#	@echo ""
#	@echo "Kokkos OpenMP Host:"
#	OMP_PROC_BIND=spread OMP_PLACES=threads $(BUILD_DIR)/$(TEST_CASE)-kokkos-omp-host $(PARAMETERS)
#	@echo ""
#	@echo "Kokkos CUDA:"
#	$(BUILD_DIR)/$(TEST_CASE)-kokkos-cuda $(PARAMETERS)
#	@echo ""
#	@echo "SYCL implicit:"
#	ONEAPI_DEVICE_SELECTOR=cuda:* $(BUILD_DIR)/$(TEST_CASE)-sycl-implicit $(PARAMETERS)
#	@echo ""
#	@echo "SYCL explicit:"
#	ONEAPI_DEVICE_SELECTOR=cuda:* $(BUILD_DIR)/$(TEST_CASE)-sycl-explicit $(PARAMETERS)
#	@echo ""
##	@echo "HIP:"
##	$(buildDir)/$(testcase)-hip $(PARAMETERS)
##	@echo ""


# clean target

.PHONY: clean
clean:
	rm $(targets)
