# configuration

BUILD_DIR = ../build


# default benchmark parameters

NUMLVL      = 10
SUPERSAMPLE = 2
NUMMG       = 10

PARAMETERS = ${NUMLVL} ${SUPERSAMPLE} ${NUMMG}


# all

targets = main


.PHONY: all
all: mk-target-dir $(targets)


mk-target-dir:
	mkdir -p $(BUILD_DIR)


# build rules

$(BUILD_DIR)/main: main.cpp ray-tracing.h ray-tracing-util.h intersect.h vec3.h optical-flow.h
	g++ -O3 -march=native -std=c++17 -o $(BUILD_DIR)/ray-tracing-optical-flow main.cpp  -lm -lpthread -lX11 -ljpeg -lpng


# aliases without build directory

.PHONY: main
main: $(BUILD_DIR)/main


# automated benchmark target

.PHONY: bench
bench: all
	@echo "Base:"
	$(BUILD_DIR)/ray-tracing-optical-flow $(PARAMETERS)
	#pnmtopng result.pnm > result-base.png
	#rm result.pnm
	#python3 convertToImage.py --input result.raw --output result-base-color.png --width $(NX) --height $(NY)
	@echo ""


# clean target

.PHONY: clean
clean:
	rm $(targets)
