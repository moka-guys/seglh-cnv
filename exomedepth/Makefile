VERSION  := $(shell git describe --tags --always --dirty)

# define image names
APP      := exomedepth
REGISTRY := seglh

# build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(VERSION)
.PHONY: push build

push: build
	docker push $(IMG_VERSIONED)

build: version
	docker buildx build --platform linux/amd64 -t $(IMG_VERSIONED) . || docker build -t $(IMG_VERSIONED) .
	docker save $(IMG_VERSIONED) | gzip > $(REGISTRY)_$(APP)_$(VERSION).tar.gz

version:
	echo $(VERSION) > VERSION
