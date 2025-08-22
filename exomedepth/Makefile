BUILD   := $(shell git log -1 --pretty=%h)

# define image names
APP      := exomedepth
REGISTRY := seglh
DIR := $(shell pwd)

# build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(BUILD)
IMG_DEV       := $(IMG):dev

.PHONY: push build version cleanbuild

push: build
	docker push $(IMG_VERSIONED)

build: version
	docker buildx build --platform linux/amd64 -t $(IMG_VERSIONED) . || docker build -t $(IMG_VERSIONED) .
	docker save $(IMG_VERSIONED) | gzip > $(DIR)/$(REGISTRY)-$(APP):$(BUILD).tar.gz

cleanbuild: version
	docker buildx build --platform linux/amd64 --no-cache -t $(IMG_VERSIONED) . || docker build -t $(IMG_VERSIONED) .

devbuild: version
	docker build -t $(IMG_DEV) .

version:
	echo $(BUILD) > VERSION
