BUILD   := $(shell git log -1 --pretty=%h)

# define image names
APP      := exomedepth
REGISTRY := seglh

# build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(BUILD)
IMG_LATEST    := $(IMG):latest
IMG_DEV       := $(IMG):dev

.PHONY: push build version cleanbuild

push: build
	docker push $(IMG_VERSIONED)
	docker push $(IMG_LATEST)

build: version
	docker buildx build --platform linux/amd64 -t $(IMG_VERSIONED) . || docker build -t $(IMG_VERSIONED) .
	docker tag $(IMG_VERSIONED) $(IMG_LATEST)

devbuild: version
	docker build -t $(IMG_DEV) .

version:
	echo $(BUILD) > VERSION

cleanbuild:
	docker buildx build --platform linux/amd64 --no-cache -t $(IMG_VERSIONED) .
