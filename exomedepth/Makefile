BUILD   := $(shell git log -1 --pretty=%h)

# define image names
APP      := exomedepth
REGISTRY := seglh

# build tags
IMG           := $(REGISTRY)/$(APP)
IMG_VERSIONED := $(IMG):$(BUILD)
IMG_LATEST    := $(IMG):latest

push: build
	docker push $(IMG_VERSIONED)
	docker push $(IMG_LATEST)

build:
	docker build -t $(IMG_VERSIONED) .
	docker tag $(IMG_VERSIONED) $(IMG_LATEST)

cleanbuild:
	docker build --no-cache -t $(IMG_VERSIONED) .
