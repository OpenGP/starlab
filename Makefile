help:
	@echo "This has only been tested for MACOSX!!!!!"
	@echo "- help 			this help"
	@echo "- reinitialize	***deletes*** all sub-repo folders and the .user files"
	@echo "- fresh			reinitialize + download + open"
	@echo "- download 		downloads all sub-repositories"
	@echo "- open	 		opens the starlab_mini.pro in QtCreator"
	@echo "- clean			deletes the deployed application folder e.g. /Applications/Starlab.app"
	@echo "- mini 			builds starlab_mini.pro from sources into ../starlab_mini_build"
	@echo "- full 			builds starlab_full.pro from sources into ../starlab_full_build"

initialize:  reinitialize download open

# CROSS PLATFORM VARIBLE DEFINITIONS
ifeq ($(OS),Windows_NT)
  SPEC=
  OPEN=explorer.exe
  DELETEREC=
else  
  SPEC=macx-g++
  OPEN=open
  DELETEREC=rm -rf
  MKDIR=mkdir -p
endif

MINI_BUILD = ../starlab_mini_build
FULL_BUILD = ../starlab_full_build

download: 
	git clone https://code.google.com/p/starlab.core core
	git clone https://code.google.com/p/starlab.surfacemesh surfacemesh
	git clone https://code.google.com/p/starlab.examples examples

reinitialize:
	$(DELETEREC) /Applications/Starlab.app
	$(DELETEREC) ./core
	$(DELETEREC) ./surfacemesh
	$(DELETEREC) ./examples
	$(DELETEREC) *.pro.user
	$(DELETEREC) $(MINI_BUILD)
	$(DELETEREC) $(FULL_BUILD)

open:
	$(OPEN) starlab_mini.pro

mini:
	$(MKDIR) $(MINI_BUILD)
	cd $(MINI_BUILD); qmake -spec $(SPEC) ../starlab/starlab_mini.pro
	cd $(MINI_BUILD); make -j8

full:
	$(MKDIR) $(FULL_BUILD)
	cd $(FULL_BUILD); qmake -spec $(SPEC) ../starlab/starlab_full.pro
	cd $(FULL_BUILD); make -j8

clean:
	$(DELETEREC) /Applications/Starlab.app
	$(DELETEREC) ../starlab_*_build