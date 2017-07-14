prefix=${HOME}/sf/rnadiff/current

all:
	mkdir -p bin
	mkdir -p doc
	cp rnadiff bin/
	cp rnaseq.R bin/
	cp README.md doc/
	cp LICENSE doc/
	chmod 755 bin/rnadiff
	chmod 644 bin/rnaseq.R
	chmod 644 doc/*
install:
	mkdir -p $(prefix)/bin
	mkdir -p $(prefix)/demo
	mkdir -p $(prefix)/doc
	install -p -m 0755 bin/rnadiff $(prefix)/bin/rnadiff
	install -p -m 0644 bin/rnaseq.R $(prefix)/bin/rnaseq.R
	install -p -m 0644 demo/* $(prefix)/demo/
	install -p -m 0644 doc/* $(prefix)/doc/
clean:
	rm -rf bin
	rm -rf doc
