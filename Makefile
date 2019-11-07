all:mrcs2jpeg
mrcs2jpeg:mrcs2jpeg.c
	g++ -Ijpeg/include -o $@ $^ jpeg/lib/libjpeg.a -lm
.phony:clean
clean:
	rm -rf mrcs2jpeg
