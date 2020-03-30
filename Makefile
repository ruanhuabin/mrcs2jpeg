all:mrcs2jpeg cDisp_3DMap
OBJ=mrcs2jpeg cDisp_3DMap
mrcs2jpeg:mrcs2jpeg.c
	g++ -Ijpeg/include -o $@ $^ jpeg/lib/libjpeg.a -lm
cDisp_3DMap:cDisp_3DMap.c mrc.h mrc.cpp readsection.h readsection.cpp
	g++ -Ijpeg/include -o $@ $^ jpeg/lib/libjpeg.a -lm
.phony:clean
clean:
	rm -rf $(OBJ)
