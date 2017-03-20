/*
Nexus

Copyright(C) 2012 - Federico Ponchio
ISTI - Italian National Research Council - Visual Computing Lab

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
for more details.
*/
#ifndef NX_CSTREAM_H
#define NX_CSTREAM_H

#include <string.h>
#include <iostream>

#include "bitstream.h"

typedef unsigned char uchar;

namespace nx {


class Stream {
protected:

	uchar *buffer;
	uchar *pos; //for reading.
	int allocated;
	int stopwatch; //used to measure stream partial size.

public:
	enum Entropy { NONE = 0, TUNSTALL = 1, HUFFMAN = 2, ZLIB = 3, LZ4 = 4 };
	Entropy entropy;

	Stream(): buffer(NULL), pos(NULL), allocated(0), stopwatch(0), entropy(TUNSTALL) {}

	Stream(int _size, uchar *_buffer) {
		init(_size, _buffer);
	}

	~Stream() {
		if(allocated)
			delete []buffer;
	}

	int size() { return pos - buffer; }
	uchar *data() { return buffer; }
	void restart() { stopwatch = size(); }
	int elapsed() {
		int e = size() - stopwatch; stopwatch = size();
		return e;
	}

	int  compress(uint32_t size, uchar *data);
	void decompress(std::vector<uchar> &data);
	int  tunstall_compress(unsigned char *data, int size);
	void tunstall_decompress(std::vector<uchar> &data);

#ifdef ENTROPY_TESTS
	int  zlib_compress(uchar *data, int size);
	void zlib_decompress(std::vector<uchar> &data);
	int  lz4_compress(uchar *data, int size);
	void lz4_decompress(std::vector<uchar> &data);
#endif

	void reserve(int reserved) {
		allocated = reserved;
		stopwatch = 0;
		pos = buffer = new uchar[allocated];
	}

	void init(int /*_size*/, uchar *_buffer) {
		buffer = _buffer;
		pos = buffer;
		allocated = 0;
	}

	void rewind() { pos = buffer; }

	template<class T> void write(T c) {
		grow(sizeof(T));
		*(T *)pos = c;
		pos += sizeof(T);
	}

	template<class T> void writeArray(int s, T *c) {
		int bytes = s*sizeof(T);
		push(c, bytes);
	}

	void write(BitStream &stream) {
		stream.flush();
		//padding to 32 bit is needed for javascript reading (which uses int words.), mem needs to be aligned.
		write<int>((int)stream.size);

		int pad = (pos - buffer) & 0x3;
		if(pad != 0)
			pad = 4 - pad;
		grow(pad);
		pos += pad;
		push(stream.buffer, stream.size*sizeof(uint64_t));
	}

	template<class T> T read() {
		T c;
		c = *(T *)pos;
		pos += sizeof(T);
		return c;
	}

	template<class T> T *readArray(int s) {
		int bytes = s*sizeof(T);
		T *buffer = (T *)pos;
		pos += bytes;
		return buffer;
	}

	void read(BitStream &stream) {
		int s = read<int>();
		//padding to 32 bit is needed for javascript reading (which uses int words.), mem needs to be aligned.
		int pad = (pos - buffer) & 0x3;
		if(pad != 0)
			pos += 4 - pad;
		stream.init(s, (uint64_t *)pos);
		pos += s*sizeof(uint64_t);
	}

	void grow(int s) {
		if(allocated == 0)
			reserve(1024);
		int size = pos - buffer;
		if(size + s > allocated) { //needs more spac
			int new_size = allocated*2;
			while(new_size < size + s)
				new_size *= 2;
			uchar *b = new uchar[new_size];
			memcpy(b, buffer, allocated);
			delete []buffer;
			buffer = b;
			pos = buffer + size;
			allocated = new_size;
		}
	}

	void push(void *b, int s) {
		grow(s);
		memcpy(pos, b, s);
		pos += s;
	}


	static int needed(int a) {
		if(a == 0) return 0;
		if(a == -1) return 1;
		if(a < 0) a = -a - 1;
		int n = 2;
		while(a >>= 1) n++;
		return n;
	}

};

} //namespace
#endif // NX_CSTREAM_H
