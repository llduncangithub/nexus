#include "nxz.h"

using namespace nx;

static int ilog2(uint64_t p) {
	int k = 0;
	while ( p>>=1 ) { ++k; }
	return k;
}

static int needed(int a) {
	if(a == 0) return 0;
	if(a == -1) return 1;
	if(a < 0) a = -a - 1;
	int n = 2;
	while(a >>= 1) n++;
	return n;
}

//T can be char, short, int
void encodeDiffs(Stream &stream, uint32_t size, int *values) {
	std::vector<uchar> logs(size);
	BitStream bitstream(size);

	for(uint32_t i = 0; i < size; i++) {
		int &val = values[i];
		if(val == 0) {
			logs[i] = 0;
			continue;
		}
		int ret = ilog2(abs(val)) + 1;  //0 -> 0, [1,-1] -> 1 [-2,-3,2,3] -> 2
		logs[i] = ret;
		int middle = (1<<ret)>>1;
		if(val < 0) val = -val -middle;
		bitstream.writeUint(val, ret);
	}
	stream.compress(logs.size(), &*logs.begin());
	stream.write(bitstream);
}

//T can be char, short, int
void encodeDiff(Stream &stream, uint32_t size, short *values) {
	std::vector<uchar> logs(size);
	BitStream bitstream(size);

	for(uint32_t i = 0; i < size; i++) {
		short &val = values[i];
		if(val == 0) {
			logs[i] = 0;
			continue;
		}
		int ret = ilog2(abs(val)) + 1;  //0 -> 0, [1,-1] -> 1 [-2,-3,2,3] -> 2
		logs[i] = ret;
		int middle = (1<<ret)>>1;
		if(val < 0) val = -val -middle;
		bitstream.writeUint(val, ret);
	}
	stream.compress(logs.size(), &*logs.begin());
	stream.write(bitstream);
}


void encodeDiff(Stream &stream, uint32_t size, Point2i *values) {
	std::vector<uchar> logs(size);
	BitStream bitstream(size);

	for(uint32_t i = 0; i < size; i++) {
		Point2i p = values[i];

		int diff = std::max(needed(p[0]), needed(p[1]));
		logs[i] = diff;
		if(diff == 0) continue;

		int max = 1<<(diff-1);
		bitstream.writeUint(p[0] + max, diff);
		bitstream.writeUint(p[1] + max, diff);
	}
	stream.compress(logs.size(), &*logs.begin());
	stream.write(bitstream);
}

void encodeDiff(Stream &stream, uint32_t size, Point2s *values) {
	std::vector<uchar> logs(size);
	BitStream bitstream(size);

	for(uint32_t i = 0; i < size; i++) {
		Point2s &p = values[i];

		int diff = std::max(needed(p[0]), needed(p[1]));
		logs[i] = diff;
		if(diff == 0) continue;

		int max = 1<<(diff-1);
		bitstream.writeUint(p[0] + max, diff);
		bitstream.writeUint(p[1] + max, diff);
	}
	stream.compress(logs.size(), &*logs.begin());
	stream.write(bitstream);
}

void encodeDiff(Stream &stream, uint32_t size, Point3i *values) {
	std::vector<uchar> logs(size);
	BitStream bitstream(size);

	for(uint32_t i = 0; i < size; i++) {
		Point3i &p = values[i];
		int diff = std::max(std::max(needed(p[0]), needed(p[1])), needed(p[2]));
		logs[i] = diff;
		if(diff == 0) continue;

		int max = 1<<(diff-1);
		bitstream.writeUint(p[0] + max, diff);
		bitstream.writeUint(p[1] + max, diff);
		bitstream.writeUint(p[2] + max, diff);
	}
	stream.compress(logs.size(), &*logs.begin());
	stream.write(bitstream);
}



const uint64_t bmask[] = {
	0x00,     0x01,       0x03,        0x07,        0x0f,       0x01f,       0x03f,       0x07f,
	0xff,     0x01ff,     0x03ff,      0x07ff,      0x0fff,     0x01fff,     0x03fff,     0x07fff,
	0xffff,   0x01ffff,   0x03ffff,    0x07ffff,    0x0fffff,   0x01fffff,   0x03fffff,   0x07fffff,
	0xffffff, 0x01ffffff, 0x03ffffff,  0x07ffffff,  0x0fffffff, 0x01fffffff, 0x03fffffff, 0x7fffffff,

	0xffffffff,       0x01ffffffff,       0x03ffffffff,        0x07ffffffff,        0x0fffffffff,       0x01fffffffff,       0x03fffffffff,       0x07fffffffff,
	0xffffffffff,     0x01ffffffffff,     0x03ffffffffff,      0x07ffffffffff,      0x0fffffffffff,     0x01fffffffffff,     0x03fffffffffff,     0x07fffffffffff,
	0xffffffffffff,   0x01ffffffffffff,   0x03ffffffffffff,    0x07ffffffffffff,    0x0fffffffffffff,   0x01fffffffffffff,   0x03fffffffffffff,   0x07fffffffffffff,
	0xffffffffffffff, 0x01ffffffffffffff, 0x03ffffffffffffff,  0x07ffffffffffffff,  0x0fffffffffffffff, 0x01fffffffffffffff, 0x03fffffffffffffff, 0x07fffffffffffffff,

	0xffffffffffffffff };

const uint64_t bmax[] = { 0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
						  1<<11, 1<<12, 1<<13, 1<<14, 1<<15, 1<<16, 1<<17, 1<<18, 1<<19, 1<< 20,
						  1<<21, 1<<22, 1<<23, 1<<24, 1<<25, 1<<26, 1<<27, 1<<28, 1<<29, 1<<30 };


int decodeDiff(Stream &stream, int *values) {
	std::vector<uchar> diffs;
	stream.decompress(diffs);
	BitStream bitstream;
	stream.read(bitstream);

	for(uint32_t i = 0; i < diffs.size(); i++) {
		uchar &diff = diffs[i];
		if(diff == 0) {
			values[i] = 0;
			continue;
		}

		int val = bitstream.readUint(diff);
		int middle = 1<<(diff-1);
		if(val < middle)
			val = -val -middle;
		values[i] = val;
	}
	return diffs.size();
}

int decodeDiff(Stream &stream, Point2i *values) {
	std::vector<uchar> diffs;
	stream.decompress(diffs);
	BitStream bitstream;
	stream.read(bitstream);

	for(uint32_t i =0; i < diffs.size(); i++) {
		Point2i &p = values[i];
		uchar &diff = diffs[i];
		if(diff == 0) {
			p[0] = p[1] =0;
			continue;
		}

		const uint64_t &max = bmax[diff];
//		uint64_t max = (1<<diff)>>1;
		if(diff < 22) {
			const int64_t &mask = bmask[diff];
			//uint64_t mask = (1<<diff)-1;
			uint64_t bits = bitstream.readUint(2*diff);
			p[1] = (bits & mask) - max;
			bits >>= diff;
			p[0] = bits - max;
		} else {
			p[0] = bitstream.readUint(diff) - max;
			p[1] = bitstream.readUint(diff) - max;
		}
	}
	return diffs.size();
}

int decodeDiff(Stream &stream, Point2s *values) {
	std::vector<uchar> diffs;
	stream.decompress(diffs);
	BitStream bitstream;
	stream.read(bitstream);

	for(uint32_t i =0; i < diffs.size(); i++) {
		Point2s &p = values[i];
		uchar &diff = diffs[i];
		if(diff == 0) {
			p[0] = p[1] =0;
			continue;
		}

		const uint64_t &max = bmax[diff];
//		uint64_t max = (1<<diff)>>1;
		if(diff < 22) {
			const uint64_t &mask = bmask[diff];
			//uint64_t mask = (1<<diff)-1;
			uint64_t bits = bitstream.readUint(2*diff);
			p[1] = (bits & mask) - max;
			bits >>= diff;
			p[0] = bits - max;
		} else {
			p[0] = bitstream.readUint(diff) - max;
			p[1] = bitstream.readUint(diff) - max;
		}
	}
	return diffs.size();
}


int decodeDiff(Stream &stream, Point3i *values) {
	std::vector<uchar> diffs;
	stream.decompress(diffs);
	BitStream bitstream;
	stream.read(bitstream);

	for(uint32_t i =0; i < diffs.size(); i++) {
		Point3i &p = values[i];
		uchar &diff = diffs[i];
		if(diff == 0) {
			p[0] = p[1] = p[2] = 0;
			continue;
		}
		//making a single read is 2/3 faster
		//uint64_t &max = bmax[diff];
		const uint64_t max = (1<<diff)>>1;
		if(diff < 22) {
			//uint64_t &mask = bmask[diff]; //using table is 4% faster
			const uint64_t mask = (1<<diff)-1;
			uint64_t bits = bitstream.readUint(3*diff);
			p[2] = (bits & mask) - max;
			bits >>= diff;
			p[1] = (bits & mask) - max;
			bits >>= diff;
			p[0] = bits - max;
		} else {
			p[0] = bitstream.readUint(diff) - max;
			p[1] = bitstream.readUint(diff) - max;
			p[2] = bitstream.readUint(diff) - max;
		}
	}
	return diffs.size();
}

int decodeDiff(Stream &stream, Point3s *values) {
	std::vector<uchar> diffs;
	stream.decompress(diffs);
	BitStream bitstream;
	stream.read(bitstream);

	for(uint32_t i =0; i < diffs.size(); i++) {
		Point3s &p = values[i];
		uchar &diff = diffs[i];
		if(diff == 0) {
			p[0] = p[1] = p[2] = 0;
			continue;
		}
		//making a single read is 2/3 faster
		//uint64_t &max = bmax[diff];
		const uint64_t max = (1<<diff)>>1;
		if(diff < 22) {
			//uint64_t &mask = bmask[diff]; //using table is 4% faster
			const uint64_t mask = (1<<diff)-1;
			uint64_t bits = bitstream.readUint(3*diff);
			p[2] = (bits & mask) - max;
			bits >>= diff;
			p[1] = (bits & mask) - max;
			bits >>= diff;
			p[0] = bits - max;
		} else {
			p[0] = bitstream.readUint(diff) - max;
			p[1] = bitstream.readUint(diff) - max;
			p[2] = bitstream.readUint(diff) - max;
		}
	}
	return diffs.size();
}

