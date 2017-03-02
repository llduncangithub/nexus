NexusCoder

	float error
	int coord_q;
	int norm_q;
	int coord_q;
	int tex_q;

NexusEncoder(



NexusCoder(nvert, nface);

coder.addCoords(float *buffer, int bits); //3x
coder.addColors(uchar *buffer, int luma_bits, int croma_bits, int alpha_bits); //4x
coder.addNormals(float *buffer, int bits); //3x
coder.addTex(float *buffer, int bits); //2x
coder.addData(float *buffer, int bits); //1x

coder.addIndex(uint32_t *buff);
coder.addIndex(uint16_t *buff);
coder.addGroups(int len, uint32_t *buffer); //30, 120, 255; //split into groups the index 0-29, 30-119, in TRIANGLES, last should nvert.
//external mechanism for groups passing

decoder.encode();



NexusDecoder decoder(len, uchar *buffer); //decodes first part of the mesh
int nvert;
int nface;
bool hasindex, normals etc.

decoder.setCoords(float *buffer); //nvert*3 expected
decoder.setColors(uchar *buffer);
decoder.setNormals(float *buffer);
decoder.setNormals(int16_t *buffer);
decoder.setTex(float *buffer);
decoder.setIndex(uint32_t *buffer);
decoder.setIndex(uint16_t *buffer);


decoder.decode();

float *getCoords(); //in case setcoords not enabled
float *getNormals();




