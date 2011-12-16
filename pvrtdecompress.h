#import <stdint.h>
 
typedef struct {
    uint32_t filesz;
    uint16_t creator1;
    uint16_t creator2;
    uint32_t bmp_offset;
} BITMAPFILEHEADER;

typedef struct {
    uint32_t header_sz;
    int32_t width;
    int32_t height;
    uint16_t nplanes;
    uint16_t bitspp;
    uint32_t compress_type;
    uint32_t bmp_bytesz;
    int32_t hres;
    int32_t vres;
    uint32_t ncolors;
    uint32_t nimpcolors;
} BITMAPINFOHEADER;

typedef struct {
	uint8_t b;
	uint8_t g;
	uint8_t r;
	uint8_t a;
} RGBQUAD;

typedef struct
{
    uint8_t blue;
    uint8_t green;
    uint8_t red;
} __attribute__((__packed__))
RGBTRIPLE;

/*!***********************************************************************
 @Function		PVRTDecompress
 @Input			pCompressedData The PVRTC texture data to decompress
 @Input			Do2bitMode Signifies whether the data is PVRTC2 or PVRTC4
 @Input			XDim X dimension of the texture
 @Input			YDim Y dimension of the texture
 @Modified		pResultImage The decompressed texture data
 @Description	Decompresses PVRTC to RGBA 8888
*************************************************************************/

void pvrtdecompress(const void *input_buf, const int is_2bpp,
                    const int xmax, const int ymax,
                    unsigned char *result_buf);
