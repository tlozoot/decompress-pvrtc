/*
 * fortseige_pvr.c
 * Adapted from https://github.com/abutcher/FortSiege/blob/\
 *      b04f5738c8d8991b7e98af21dac956cf1fdd8d79/FortSiege/\
 *      cocos3d/cc3PVR/PVRT%202.07/PVRTDecompress.cpp
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <assert.h>
 
#define PT_INDEX    2  // The Punch-through index
#define BLK_Y_SIZE  4  // always 4 for all 2D block types
#define BLK_X_MAX   8  // Max X dimension for blocks
#define BLK_X_2BPP  8  // dimensions for the two formats
#define BLK_X_4BPP  4  

#define WRAP_COORD(val, size) ((val) & ((size) - 1))
#define PVRT_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define PVRT_MAX(a ,b) (((a) > (b)) ? (a) : (b))
#define PVRT_CLAMP(x, l, h) (PVRT_MIN((h), PVRT_MAX((x), (l))))

#define POWER_OF_2(x) number_is_power_2(x)

// Wrap or clamp large or small vals to the legal coordinate range
#define LIMIT_COORD(val, size, tiles) ((tiles) ? WRAP_COORD((val), (size)) : PVRT_CLAMP((val), 0, (size) - 1))

typedef struct
{
	// Uses 64 bits per block
	uint32_t PackedData[2];
} AMTC_BLOCK_STRUCT;

/***********************************************************
				DECOMPRESSION ROUTINES
************************************************************/

// Returns the n left bits of x and shifts x by n bits
uint8_t shiftl_16(uint16_t *x, int n)
{
    uint8_t result = *x >> (16 - n);
    *x <<= n;
    return result;
}

uint8_t replicate_top_bit(uint8_t x)
{
    return x | x >> 4;
}


 /*!***********************************************************************
  @Function		number_is_power_2
  @Input		input A number
  @Returns		1 if the number is an integer power of two, else FALSE.
  @Description	Check that a number is an integer power of two, i.e.
             	1, 2, 4, 8, ... etc.
				Returns 0 for zero.
*************************************************************************/
int number_is_power_2(unsigned input)
{
    if (!input) return 0;
    
    unsigned minus1 = input - 1;
    return ((input | minus1) == (input ^ minus1)) ? 1 : 0;
}


/*!***********************************************************************
 @Function		Unpack5554Colour
 @Input			pBlock
 @Input			ABColours
 @Description	Given a block, extract the colour information and convert
 				to 5554 formats
*************************************************************************/
static void Unpack5554Colour(const AMTC_BLOCK_STRUCT *pBlock,
							 int   ABColours[2][4])
{    
    uint16_t raw_bits[2];
    // Extract A and B
    
    raw_bits[0] = pBlock->PackedData[1] & (0xFFFE); // 15 bits (shifted up by one)
    raw_bits[1] = pBlock->PackedData[1] >> 16;       // 16 bits
        
    // step through both colours
    for (int i = 0; i < 2; i++)
    {
        uint16_t raw_pixel = raw_bits[i];
        
        int is_opaque = shiftl_16(&raw_pixel, 1);
    
        // If completely opaque
        if (is_opaque) {
            ABColours[i][0] = shiftl_16(&raw_pixel, 5);
            ABColours[i][1] = shiftl_16(&raw_pixel, 5);            
            ABColours[i][2] = (i == 0) ? replicate_top_bit(shiftl_16(&raw_pixel, 5))
                                       : shiftl_16(&raw_pixel, 5);          
            ABColours[i][3] = 0xF;
        }
    
        else {
            ABColours[i][3] = shiftl_16(&raw_pixel, 3) << 1;
            ABColours[i][0] = replicate_top_bit(shiftl_16(&raw_pixel, 4) << 1);
            ABColours[i][1] = replicate_top_bit(shiftl_16(&raw_pixel, 4) << 1);

            ABColours[i][2] = shiftl_16(&raw_pixel, 4) << 1;
            if (i == 0)
                ABColours[0][2] |= ABColours[0][2] >> 3;
            else
                ABColours[0][2] = replicate_top_bit(ABColours[0][2]);
    
        }
    }
}

/*!***********************************************************************
 @Function		UnpackModulations
 @Input			pBlock
 @Input			Do2bitMode
 @Input			ModulationVals
 @Input			ModulationModes
 @Input			StartX
 @Input			StartY
 @Description	Given the block and the texture type and it's relative
 				position in the 2x2 group of blocks, extract the bit
 				patterns for the fully defined pixels.
*************************************************************************/
static void	UnpackModulations(const AMTC_BLOCK_STRUCT *pBlock,
							  const int Do2bitMode,
							  int ModulationVals[8][16],
							  int ModulationModes[8][16],
							  int StartX,
							  int StartY)
{
	int BlockModMode;
	uint32_t ModulationBits;

	int x, y;

	BlockModMode= pBlock->PackedData[1] & 1;
	ModulationBits	= pBlock->PackedData[0];

	// if it's in an interpolated mode
	if(Do2bitMode && BlockModMode)
	{
		/*
			run through all the pixels in the block. Note we can now treat all the
			"stored" values as if they have 2bits (even when they didn't!)
		*/
		for(y = 0; y < BLK_Y_SIZE; y++)
		{
			for(x = 0; x < BLK_X_2BPP; x++)
			{
				ModulationModes[y+StartY][x+StartX] = BlockModMode;

				// if this is a stored value...
				if(((x^y)&1) == 0)
				{
					ModulationVals[y+StartY][x+StartX] = ModulationBits & 3;
					ModulationBits >>= 2;
				}
			}
		}
	}
	else if(Do2bitMode) // else if direct encoded 2bit mode - i.e. 1 mode bit per pixel
	{
		for(y = 0; y < BLK_Y_SIZE; y++)
		{
			for(x = 0; x < BLK_X_2BPP; x++)
			{
				ModulationModes[y+StartY][x+StartX] = BlockModMode;

				// double the bits so 0=> 00, and 1=>11
				if(ModulationBits & 1)
				{
					ModulationVals[y+StartY][x+StartX] = 0x3;
				}
				else
				{
					ModulationVals[y+StartY][x+StartX] = 0x0;
				}
				ModulationBits >>= 1;
			}
		}
	}
	else // else its the 4bpp mode so each value has 2 bits
	{
		for(y = 0; y < BLK_Y_SIZE; y++)
		{
			for(x = 0; x < BLK_X_4BPP; x++)
			{
				ModulationModes[y+StartY][x+StartX] = BlockModMode;

				ModulationVals[y+StartY][x+StartX] = ModulationBits & 3;
				ModulationBits >>= 2;
			}
		}
	}

	// make sure nothing is left over
	assert((ModulationBits == 0));
}

/*!***********************************************************************
 @Function		InterpolateColours
 @Input			ColourP
 @Input			ColourQ
 @Input			ColourR
 @Input			ColourS
 @Input			Do2bitMode
 @Input			x
 @Input			y
 @Modified		result
 @Description	This performs a HW bit accurate interpolation of either the
				A or B colours for a particular pixel.

				NOTE: It is assumed that the source colours are in ARGB 5554
				format - This means that some "preparation" of the values will
				be necessary.
*************************************************************************/
static void InterpolateColours(const int ColourP[4],
						  const int ColourQ[4],
						  const int ColourR[4],
						  const int ColourS[4],
						  const int Do2bitMode,
						  const int x,
						  const int y,
						  int result[4])
{
	int u, v, uscale;
	int k;

	int tmp1, tmp2;

	int P[4], Q[4], R[4], S[4];

	// Copy the colours
	for (k = 0; k < 4; k++)
	{
		P[k] = ColourP[k];
		Q[k] = ColourQ[k];
		R[k] = ColourR[k];
		S[k] = ColourS[k];
	}

	// put the x and y values into the right range
	v = (y & 0x3) | ((~y & 0x2) << 1);

	if (Do2bitMode)
		u = (x & 0x7) | ((~x & 0x4) << 1);
	else
		u = (x & 0x3) | ((~x & 0x2) << 1);

	// get the u and v scale amounts
	v  = v - BLK_Y_SIZE/2;

	if (Do2bitMode) {
		u = u - BLK_X_2BPP/2;
		uscale = 8;
	}
	 
	else {
		u = u - BLK_X_4BPP/2;
		uscale = 4;
	}

	for (k = 0; k < 4; k++) {
		tmp1 = P[k] * uscale + u * (Q[k] - P[k]);
		tmp2 = R[k] * uscale + u * (S[k] - R[k]);

		tmp1 = tmp1 * 4 + v * (tmp2 - tmp1);

		result[k] = tmp1;
	}

	// Lop off the appropriate number of bits to get us to 8 bit precision
	if (Do2bitMode) {
		// do RGB
		for (k = 0; k < 3; k++) {
			result[k] >>= 2;
		}

		result[3] >>= 1;
	}
	
	else {
		// do RGB  (A is ok)
		for (k = 0; k < 3; k++) {
			result[k] >>= 1;
		}
	}

	// sanity check
	for (k = 0; k < 4; k++) {
		assert(result[k] < 256);
	}


	/*
		Convert from 5554 to 8888

		do RGB 5.3 => 8
	*/
	
	for (k = 0; k < 3; k++) {
		result[k] += result[k] >> 5;
	}

	result[3] += result[3] >> 4;

	// 2nd sanity check
	for (k = 0; k < 4; k++) {
		assert(result[k] < 256);
	}

}

/*!***********************************************************************
 @Function		GetModulationValue
 @Input			x
 @Input			y
 @Input			Do2bitMode
 @Input			ModulationVals
 @Input			ModulationModes
 @Input			Mod
 @Input			DoPT
 @Description	Get the modulation value as a numerator of a fraction of 8ths
*************************************************************************/
static void GetModulationValue(int x,
							   int y,
							   const int Do2bitMode,
							   const int ModulationVals[8][16],
							   const int ModulationModes[8][16],
							   int *Mod,
							   int *DoPT)
{
	static const int RepVals0[4] = {0, 3, 5, 8};
	static const int RepVals1[4] = {0, 4, 4, 8};

	int ModVal;

	// Map X and Y into the local 2x2 block
	y = (y & 0x3) | ((~y & 0x2) << 1);

	if(Do2bitMode)
		x = (x & 0x7) | ((~x & 0x4) << 1);
	else
		x = (x & 0x3) | ((~x & 0x2) << 1);

	// assume no PT for now
	*DoPT = 0;

	// extract the modulation value. If a simple encoding
	if(ModulationModes[y][x]==0)
	{
		ModVal = RepVals0[ModulationVals[y][x]];
	}
	else if(Do2bitMode)
	{
		// if this is a stored value
		if(((x^y)&1)==0)
			ModVal = RepVals0[ModulationVals[y][x]];
		else if(ModulationModes[y][x] == 1) // else average from the neighbours if H&V interpolation..
		{
			ModVal = (RepVals0[ModulationVals[y-1][x]] +
					  RepVals0[ModulationVals[y+1][x]] +
					  RepVals0[ModulationVals[y][x-1]] +
					  RepVals0[ModulationVals[y][x+1]] + 2) / 4;
		}
		else if(ModulationModes[y][x] == 2) // else if H-Only
		{
			ModVal = (RepVals0[ModulationVals[y][x-1]] +
					  RepVals0[ModulationVals[y][x+1]] + 1) / 2;
		}
		else // else it's V-Only
		{
			ModVal = (RepVals0[ModulationVals[y-1][x]] +
					  RepVals0[ModulationVals[y+1][x]] + 1) / 2;
		}
	}
	else // else it's 4BPP and PT encoding
	{
		ModVal = RepVals1[ModulationVals[y][x]];

		*DoPT = ModulationVals[y][x] == PT_INDEX;
	}

	*Mod = ModVal;
}

/*!***********************************************************************
 @Function		TwiddleUV
 @Input			YSize	Y dimension of the texture in pixels
 @Input			XSize	X dimension of the texture in pixels
 @Input			YPos	Pixel Y position
 @Input			XPos	Pixel X position
 @Returns		The twiddled offset of the pixel
 @Description	Given the Block (or pixel) coordinates and the dimension of
 				the texture in blocks (or pixels) this returns the twiddled
 				offset of the block (or pixel) from the start of the map.

				NOTE the dimensions of the texture must be a power of 2
*************************************************************************/
static int DisableTwiddlingRoutine = 0;

static uint32_t TwiddleUV(uint32_t YSize, uint32_t XSize, uint32_t YPos, uint32_t XPos)
{
	uint32_t Twiddled;

	uint32_t MinDimension;
	uint32_t MaxValue;

	uint32_t SrcBitPos;
	uint32_t DstBitPos;

	int ShiftCount;

	assert(YPos < YSize);
	assert(XPos < XSize);

	assert(POWER_OF_2(YSize));
	assert(POWER_OF_2(XSize));

	if (YSize < XSize) {
		MinDimension = YSize;
		MaxValue	 = XPos;
	} else {
		MinDimension = XSize;
		MaxValue	 = YPos;
	}

	// Nasty hack to disable twiddling
	if (DisableTwiddlingRoutine)
		return (YPos * XSize + XPos);

	// Step through all the bits in the "minimum" dimension
	SrcBitPos = 1;
	DstBitPos = 1;
	Twiddled = 0;
	ShiftCount = 0;

	while (SrcBitPos < MinDimension)
	{
		if (YPos & SrcBitPos) {
			Twiddled |= DstBitPos;
		}

		if (XPos & SrcBitPos) {
			Twiddled |= (DstBitPos << 1);
		}


		SrcBitPos <<= 1;
		DstBitPos <<= 2;
		ShiftCount += 1;

	}

	// prepend any unused bits
	MaxValue >>= ShiftCount;

	Twiddled |= (MaxValue << (2*ShiftCount));

	return Twiddled;
}

/*!***********************************************************************
 @Function		PVRTDecompress
 @Input			input_buf The PVRTC texture data to decompress
 @Input			Do2BitMode Signifies whether the data is PVRTC2 or PVRTC4
 @Input			x_dim X dimension of the texture
 @Input			y_dim Y dimension of the texture
 @Input			AssumeImageTiles Assume the texture data tiles
 @Modified		result_bufzThe decompressed texture data
 @Description	Decompresses PVRTC to RGBA 8888
*************************************************************************/
void pvrtdecompress(AMTC_BLOCK_STRUCT *input_buf, const int is_2bpp,
                    const int x_dim, const int y_dim,
                    unsigned char *result_buf)

{       
    int BLK_X_SIZE = (is_2bpp) ? BLK_X_2BPP : BLK_X_4BPP;
    int AssumeImageTiles = 1;
    
	int BlkX, BlkY;
	int BlkXp1, BlkYp1;
	int BlkXDim, BlkYDim;

	int StartX, StartY;

	int ModulationVals[8][16];
	int ModulationModes[8][16];

	int Mod, DoPT;
	
	// local neighbourhood of blocks
	AMTC_BLOCK_STRUCT *pBlocks[2][2];

	AMTC_BLOCK_STRUCT *pPrevious[2][2] = {{NULL, NULL}, {NULL, NULL}};

	// Low precision colours extracted from the blocks
	struct
	{
		int Reps[2][4];
	} Colours5554[2][2];

	// Interpolated A and B colours for the pixel
	int a_color[4], b_color[4];

	int result[4];

	// For MBX don't allow the sizes to get too small
	BlkXDim = PVRT_MAX(2, x_dim / BLK_X_SIZE);
	BlkYDim = PVRT_MAX(2, y_dim / BLK_Y_SIZE);

	/*
		Step through the pixels of the image decompressing each one in turn

		Note that this is a hideously inefficient way to do this!
	*/
	for (int y = 0; y < y_dim; y++)
	{
		for (int x = 0; x < x_dim; x++)
		{
			// map this pixel to the top left neighbourhood of blocks
			BlkX = (x - BLK_X_SIZE / 2);
			BlkY = (y - BLK_Y_SIZE / 2);

			BlkX = LIMIT_COORD(BlkX, x_dim, AssumeImageTiles);
			BlkY = LIMIT_COORD(BlkY, y_dim, AssumeImageTiles);

			BlkX /= BLK_X_SIZE;
			BlkY /= BLK_Y_SIZE;

			// compute the positions of the other 3 blocks
			BlkXp1 = LIMIT_COORD(BlkX + 1, BlkXDim, AssumeImageTiles);
			BlkYp1 = LIMIT_COORD(BlkY + 1, BlkYDim, AssumeImageTiles);

			// Map to block memory locations
			pBlocks[0][0] = input_buf + TwiddleUV(BlkYDim, BlkXDim, BlkY, BlkX);
			pBlocks[0][1] = input_buf + TwiddleUV(BlkYDim, BlkXDim, BlkY, BlkXp1);
			pBlocks[1][0] = input_buf + TwiddleUV(BlkYDim, BlkXDim, BlkYp1, BlkX);
			pBlocks[1][1] = input_buf + TwiddleUV(BlkYDim, BlkXDim, BlkYp1, BlkXp1);


			/*
				extract the colours and the modulation information IF the previous values
				have changed.
			*/
			if (memcmp(pPrevious, pBlocks, 4 * sizeof(void *)) != 0)
			{
				StartY = 0;
				for (int i = 0; i < 2; i++)
				{
					StartX = 0;
					for (int j = 0; j < 2; j++)
					{
						Unpack5554Colour(pBlocks[i][j], Colours5554[i][j].Reps);

						UnpackModulations(pBlocks[i][j],
										  is_2bpp,
										  ModulationVals,
										  ModulationModes,
										  StartX, StartY);

						StartX += BLK_X_SIZE;
					}

					StartY += BLK_Y_SIZE;
				}

				// make a copy of the new pointers
				memcpy(pPrevious, pBlocks, 4 * sizeof(void *));
			}

			// decompress the pixel.  First compute the interpolated A and B signals
			InterpolateColours(Colours5554[0][0].Reps[0],
							   Colours5554[0][1].Reps[0],
							   Colours5554[1][0].Reps[0],
							   Colours5554[1][1].Reps[0],
							   is_2bpp, x, y,
							   a_color);

			InterpolateColours(Colours5554[0][0].Reps[1],
							   Colours5554[0][1].Reps[1],
							   Colours5554[1][0].Reps[1],
							   Colours5554[1][1].Reps[1],
							   is_2bpp, x, y,
							   b_color);

			GetModulationValue(x, y, is_2bpp, (const int (*)[16])ModulationVals, (const int (*)[16])ModulationModes,
							   &Mod, &DoPT);

			// compute the modulated colour
			for (int i = 0; i < 4; i++) {
				result[i] = a_color[i] + (Mod * (b_color[i] - a_color[i]) >> 3);
			}

			if (DoPT)
				result[3] = 0;

			// Store the result in the output image
			unsigned int pos = (x + (y_dim - y - 1) * x_dim ) << 2;
			result_buf[pos + 0] = (uint8_t) result[2];
			result_buf[pos + 1] = (uint8_t) result[1];
			result_buf[pos + 2] = (uint8_t) result[0];
			result_buf[pos + 3] = (uint8_t) result[3];
		}
	}
}