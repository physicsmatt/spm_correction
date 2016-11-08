// This file contains an example of reading a wave file.

#if WIN32
#define _CRT_SECURE_NO_DEPRECATE	// This will prevent the VC2005 compiler from warning you about safe string functions
#endif // WIN32

#include <ctype.h>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include "CrossPlatformFileIO.h"
#include "ReadWave.h"


static void
ReorderBytes(void *p, int bytesPerPoint, long numValues)	// Reverses byte order.
{
	unsigned char ch, *p1, *p2, *pEnd;
	
	pEnd = (unsigned char *)p + numValues*bytesPerPoint;
	while (p < (void *)pEnd) {
		p1 = (unsigned char *)p;
		p2 = (unsigned char *)p + bytesPerPoint-1;
		while (p1 < p2) {
			ch = *p1;
			*p1++ = *p2;
			*p2-- = ch;
		}
		p = (unsigned char *)p + bytesPerPoint;
	}
}

static void
ReorderShort(void* sp)
{
	ReorderBytes(sp, 2, 1);
}

static void
ReorderLong(void* lp)
{
	ReorderBytes(lp, 4, 1);
}

static void
ReorderDouble(void* dp)
{
	ReorderBytes(dp, 8, 1);
}

static void
ReorderBinHeader1(BinHeader1* p)
{
	ReorderShort(&p->version);
	ReorderLong(&p->wfmSize);
	ReorderShort(&p->checksum);
}

static void
ReorderBinHeader2(BinHeader2* p)
{
	ReorderShort(&p->version);
	ReorderLong(&p->wfmSize);
	ReorderLong(&p->noteSize);
	ReorderLong(&p->pictSize);
	ReorderShort(&p->checksum);
}

static void
ReorderBinHeader3(BinHeader3* p)
{
	ReorderShort(&p->version);
	ReorderLong(&p->wfmSize);
	ReorderLong(&p->noteSize);
	ReorderLong(&p->formulaSize);
	ReorderLong(&p->pictSize);
	ReorderShort(&p->checksum);
}

static void
ReorderBinHeader5(BinHeader5* p)
{
	ReorderShort(&p->version);
	ReorderShort(&p->checksum);
	ReorderLong(&p->wfmSize);
	ReorderLong(&p->formulaSize);
	ReorderLong(&p->noteSize);
	ReorderLong(&p->dataEUnitsSize);
	ReorderBytes(&p->dimEUnitsSize, 4, 4);
	ReorderBytes(&p->dimLabelsSize, 4, 4);
	ReorderLong(&p->sIndicesSize);
	ReorderLong(&p->optionsSize1);
	ReorderLong(&p->optionsSize2);
}

static void
ReorderWaveHeader2(WaveHeader2* p)
{
	ReorderShort(&p->type);
	ReorderLong(&p->next);
	// char bname does not need to be reordered.
	ReorderShort(&p->whVersion);
	ReorderShort(&p->srcFldr);
	ReorderLong(&p->fileName);
	// char dataUnits does not need to be reordered.
	// char xUnits does not need to be reordered.
	ReorderLong(&p->npnts);
	ReorderShort(&p->aModified);
	ReorderDouble(&p->hsA);
	ReorderDouble(&p->hsB);
	ReorderShort(&p->wModified);
	ReorderShort(&p->swModified);
	ReorderShort(&p->fsValid);
	ReorderDouble(&p->topFullScale);
	ReorderDouble(&p->botFullScale);
	// char useBits does not need to be reordered.
	// char kindBits does not need to be reordered.
	ReorderLong(&p->formula);
	ReorderLong(&p->depID);
	ReorderLong(&p->creationDate);
	// char wUnused does not need to be reordered.
	ReorderLong(&p->modDate);
	ReorderLong(&p->waveNoteH);
	// The wData field marks the start of the wave data which will be reordered separately.
}

static void
ReorderWaveHeader5(WaveHeader5* p)
{
	ReorderLong(&p->next);
	ReorderLong(&p->creationDate);
	ReorderLong(&p->modDate);
	ReorderLong(&p->npnts);
	ReorderShort(&p->type);
	ReorderShort(&p->dLock);
	// char whpad1 does not need to be reordered.
	ReorderShort(&p->whVersion);
	// char bname does not need to be reordered.
	ReorderLong(&p->whpad2);
	ReorderLong(&p->dFolder);
	ReorderBytes(&p->nDim, 4, 4);
	ReorderBytes(&p->sfA, 8, 4);
	ReorderBytes(&p->sfB, 8, 4);
	// char dataUnits does not need to be reordered.
	// char dimUnits does not need to be reordered.
	ReorderShort(&p->fsValid);
	ReorderShort(&p->whpad3);
	ReorderDouble(&p->topFullScale);
	ReorderDouble(&p->botFullScale);
	ReorderLong(&p->dataEUnits);
	ReorderBytes(&p->dimEUnits, 4, 4);
	ReorderBytes(&p->dimLabels, 4, 4);
	ReorderLong(&p->waveNoteH);
	ReorderBytes(&p->whUnused, 4, 16);
	ReorderShort(&p->aModified);
	ReorderShort(&p->wModified);
	ReorderShort(&p->swModified);
	// char useBits does not need to be reordered.
	// char kindBits does not need to be reordered.
	ReorderLong(&p->formula);
	ReorderLong(&p->depID);
	ReorderShort(&p->whpad4);
	ReorderShort(&p->srcFldr);
	ReorderLong(&p->fileName);
	ReorderLong(&p->sIndices);
	// The wData field marks the start of the wave data which will be reordered separately.
}

static int
Checksum(short *data, int needToReorderBytes, int oldcksum, int numbytes)
{
	unsigned short s;
	
	numbytes >>= 1;				// 2 bytes to a short -- ignore trailing odd byte.
	while(numbytes-- > 0) {
		s = *data++;
		if (needToReorderBytes)
			ReorderShort(&s);
		oldcksum += s;
	}
	return oldcksum&0xffff;
}

/*	NumBytesPerPoint(int type)
	
	Given a numeric wave type, returns the number of data bytes per point.
*/
static int
NumBytesPerPoint(int type)
{
	int numBytesPerPoint;
	
	// Consider the number type, not including the complex bit or the unsigned bit.
	switch(type & ~(NT_CMPLX | NT_UNSIGNED)) {
		case NT_I8:
			numBytesPerPoint = 1;		// char
			break;
		case NT_I16:
			numBytesPerPoint = 2;		// short
			break;
		case NT_I32:
			numBytesPerPoint = 4;		// long
			break;
		case NT_FP32:
			numBytesPerPoint = 4;		// float
			break;
		case NT_FP64:
			numBytesPerPoint = 8;		// double
			break;
		default:
			return 0;
			break;
	}

	if (type & NT_CMPLX)
		numBytesPerPoint *= 2;			// Complex wave - twice as many points.
	
	return numBytesPerPoint;
}

/*	LoadNumericWaveData(fr, type, npnts, waveDataSize, needToReorderBytes, pp)

	fr is a file reference.
	type is the Igor number type.
	npnts is the total number of elements in all dimensions.
	waveDataSize is the number of data bytes stored in the file.
	needToReorderBytes if the byte ordering of the file is not the byte ordering of the current platform.
	pp is a pointer to a pointer.
	
	If an error occurs, LoadWaveData returns a non-zero error code and sets
	*pp to NULL.
	
	If no error occurs, LoadWaveData returns 0 and sets *pp to a pointer allocated
	via malloc. This pointer must be freed by the calling routine.
*/
static int
LoadNumericWaveData(CP_FILE_REF fr, int type, long npnts, unsigned long waveDataSize, int needToReorderBytes, void**pp)
{
	int numBytesPerPoint;
	unsigned long numBytesToRead, numBytesToAllocate;
	unsigned long numBytesRead;
	void* p;
	int err;
	
	*pp = NULL;							// Assume that we can not allocate memory.

	numBytesPerPoint = NumBytesPerPoint(type);
	if (numBytesPerPoint <= 0) {
		printf("Invalid wave type (0x%x).\n", type);
		return -1;
	}
	numBytesToRead = npnts * numBytesPerPoint;

	numBytesToAllocate = numBytesToRead;
	if (numBytesToAllocate == 0)
		numBytesToAllocate = 8;			// This is just because malloc refuses to allocate a zero byte block.
	p = malloc(numBytesToAllocate);		// Allocate memory to store the wave data.
	if (p == NULL) {
		printf("Unable to allocate %ld bytes to store data.\n", numBytesToAllocate);
		return -1;
	}
	if (numBytesToRead > 0) {
		if (waveDataSize < numBytesToRead) {
			/*	If here, this should be a wave governed by a dependency formula
				for which no wave data was written to the file. Since we can't
				execute the dependency formula we have no way to recreate the wave's
				data. Therefore, we return 0 for all points in the wave.
			*/
			memset(p, 0, numBytesToRead);
		}
		else {
			if (err = CPReadFile(fr, numBytesToRead, p, &numBytesRead)) {
				free(p);
				printf("Error %d occurred while reading the wave data.\n", err);
				return err;
			}
			if (needToReorderBytes) {
				if (type != 0)				// Text wave data does not need to be reordered.
					ReorderBytes(p, numBytesPerPoint, numBytesToRead/numBytesPerPoint);
			}
		}
	}
	
	*pp = p;							// Return the pointer to the calling routine.
	return 0;
}


//new stuff, GH2015
static int
SetNumericWaveData(CP_FILE_REF fr, int type, long npnts, unsigned long waveDataSize, int needToReorderBytes, float rawData[])//changed because now rawData is a float array
{
	int numBytesPerPoint;
	unsigned long numBytesToRead, numBytesToAllocate;
	unsigned long numBytesRead;
	int err;



	numBytesPerPoint = NumBytesPerPoint(type);
	if (numBytesPerPoint <= 0) {
		printf("Invalid wave type (0x%x).\n", type);
		return -1;
	}
	numBytesToRead = npnts * numBytesPerPoint;

	numBytesToAllocate = numBytesToRead;
	if (numBytesToAllocate == 0)
		numBytesToAllocate = 8;			// This is just because malloc refuses to allocate a zero byte block.
// Allocate memory to store the wave data.

	if (numBytesToRead > 0) {
		if (waveDataSize < numBytesToRead) {
			/*	If here, this should be a wave governed by a dependency formula
			for which no wave data was written to the file. Since we can't
			execute the dependency formula we have no way to recreate the wave's
			data. Therefore, we return 0 for all points in the wave.
			*/
		}
		else {
			if (err = CPWriteFile(fr, waveDataSize, rawData, &numBytesRead)) {

				printf("Error %d occurred while reading the wave data.\n", err);
				return err;
			}

		}
	}


	return 0;
}



/*	ReadWaveAndBinHeader(CP_FILE_REF fr, WaveHeader5 *outWaveHeader5, BinHeader5 *outBinHeader5, int *pReorderBytes )

	This routine is used by several functions to read and return the information 
	about the IBW file.  Byte swapping will be applied to data that is being read
	on disimilar Indian processors.

	Return 0 to indicate success otherwise there was an error.
*/
int ReadWaveAndBinHeader(CP_FILE_REF fr, WaveHeader5 *outWaveHeader5, BinHeader5 *outBinHeader5, int *pReorderBytes )
{
	unsigned long startFilePos;
	short version;
	short check;
	int binHeaderSize, waveHeaderSize, checkSumSize;

	unsigned long numBytesRead;
	int needToReorderBytes;
	char buffer[512];
	int err;

	BinHeader5* b5 = NULL;
	WaveHeader5* w5 = NULL;

	// Set file position to head
	CPSetFilePosition(fr, 0, -1);

	if (err = CPGetFilePosition(fr, &startFilePos))
		return err;

	// Read the file version field.
	if (err = CPReadFile(fr, 2, &version, &numBytesRead)) {
		printf("Error %d occurred while reading the file version.\n", err);
		return err;
	}
	
	/*	Reorder version field bytes if necessary.
		If the low order byte of the version field of the BinHeader structure
		is zero then the file is from a platform that uses different byte-ordering
		and therefore all data will need to be reordered.
	*/
	needToReorderBytes = (version & 0xFF) == 0;
	if (needToReorderBytes)
		ReorderShort(&version);

	if( version < 5 )
	{
		printf("This does not appear to be a valid Igor binary wave file or its older than what we want to support. The version field = %d.\n", version);
		return -1;	
	}

	binHeaderSize = sizeof(BinHeader5);
	waveHeaderSize = sizeof(WaveHeader5);
	checkSumSize = binHeaderSize + waveHeaderSize - 4;	// Version 5 checksum does not include the wData field.

	// Load the BinHeader and the WaveHeader into memory.
	CPSetFilePosition(fr, startFilePos, -1);
	if (err = CPReadFile(fr, binHeaderSize+waveHeaderSize, buffer, &numBytesRead)) {
		printf("Error %d occurred while reading the file headers.\n", err);
		return err;
	}
	
	// Check the checksum.	
	check = Checksum((short*)buffer, needToReorderBytes, 0, checkSumSize);
	if (check != 0) {
		printf("Error in checksum - should be 0, is %d.\n", check);
		printf("This does not appear to be a valid Igor binary wave file.\n");
		return -1;
	}

	b5 = (BinHeader5*)buffer;
	w5 = (WaveHeader5*)(buffer+binHeaderSize);

	// Do byte reordering if the file is from another platform.	
	if (needToReorderBytes) 
	{
		ReorderBinHeader5(b5);
		ReorderWaveHeader5(w5);
	}

	// Return information to the calling routine.
	if( outBinHeader5 )
		memcpy( outBinHeader5, b5, sizeof(BinHeader5) );
	
	if( outWaveHeader5 )
		memcpy( outWaveHeader5, w5, sizeof(WaveHeader5) );

	if( pReorderBytes ) 
		*pReorderBytes = needToReorderBytes;
	
	return 0;
}

/*	ReadARImageType(CP_FILE_REF fr)

	Return the image type of the IBW file.  If it is an Asylum branded file,
	the last 6 bytes of the file will contain either "MFP3D" or "Force" to 
	indicate an image or forceplot, respectively.

	See cARType_* constants in the header file for a description of the types
	that can be returned.
	
*/

int ReadARImageType(CP_FILE_REF fr)
{
	char szType[6]={0};
	unsigned long fileSize = 0;
	
	CPNumberOfBytesInFile(fr, &fileSize);
	CPSetFilePosition(fr, fileSize - 5, -1);	// Walk 6 bytes back from the end of the file

	CPReadFile(fr, 5, szType, NULL);

	if( !strcmp(szType,"MFP3D") )
		return cARType_Image;
	else if( !strcmp(szType,"Force") )
		return cARType_ForcePlot;

	return cARType_Unknown;
}

/*	ReadDimLabels(CP_FILE_REF fr, tDimLabelInfo ** )

	This will return an information structure that describes the dimension
	labels in the wave.  The tDimLabelInfo block is created by the function
	and must be destroyed by a call to DestroyDimLabels(...) after use.

	return 0 means success otherwise there was a failure creating the block
	and the pointer passed back by reference will not be valid.
*/

int ReadDimLabels(CP_FILE_REF fr, tDimLabelInfo **pOutDimInfo )
{
	WaveHeader5 w5;
	BinHeader5 b5;
	int i, j;
	int noteStart, noteLength;
	int curLabelsPerDim = 0;
	int needToReorderBytes = 0;
	int totalIndividualDimLabels = 0;
	int blockSize = 0;
	tDimLabelInfo *pDimInfo = NULL;
	tDimLabel *pLabelLocation = NULL;

	if( ReadWaveAndBinHeader(fr, &w5, &b5, &needToReorderBytes ) )
		return -1;

	if( pOutDimInfo == NULL )
		return 1; // can't do anything if there is no pointer

	noteStart = b5.wfmSize + sizeof(BinHeader5);
	noteLength = b5.noteSize;

	// Set file position to just after the note.  That's were the dimlabels start
	CPSetFilePosition(fr, noteStart + noteLength, -1);

	// Now allocate the structure for holding the dim labels.
	*pOutDimInfo = (tDimLabelInfo *)malloc(sizeof(tDimLabelInfo));
	pDimInfo = *pOutDimInfo;

	memset((char *)pDimInfo, 0, sizeof(tDimLabelInfo));	// zero out the memory

	for( i=0; i<MAXDIMS; i++ )
	{
		if( b5.dimLabelsSize[i] )
		{
			// Since there are always 32 chars per label, we can just divide the total by 32 to get
			// the number of labels.  When there are labels for a dimension, there is always
			// a reserved name of a global label for that dimension in the first spot.  Since
			// we do not count global dimension labels in the total, we subract 1 from the total.
			curLabelsPerDim = b5.dimLabelsSize[i]/32-1;
			pDimInfo->nLabels[i] = curLabelsPerDim;

			// Read the whole dimension label
			CPReadFile(fr, sizeof(tDimLabel), (void *)pDimInfo->dimLabelWhole[i].label, NULL);

			// If there are labels for the specified dimension, allocate the block and load them.
			if( curLabelsPerDim )
			{
				pDimInfo->dimLabels[i] = (tDimLabel*)malloc(sizeof(tDimLabel)*curLabelsPerDim);

				for( j=0; j<curLabelsPerDim; j++)
					CPReadFile(fr, sizeof(tDimLabel), (void *)&pDimInfo->dimLabels[i][j], NULL);
			}
		}
	}

	return 0;
}

/*	DestroyDimLabels( tDimLabelInfo * )

	This call is required to properly destroy a tDimLabelInfo memory block
	that was created by a call to ReadDimLabels.  Since the memory is dynamic
	we need to destroy it selectively

*/

void DestroyDimLabels( tDimLabelInfo *pDimLabels )
{
	int i;

	if( pDimLabels )
	{
		// Free the individual label chunks first
		for(i=0; i<MAXDIMS; i++)
			free( (char *)pDimLabels->dimLabels[i] );

		// Now free the struct as a whole
		free( (char *)pDimLabels );
	}
}

/*	ReadWaveNote(CP_FILE_REF fr, char **pNoteBuffer, int *nNoteBufferSize)

	Allocate a buffer and return its pointer and size into the passed
	in parameters.  It is up to the caller to destroy the memory
	holding the note.

	The separator between name and value is ":" while each note entry
	is delimeted with '\x0d'.  There is also a ' ' preceding the
	value that can be trimmed.

	0 is returned for success otherwise there was an error.
*/
int ReadWaveNote(CP_FILE_REF fr, char **pNoteBuffer, int *nNoteSize)
{
	WaveHeader5 w5;
	BinHeader5 b5;
	int needToReorderBytes = 0;
	int	noteStart = 0;
	int	dependLength = 0;
	int	noteLength = 0;
	unsigned long numBytesRead = 0;

	ReadWaveAndBinHeader(fr, &w5, &b5, &needToReorderBytes );

	noteStart = b5.wfmSize;
	dependLength = b5.formulaSize;
	noteLength = b5.noteSize;
	noteStart += dependLength + sizeof(BinHeader5);

	// Just requesting the size.
	if( pNoteBuffer == NULL )
	{
		if( nNoteSize )
			*nNoteSize = noteLength + 2; // Adding 2 bytes of NULL for string padding.
		return 0;
	}
	
	CPSetFilePosition(fr, noteStart, -1);

	*pNoteBuffer = (char *)malloc(noteLength + 2);

	if( *pNoteBuffer == NULL )
		return -1; // problem with allocation

	*(*pNoteBuffer+noteLength) = 0x0d; // pad with 0x0d for proper parsing used in GetNameValue()
	*(*pNoteBuffer+noteLength+1) = 0x0d;

	CPReadFile(fr, noteLength, *pNoteBuffer, &numBytesRead);

	*nNoteSize = numBytesRead;

	return 0;
}

/*	ReadWave(fr, typePtr, npntsPtr, waveDataPtrPtr)

	Reads the wave file and prints some information about it.
	
	Returns to the calling routine the wave's type, number of points, and the
	wave data. The calling routine must free *waveDataPtrPtr if it is
	not null.
	
	Returns 0 or an error code.
	
	This routine is written such that it could be used to read waves
	from an Igor packed experiment file as well as from a standalone
	Igor binary wave file. In order to achieve this, we must not assume
	that the wave is at the start of the file. We do assume that, on entry
	to this routine, the file position is at the start of the wave.
*/
int ReadWave(CP_FILE_REF fr, int* typePtr, long* npntsPtr, void** waveDataPtrPtr, struct WaveHeader5 *pWaveDesc /*= NULL*/)
{
	unsigned long startFilePos;

	unsigned long waveDataSize;
	int needToReorderBytes;
	long wfmSize;
	int err;
	struct WaveHeader5 w5;
	struct BinHeader5 b5;

	WaveHeader5 *pWaveHeader5Instance = NULL;
	
	*waveDataPtrPtr = NULL;
	*typePtr = 0;
	*npntsPtr = 0;
	
	ReadWaveAndBinHeader(fr, &w5, &b5, &needToReorderBytes );

	// Return information to the calling routine.
	if( typePtr )
		*typePtr = w5.type;
	if( npntsPtr )
		*npntsPtr = w5.npnts;

	wfmSize = b5.wfmSize;

	// Determine the number of bytes of wave data in the file.
	waveDataSize = wfmSize - offsetof(WaveHeader5, wData);

	// Set file position to head
	CPSetFilePosition(fr, 0, -1);
	if (err = CPGetFilePosition(fr, &startFilePos))
		return err;

	// Position the file pointer to the start of the wData field.
	CPSetFilePosition(fr, startFilePos+sizeof(BinHeader5)+sizeof(WaveHeader5)-4, -1);		// 4 = size of wData field in WaveHeader2 structure.
	
	if (w5.type == 0) {
		// For simplicity, we don't load text wave data in this example program.
		printf("This is a text wave.\n");
		return 0;
	}

	// Load the data and allocates memory to store it.
	if (err = LoadNumericWaveData(fr, w5.type, w5.npnts, waveDataSize, needToReorderBytes, waveDataPtrPtr))
		return err;
	
	if( pWaveDesc )
	{
		memcpy(pWaveDesc, &w5, sizeof(WaveHeader5));
	}

	return 0;
}


//GH, 2015
int ReadWaveArgo(CP_FILE_REF fr, int* typePtr, long* npntsPtr, void** waveDataPtrPtr, struct WaveHeader5 *pWaveDesc, float rawData[] /*= NULL*/)
{
	unsigned long startFilePos;

	unsigned long waveDataSize;
	int needToReorderBytes;
	long wfmSize;
	int err;
	struct WaveHeader5 w5;
	struct BinHeader5 b5;

	WaveHeader5 *pWaveHeader5Instance = NULL;

	*waveDataPtrPtr = NULL;
	*typePtr = 0;
	*npntsPtr = 0;

	ReadWaveAndBinHeader(fr, &w5, &b5, &needToReorderBytes);

	// Return information to the calling routine.
	if (typePtr)
		*typePtr = w5.type;
	if (npntsPtr)
		*npntsPtr = w5.npnts;

	wfmSize = b5.wfmSize;

	// Determine the number of bytes of wave data in the file.
	waveDataSize = wfmSize - offsetof(WaveHeader5, wData);

	// Set file position to head
	CPSetFilePosition(fr, 0, -1);
	if (err = CPGetFilePosition(fr, &startFilePos))
		return err;

	// Position the file pointer to the start of the wData field.
	CPSetFilePosition(fr, startFilePos + sizeof(BinHeader5) + sizeof(WaveHeader5) - 4, -1);		// 4 = size of wData field in WaveHeader2 structure.

	if (w5.type == 0) {
		// For simplicity, we don't load text wave data in this example program.
		printf("This is a text wave.\n");
		return 0;
	}

	// Load the data and allocates memory to store it.
	if (err = SetNumericWaveData(fr, w5.type, w5.npnts, waveDataSize, needToReorderBytes, rawData))
		return err;

	if (pWaveDesc)
	{
		memcpy(pWaveDesc, &w5, sizeof(WaveHeader5));
	}

	return 0;
}

/*	GetUnitsFromTitle(const char *pChannelTitle)

	The units for the individual channels are not stored in the
	note or IBW headers.  To retrieve the units, the title
	of the channel may be parsed and cross referenced among the
	titles as shown within this function.

	If the keyword is found in the channel title, the 
	respective units will be returned otherwise "V" will
	be returned as a default.
*/
const char * GetUnitsFromTitle(const char *pChannelTitle)
{
	char szChannelBuf[32]={0};
	strcpy(szChannelBuf, pChannelTitle);
	_strlwr(szChannelBuf);

	if( strstr(szChannelBuf, "height") != NULL )
		return "m";

	if( strstr(szChannelBuf, "zsensor") != NULL )
		return "m";

	if( strstr(szChannelBuf, "amplitude") != NULL )
		return "m";

	if( strstr(szChannelBuf, "phase") != NULL )
		return "\x0b0"; // degrees symbol

	if( strstr(szChannelBuf, "potential") != NULL )
		return "V";

	if( strstr(szChannelBuf, "current") != NULL )
		return "A";

	if( strstr(szChannelBuf, "frequency") != NULL )
		return "Hz";

	if( strstr(szChannelBuf, "capacitance") != NULL )
		return "f";

	if( strstr(szChannelBuf, "count") != NULL )
		return "";

	return "V"; // default for now
}


int DoReadTest(const char* filePath);

int DoReadTest(const char* filePath)
{
	CP_FILE_REF fr;
	int type;
	long npnts;
	void* waveDataPtr;
	int err;
	
	if (err = CPOpenFile(filePath, 0, &fr)) {
		printf("Error %d occurred while opening the file.\n", err);
		return err;
	}
	
	// Here you would do something with the data.


	err = ReadWave(fr, &type, &npnts, &waveDataPtr, NULL);
	
	CPCloseFile(fr);
	
	if (waveDataPtr != NULL)
		free(waveDataPtr);
	
	printf("End of read test.\n");
	
	return err;
}

extern "C" __declspec( dllexport ) int ReleaseMemory( int* pArray ) {
	delete[] pArray;
	return 0;
}

extern "C" __declspec( dllexport )
int GetNumberOfLayers( const char * filePath ) {
	CP_FILE_REF fr;
	int err;
	int arImageType;
	tDimLabelInfo *pDimLabelInfo = NULL;

	if ( err = CPOpenFile( filePath, 0, &fr ) ) {
		printf( "Error %d occurred while opening the file.\n", err );
		return err;
	}

	// Determine the Asylum type of image (i.e. Image or forceplot)
	arImageType = ReadARImageType( fr );

	switch ( arImageType ) {
		case cARType_Image:
			break;

		case cARType_ForcePlot:
			printf( "Image type:	Force plot\n\n" );
			printf( "*** Not processing force plot files in this example - exiting.\n" );
			goto end;

		default:
			printf( "*** Not an Asylum type IBW file.\n\n" );
			goto end;
	}

	//Let's dump stats on the file and get the channels

	if ( ReadDimLabels( fr, &pDimLabelInfo ) == 0 ) {
		err = pDimLabelInfo->nLabels[ 2 ];
		DestroyDimLabels( pDimLabelInfo );
	}

end:

	// Done reading the file and gathering info
	CPCloseFile( fr );
	return err;
}

extern "C" __declspec( dllexport )
char * GetLayerName( const char * filePath, int index ) {
	CP_FILE_REF fr;
	int err;
	int arImageType;
	char * name = new char[ 32 ];
	tDimLabelInfo *pDimLabelInfo = NULL;

	if ( err = CPOpenFile( filePath, 0, &fr ) ) {
		printf( "Error %d occurred while opening the file.\n", err );
		return "Error";
	}

	// Determine the Asylum type of image (i.e. Image or forceplot)
	arImageType = ReadARImageType( fr );

	switch ( arImageType ) {
		case cARType_Image:
			break;

		case cARType_ForcePlot:
			printf( "Image type:	Force plot\n\n" );
			printf( "*** Not processing force plot files in this example - exiting.\n" );
			goto end;

		default:
			printf( "*** Not an Asylum type IBW file.\n\n" );
			goto end;
	}

	//Let's dump stats on the file and get the channels

	if ( ReadDimLabels( fr, &pDimLabelInfo ) == 0 ) {
		name = strncpy( name, pDimLabelInfo->dimLabels[ 2 ][ index ].label, 32 );
		DestroyDimLabels( pDimLabelInfo );
	}

end:

	// Done reading the file and gathering info
	CPCloseFile( fr );
	return name;
}

extern "C" __declspec( dllexport )
int * GetDimensions( const char * filePath ) {
	CP_FILE_REF fr;
	int err;
	int arImageType;
	int type;
	long npnts;
	void* waveDataPtr;
	WaveHeader5 w5;
	int* dimensions = new int[ 3 ];

	if ( err = CPOpenFile( filePath, 0, &fr ) ) {
		printf( "Error %d occurred while opening the file.\n", err );
		return NULL;
	}

	// Determine the Asylum type of image (i.e. Image or forceplot)
	arImageType = ReadARImageType( fr );

	switch ( arImageType ) {
		case cARType_Image:
			break;

		case cARType_ForcePlot:
			printf( "Image type:	Force plot\n\n" );
			printf( "*** Not processing force plot files in this example - exiting.\n" );
			goto end;

		default:
			printf( "*** Not an Asylum type IBW file.\n\n" );
			goto end;
	}

	// Read the data out of the IBW.
	if ( ( err = ReadWave( fr, &type, &npnts, &waveDataPtr, &w5 ) ) != 0 ) {
		printf( "*** ReadWave(...) error, exiting.\n\n" );
		goto end;
	}

	dimensions[ 0 ] = w5.nDim[ 0 ];
	dimensions[ 1 ] = w5.nDim[ 1 ];
	dimensions[ 2 ] = w5.nDim[ 2 ];
end:

	// Done reading the file and gathering info

	if ( waveDataPtr != NULL )
		free( waveDataPtr );
	CPCloseFile( fr );
	return dimensions;
}

extern "C" __declspec( dllexport )
float * GetRawData( const char * filePath, unsigned int size ) {
	CP_FILE_REF fr;
	int type;
	long npnts;
	void* waveDataPtr;
	int err;
	WaveHeader5 w5;

	if ( err = CPOpenFile( filePath, 0, &fr ) ) {
		printf( "Error %d occurred while opening the file.\n", err );
		return NULL;
	}

	// Read the data out of the IBW.
	if ( ( err = ReadWave( fr, &type, &npnts, &waveDataPtr, &w5 ) ) != 0 ) {
		printf( "*** ReadWave(...) error, exiting.\n\n" );
		goto end;
	}

	float * rawData = new float[ size ];
	memcpy( rawData, waveDataPtr, size * sizeof(float) );
end:

	// Done reading the file and gathering info
	CPCloseFile( fr );

	if ( waveDataPtr != NULL )
		free( waveDataPtr );

	return rawData;
}


//this is new stuff, Setting the raw data for the new IBW file after Argo, GH 2015
//here I'd need to ReadWave, get the data, then put in the new data and write a new file using WriteWave 
//instead of being a float pointer, it now references a float array
extern "C" __declspec( dllexport )
int SetRawData( const char * filePath, unsigned int size, float rawData[] ) {
	CP_FILE_REF fr;
	int type;
	long npnts;
	void* waveDataPtr;
	int err;
	WaveHeader5 w5;
	//was 0, but we need to write to the file GH 2016
	if ( err = CPOpenFile( filePath, 1, &fr ) ) {
		printf( "Error %d occurred while opening the file.\n", err );
		return NULL;
	}

	// Set the data into the IBW.
	if ( ( err = ReadWaveArgo( fr, &type, &npnts, &waveDataPtr, &w5, rawData ) ) != 0 ) {
		printf( "*** ReadWave(...) error, exiting.\n\n" );
		goto end;
	}
	
end:

	// Done reading the file and gathering info
	CPCloseFile( fr );

	if ( waveDataPtr != NULL )
		free( waveDataPtr );

	return 0;
}


extern "C" __declspec( dllexport ) 
int DoAsylumReadTest( const char * filePath )
{
	CP_FILE_REF fr;
	int type;
	long npnts;
	void* waveDataPtr;
	int err;
	int i;
	const char *pUnits;
	int arImageType;
	tDimLabelInfo *pDimLabelInfo = NULL;
	WaveHeader5 w5;
	char dataUnits[2][MAX_UNIT_CHARS+1];

	int xDim, yDim, zLayers;
	int xyPoints;
	double xStart, xScale, yStart, yScale;

	float *pRawData=NULL;
	
	if (err = CPOpenFile(filePath, 0, &fr)) {
		printf("Error %d occurred while opening the file.\n", err);
		return err;
	}
	
	printf("Asylum Research File Test:\n");

	// Determine the Asylum type of image (i.e. Image or forceplot)
	arImageType = ReadARImageType(fr);

	switch(arImageType)
	{
	case cARType_Image:
		printf("Image type:	Image\n\n");
		break;

	case cARType_ForcePlot:
		printf("Image type:	Force plot\n\n");
		printf("*** Not processing force plot files in this example - exiting.\n");
		goto end;

	default:
		printf("*** Not an Asylum type IBW file.\n\n");
		goto end;
	}

	// Read the data out of the IBW.
	if( (err = ReadWave(fr, &type, &npnts, &waveDataPtr, &w5)) != 0 )
	{
		printf("*** ReadWave(...) error, exiting.\n\n");
		goto end;
	}
	
	
	//Let's dump stats on the file and get the channels

	if( ReadDimLabels(fr, &pDimLabelInfo ) == 0 )
	{
		printf( "Channel listing by name and units - %d channels\n", pDimLabelInfo->nLabels[2]);

		// The layer dimension labels hold the channel title (0 = row, 1 = column, 2 = layer)
		for( i=0; i<pDimLabelInfo->nLabels[2]; i++ )
		{
			pUnits = GetUnitsFromTitle(pDimLabelInfo->dimLabels[2][i].label);
			printf( "\t%d - %s ( %s )\n",i+1, pDimLabelInfo->dimLabels[2][i].label, pUnits);
		}
		// Clean up the allocated block of memory describing the labels
		DestroyDimLabels( pDimLabelInfo );
	}


	// Dump the x/y dimensions of the file.
	xDim = w5.nDim[0];
	yDim = w5.nDim[1];
	
	xStart = w5.sfB[0];
	yStart = w5.sfB[1];

	xScale = w5.sfA[0]*(xDim-1);
	yScale = w5.sfA[1]*(yDim-1);

	strcpy(dataUnits[0], w5.dimUnits[0] );
	strcpy(dataUnits[1], w5.dimUnits[1] );

	printf("\n\n");
	printf("Image dimension information:\n");
	printf("\tX/Y points: %d, %d\n", xDim, yDim);
	printf("\tX/Y units: %s, %s\n", dataUnits[0], dataUnits[1] );
	printf("\tX Start: %g %s\n", xStart, dataUnits[0]);
	printf("\tY Start: %g %s\n", yStart, dataUnits[1]);
	printf("\tX Scale: %g %s\n", xScale, dataUnits[0]);
	printf("\tY Scale: %g %s\n", yScale, dataUnits[1]);

	// Show the offset of where the data starts for each layer and several points of data.
	printf("\n\n");
	printf("Display offset into wave data for each channel plus a few points of data.\n");

	zLayers = w5.nDim[2];
	xyPoints = xDim*yDim;

	printf("Channel\tOffset\tData\n");
	for( i=0; i<zLayers; i++ )
	{
		pRawData = (float*)waveDataPtr + xDim*yDim*i;
		printf("%d - %7d - %g, %g, %g, %g\n", i+1, xDim*yDim*i, 
			*(pRawData+0), *(pRawData+1), *(pRawData+2), *(pRawData+3) );
	}
	

	// NOTE: Last saved range and offset values
	// To find the preferred range and offset of each channel, 
	// the note will need to be read and look for "Display range #"
	// and "Display offset #" where # is the channel starting at 0.
	// The value associated with the name/value pair can be used as
	// a starting point for the range and offset when displaying the data.

	// An example of this will be added at a later time but for now you
	// can read the note block with 
	// ReadWaveNote(CP_FILE_REF fr, char **pNoteBuffer, int *nNoteBufferSize)
	// and parse the name/value pairs into your own structure and search.


	printf("\nNOTE:\tTo find the preferred range and offset of each channel,\n \
	the note will need to be read and look for \"Display range #\" \n \
	and \"Display offset #\" where # is the channel starting at 0. \n \
	The value associated with the name/value pair can be used as \n \
	a starting point for the range and offset when displaying the data.\n\n \
	An example of this will be added at a later time but for now you \n \
	can read the note block with \n \
	ReadWaveNote(CP_FILE_REF fr, char **pNoteBuffer, int *nNoteBufferSize) \n \
	and parse the name/value pairs into your own structure and search.\n\n" );

end:

	// Done reading the file and gathering info
	CPCloseFile(fr);
	
	if (waveDataPtr != NULL)
		free(waveDataPtr);
	
	printf("End of read test.\n");
	
	return err;
}

