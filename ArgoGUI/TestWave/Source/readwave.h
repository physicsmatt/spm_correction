#ifndef READWAVE_H
#define READWAVE_H

#include "IgorBin.h"
#include "CrossPlatformFileIO.h"

typedef struct {
	char label[32];
} tDimLabel;

typedef struct {
	int	nLabels[MAXDIMS];
	tDimLabel dimLabelWhole[MAXDIMS];
	tDimLabel *dimLabels[MAXDIMS];	// This block is dynamically allocated to hold the correct
									// amount of labels based on the number per dimension as 
									// specified by nLabels.

} tDimLabelInfo;



#define cARType_Unknown		0
#define cARType_Image		1
#define cARType_ForcePlot	2

int ReadWaveAndBinHeader(CP_FILE_REF fr, WaveHeader5 *outWaveHeader5, BinHeader5 *outBinHeader5, int *pReorderBytes );

int ReadARImageType(CP_FILE_REF fr);
int ReadDimLabels(CP_FILE_REF fr, tDimLabelInfo ** );
void DestroyDimLabels( tDimLabelInfo * );
int ReadWaveNote(CP_FILE_REF fr, char **pNoteBuffer, int *nNoteBufferSize);
int ReadWave(CP_FILE_REF fr, int* typePtr, long* npntsPtr, void** waveDataPtrPtr, struct WaveHeader5 *pWaveDesc);
const char * GetUnitsFromTitle(const char *pChannelTitle);

#endif // READWAVE_H
