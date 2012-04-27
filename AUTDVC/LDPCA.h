
void LDPCAencodeBits(char *InitByLadderFile, double *source, double *accumulatedSyndrome);

void LDPCAdecodeBits(char* InitByLadderFile, double *LLR_intrinsic, double *accumulatedSyndrome, double *source,
                double *decoded, double *rate, double *numErrors);