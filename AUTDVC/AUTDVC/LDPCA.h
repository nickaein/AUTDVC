
void LDPCAencodeBits(const char *InitByLadderFile, double *source, double *accumulatedSyndrome);

void LDPCAdecodeBits(const char* InitByLadderFile, double *LLR_intrinsic, double *accumulatedSyndrome, double *source,
                double *decoded, double *rate, double *numErrors);