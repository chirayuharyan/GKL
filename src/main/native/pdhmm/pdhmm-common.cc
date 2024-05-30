#include "pdhmm-common.h"

double g_matchToMatchLog10[(((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1)] __attribute__((aligned(64)));
double g_matchToMatchProb[(((MAX_QUAL + 1) * (MAX_QUAL + 2)) >> 1)] __attribute__((aligned(64)));
double g_qualToErrorProbCache[(MAX_QUAL + 1)] __attribute__((aligned(64)));
double g_qualToProbLog10Cache[(MAX_QUAL + 1)] __attribute__((aligned(64)));

inline double qualToErrorProb(double qual, int32_t &status)
{
    if (qual < 0.0)
    {
        status = PDHMM_INPUT_DATA_ERROR;
        DBG("deletion quality cannot be less than 0 \n");
    }
    return pow(10.0, qual / -10.0);
}

int32_t initializeCache()
{
    int32_t status = PDHMM_SUCCESS;
    /* Step 1  */
    JacobianLogTable::initCache();

    for (int32_t i = 0, offset = 0; i <= MAX_QUAL; offset += ++i)
    {
        for (int32_t j = 0; j <= i; j++)
        {
            double log10Sum = approximateLog10SumLog10(-0.1 * i, -0.1 * j);
            g_matchToMatchLog10[offset + j] =
                log1p(-std::min(1.0, pow(10, log10Sum))) * INV_LN10;
            g_matchToMatchProb[offset + j] = pow(10, g_matchToMatchLog10[offset + j]);
        }
    }

    /* Step 3 */

    for (int32_t i = 0; i <= MAX_QUAL; i++)
    {
        g_qualToErrorProbCache[i] = qualToErrorProb((double)i, status);
        if (status != PDHMM_SUCCESS)
        {
            return status;
        }
        g_qualToProbLog10Cache[i] = log10(1.0 - g_qualToErrorProbCache[i]);
    }
    return status;
}

int32_t init(double *&matchToMatchLog10, double *&matchToMatchProb,
             double *&qualToErrorProbCache,
             double *&qualToProbLog10Cache)
{
    matchToMatchLog10 = g_matchToMatchLog10;
    matchToMatchProb = g_matchToMatchProb;
    qualToErrorProbCache = g_qualToErrorProbCache;
    qualToProbLog10Cache = g_qualToProbLog10Cache;
    return PDHMM_SUCCESS;
}
