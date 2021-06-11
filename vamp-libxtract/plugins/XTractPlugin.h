/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugins using Jamie Bullock's
    libxtract audio feature extraction library.

    Centre for Digital Music, Queen Mary, University of London.
    This file copyright 2006 Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef _XTRACT_PLUGIN_H_
#define _XTRACT_PLUGIN_H_

#include <vamp-sdk/Plugin.h>
#include "../LibXtract/include/xtract/libxtract.h"

class XTractPlugin : public Vamp::Plugin
{
public:
    XTractPlugin(unsigned int xtFeature, float inputSampleRate);
    virtual ~XTractPlugin();

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    InputDomain getInputDomain() const;

    std::string getIdentifier() const;
    std::string getName() const;
    std::string getDescription() const;
    std::string getMaker() const;
    int getPluginVersion() const;
    std::string getCopyright() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(std::string) const;
    void setParameter(std::string, float);

    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    size_t getPreferredStepSize() const;
    size_t getPreferredBlockSize() const;

    OutputList getOutputDescriptors() const;

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    bool needPeakThreshold() const;
    bool needHarmonicThreshold() const;
    bool needRolloffThreshold() const;

    mutable OutputList m_outputDescriptors;
    void setupOutputDescriptors() const;

    bool processSPF0(const double *data);

    const unsigned int m_xtFeature;
    size_t m_channels;
    size_t m_stepSize;
    size_t m_blockSize;

    double *m_resultBuffer;

    double m_peakThreshold;
    double m_rolloffThreshold;
    double m_harmonicThreshold;

    double m_minFreq;
    double m_maxFreq;

    int m_coeffCount;
    int m_highestCoef;
    int m_lowestCoef;
    double **m_mfccFilters;
    int m_mfccStyle;

    int m_spectrumType;
    int m_dc;
    int m_normalise;

    int *m_barkBandLimits;

    static xtract_function_descriptor_t *m_xtDescriptors;
    static int m_xtDescRefCount;
    xtract_function_descriptor_t *xtDescriptor() {
        return &m_xtDescriptors[m_xtFeature];
    }
    const xtract_function_descriptor_t *xtDescriptor() const {
        return &m_xtDescriptors[m_xtFeature];
    }

    size_t m_outputBinCount;
    bool m_initialised;
    static bool m_anyInitialised;
};


#endif
