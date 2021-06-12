/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/*
    Vamp feature extraction plugin using the OnsetsDS onset detector.
    This file copyright (c) 2008 Chris Cannam.

    OnsetsDS - real time musical onset detection library.
    Copyright (c) 2007 Dan Stowell. All rights reserved.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _ONSETSDS_PLUGIN_H_
#define _ONSETSDS_PLUGIN_H_

#include <vamp-sdk/Plugin.h>

#include "onsetsds/onsetsds.h"

class OnsetsDSPlugin : public Vamp::Plugin
{
public:
    OnsetsDSPlugin(float inputSampleRate);
    virtual ~OnsetsDSPlugin();

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    InputDomain getInputDomain() const { return FrequencyDomain; }

    std::string getIdentifier() const;
    std::string getName() const;
    std::string getDescription() const;
    std::string getMaker() const;
    int getPluginVersion() const;
    std::string getCopyright() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(std::string) const;
    void setParameter(std::string, float);

    size_t getPreferredStepSize() const;
    size_t getPreferredBlockSize() const;

    OutputList getOutputDescriptors() const;

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    OnsetsDS *m_ods;
    float *m_odsdata;
    onsetsds_odf_types m_dfType;
    float m_threshold;
    size_t m_medspan;
    size_t m_stepSize;
    size_t m_fftSize;
};


#endif
