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

#include "onsetsdsplugin.h"
#include <vamp-sdk/PluginAdapter.h>

using std::string;
using std::vector;
using std::cerr;
using std::endl;


OnsetsDSPlugin::OnsetsDSPlugin(float inputSampleRate) :
    Vamp::Plugin(inputSampleRate),
    m_ods(0),
    m_odsdata(0),
    m_dfType(ODS_ODF_RCOMPLEX),
    m_threshold(0.5),
    m_medspan(11),
    m_stepSize(256),
    m_fftSize(512)
{
}

OnsetsDSPlugin::~OnsetsDSPlugin()
{
    delete[] m_odsdata;
    delete m_ods;
}

string
OnsetsDSPlugin::getIdentifier() const
{
    return "onsetsds";
}

string
OnsetsDSPlugin::getName() const
{
    return "OnsetsDS Onset Detector";
}

string
OnsetsDSPlugin::getDescription() const
{
    return "Detect note onsets";
}

string
OnsetsDSPlugin::getMaker() const
{
    return "Dan Stowell";
}

int
OnsetsDSPlugin::getPluginVersion() const
{
    return 1;
}

string
OnsetsDSPlugin::getCopyright() const
{
    return "Copyright (c) 2007-2008 Dan Stowell";
}

OnsetsDSPlugin::ParameterList
OnsetsDSPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor desc;
    desc.identifier = "dftype";
    desc.name = "Onset detection function";
    desc.description = "Method used to calculate the onset detection function";
    desc.minValue = 0;
    desc.maxValue = 6;
    desc.defaultValue = 3;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.push_back("Power");
    desc.valueNames.push_back("Sum of magnitudes");
    desc.valueNames.push_back("Complex-domain deviation");
    desc.valueNames.push_back("Rectified complex-domain deviation");
    desc.valueNames.push_back("Phase deviation");
    desc.valueNames.push_back("Weighted phase deviation");
    desc.valueNames.push_back("Modified Kullback-Liebler deviation");
    list.push_back(desc);

    desc.identifier = "threshold";
    desc.name = "Detection threshold";
    desc.description = "Onsets trigger when the function beats this value";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 0.5;
    desc.isQuantized = false;
    desc.valueNames.clear();
    list.push_back(desc);

    desc.identifier = "medspan";
    desc.name = "Median frame span";
    desc.description = "Number of past frames used in median calculation";
    desc.minValue = 5;
    desc.maxValue = 21;
    desc.defaultValue = 11;
    desc.isQuantized = true;
    desc.quantizeStep = 2;
    desc.valueNames.clear();
    list.push_back(desc);

    return list;
}

float
OnsetsDSPlugin::getParameter(std::string name) const
{
    if (name == "dftype") {
        switch (m_dfType) {
        case ODS_ODF_POWER:    return 0;
        case ODS_ODF_MAGSUM:   return 1;
        case ODS_ODF_COMPLEX:  return 2;
        case ODS_ODF_RCOMPLEX: return 3;
        case ODS_ODF_PHASE:    return 4;
        case ODS_ODF_WPHASE:   return 5;
        case ODS_ODF_MKL:      return 6;
        }
    } else if (name == "threshold") {
        return m_threshold;
    } else if (name == "medspan") {
        return m_medspan;
    }
    return 0.0;
}

void
OnsetsDSPlugin::setParameter(std::string name, float value)
{
    if (name == "dftype") {
        onsetsds_odf_types dfType = m_dfType;
        switch (lrintf(value)) {
        case 0: dfType = ODS_ODF_POWER; break;
        case 1: dfType = ODS_ODF_MAGSUM; break;
        case 2: dfType = ODS_ODF_COMPLEX; break;
        case 3: dfType = ODS_ODF_RCOMPLEX; break;
        case 4: dfType = ODS_ODF_PHASE; break;
        case 5: dfType = ODS_ODF_WPHASE; break;
        case 6: dfType = ODS_ODF_MKL; break;
        }
        if (dfType == m_dfType) return;
        m_dfType = dfType;
    } else if (name == "threshold") {
        m_threshold = value;
    } else if (name == "medspan") {
        m_medspan = lrintf(value);
    }
}

bool
OnsetsDSPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) {
        std::cerr << "OnsetsDSPlugin::initialise: Unsupported channel count: "
                  << channels << std::endl;
        return false;
    }

    if (stepSize != getPreferredStepSize()) {
        std::cerr << "WARNING: OnsetsDSPlugin::initialise: Using unusual step size: "
                  << stepSize << " (wanted " << (getPreferredStepSize()) << ")" << std::endl;
    }

    if (blockSize != getPreferredBlockSize()) {
        std::cerr << "WARNING: OnsetsDSPlugin::initialise: Using unusual block size: "
                  << blockSize << " (wanted " << (getPreferredBlockSize()) << ")" << std::endl;
    }

    m_stepSize = stepSize;
    m_fftSize = blockSize;

    delete[] m_odsdata;
    delete m_ods;

    m_odsdata = new float[onsetsds_memneeded(m_dfType, m_fftSize, m_medspan)];
    m_ods = new OnsetsDS;
    memset(m_ods, 0, sizeof(OnsetsDS));
    onsetsds_init(m_ods, m_odsdata, ODS_FFT_FFTW3_R2C, m_dfType, m_fftSize,
                  m_medspan, m_inputSampleRate);
	m_ods->thresh = m_threshold;
	
    return true;
}

void
OnsetsDSPlugin::reset()
{
    if (!m_ods) {
        std::cerr << "ERROR: OnsetsDSPlugin::reset: Plugin has not been initialised" << std::endl;
        return;
    }
    onsetsds_init(m_ods, m_odsdata, ODS_FFT_FFTW3_R2C, m_dfType, m_fftSize,
                  m_medspan, m_inputSampleRate);
}

size_t
OnsetsDSPlugin::getPreferredStepSize() const
{
    return 256;
}

size_t
OnsetsDSPlugin::getPreferredBlockSize() const
{
    return 512;
}

OnsetsDSPlugin::OutputList
OnsetsDSPlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor onsets;
    onsets.identifier = "onsets";
    onsets.name = "Note Onsets";
    onsets.description = "Note onset positions";
    onsets.unit = "";
    onsets.hasFixedBinCount = true;
    onsets.binCount = 0;
    onsets.sampleType = OutputDescriptor::VariableSampleRate;
    onsets.sampleRate = (m_inputSampleRate / m_stepSize);

    list.push_back(onsets);

    return list;
}

OnsetsDSPlugin::FeatureSet
OnsetsDSPlugin::process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp)
{
    if (!m_ods) {
	cerr << "ERROR: OnsetsDSPlugin::process: Plugin has not been initialised"
	     << endl;
	return FeatureSet();
    }

    // We can const_cast because we happen to know onsetsds_process
    // does not modify this buffer
    bool result = onsetsds_process(m_ods, const_cast<float *>(inputBuffers[0]));

    FeatureSet returnFeatures;

    if (result) {
        Feature feature;
        feature.hasTimestamp = true;
        feature.timestamp = timestamp;
        returnFeatures[0].push_back(feature); // onsets are output 0
    }

    return returnFeatures;
}

OnsetsDSPlugin::FeatureSet
OnsetsDSPlugin::getRemainingFeatures()
{
    return FeatureSet();
}


static Vamp::PluginAdapter<OnsetsDSPlugin> adapter;

const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int version,
                                                    unsigned int index)
{
    if (version < 1) return 0;

    switch (index) {
    case  0: return adapter.getDescriptor();
    default: return 0;
    }
}

