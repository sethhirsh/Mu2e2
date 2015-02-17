#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <vector>
#include "FitModel.C"
#include "config.hh"

typedef unsigned int * adcWaveform;

struct resultantPeakData
{
	Float_t _peakHeight; // in units of bits
	Float_t _peakTime; // time of peak relative to 140.0 ns interval

	// Default constructor - should probably be deleted
	resultantPeakData() : _peakTime(0.0), _peakHeight(0.0){};

	// True constructor
	resultantPeakData(Float_t peakTime, Float_t peakHeight) : _peakTime(peakTime), _peakHeight(peakHeight){};
};

// This is object top which will be filled by the process method 
typedef std::vector<resultantPeakData> resultantHitData;


// Virtual class providing structure for FindSinglePeak, FindDoublePeak, FindMutiplePeaks, etc. 
class FindPeakBase{
	public:
		
		// Fills result using adc waveform data
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result) = 0;

		// Destructor
		virtual ~FindPeakBase(){}

		// Default Constructor
		FindPeakBase(){}

		// FindPeakBase normal constructor with configStruct initilization parameters
		FindPeakBase(const configStruct &initParams) : _initParams(initParams){}

	protected:

		// UNDERSCORE MEMBER VARIABLES

		const configStruct _initParams; 

		// These should probably change from Float_t to Int_t (or unsigned int)


		// Fits a model function to a waveform
		void fitModel2Waveform(TF1 &fitModel, TGraphErrors &fitData, const Double_t *initialParameters, Double_t *fitParameters)
		{
 			// These lines will be replaced with the chi-square minimization
			TF1 *fitModelPtr = &fitModel; 
			TGraphErrors *fitDataPtr = &fitData;
			fitModel.SetParameters(initialParameters);
			fitDataPtr->Fit(fitModelPtr,"QN");

			const Int_t numParameters = fitModel.GetNumberFreeParameters();
			std::cout << "numParam : " << numParameters << std::endl;
			for (int i = 0; i < numParameters; ++i)
			{
				fitParameters[i] = fitModel.GetParameter(i);
			}
		}

		// Converts adcWaveform object to TGraphErrors object for easier manipulation in ROOT
		void adcWaveform2TGraphErrors(adcWaveform adcData, TGraphErrors &fitData)
		{
			Double_t adcDataTemp[_initParams._numSamplesPerHit];
			Double_t measurementTimes[_initParams._numSamplesPerHit];
			Double_t measurementTimesErrors[_initParams._numSamplesPerHit];
			Double_t adcDataErrors[_initParams._numSamplesPerHit];

			for (int i = 0; i < _initParams._numSamplesPerHit; ++i)
			{
				adcDataTemp[i] = (Double_t) adcData[i];
				measurementTimes[i] = (Double_t) i * _initParams._measurementFrequency; 
				measurementTimesErrors[i] = 0.0;
				adcDataErrors[i] = _initParams._adcError;
			}

			fitData = TGraphErrors(_initParams._numSamplesPerHit,adcDataTemp,measurementTimes,measurementTimesErrors,adcDataErrors);
		}

		// Precomputed constants
		const Double_t _bits2scalingFactor = _initParams._shapingTime * TMath::E(); // approximately 67.96
		const Double_t _scalingFactor2bits  = 1.0 / _bits2scalingFactor; // approximately 0.0147
		const Double_t _hitPeriod = _initParams._numSamplesPerHit * (_initParams._measurementFrequency - 1.0);

};

class FindSinglePeakWithDynamicPedestal : public FindPeakBase{
	public:

		// FindSinglePeakWithDynamicPedestal normal constructor with configStruct initilization parameters
		FindSinglePeakWithDynamicPedestal(const configStruct &initParams) : FindPeakBase(initParams)
		{
			_fitModel = TF1("fitModel",convolutionSinglePeakWithDynamicPedestal,0.0,_hitPeriod,4);
		}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result)
		{

			// Set initial fit parameters
			const double timeshift = 30.0;
			const double scalingFactor = result[1]._peakHeight * _bits2scalingFactor;
			const double Q = result[0]._peakHeight;
			const double sigma = 10.0;
			Double_t initialParameters[4] = {timeshift, scalingFactor, Q, sigma};

			Double_t finalParameters[4];

			adcWaveform2TGraphErrors(adcData, _fitData);
			fitModel2Waveform(_fitModel, _fitData, initialParameters, finalParameters);
			fitParams2ResultantData(finalParameters, result);
		}

	private:
		TGraphErrors _fitData;
		TF1 _fitModel;

		void fitParams2ResultantData(Double_t *fitParameters, resultantHitData &result)
		{
			resultantPeakData peakData;

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			peakData._peakTime = fitParameters[0];
			peakData._peakHeight = fitParameters[1];

			result[0] = peakData;

			// Since dynamic pedestal is not counted as peak
			result.pop_back();
		}




};

class FindSinglePeak : public FindPeakBase{
	public:
		// FindSinglePeak normal constructor with configStruct initilization parameters
		FindSinglePeak(const configStruct &initParams) : FindPeakBase(initParams)
		{
		  _fitModel = TF1("fitModel",convolutionSinglePeakWithConstantPedestal,0.0,_hitPeriod,4);
		}


		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result)
		{	
			// Set initial fit parameters
			const double timeshift = 30.0;
			const double scalingFactor = TMath::Max(result[0]._peakHeight * _bits2scalingFactor, 1000.0);
			const double sigma = 10.0;
			double verticalShift = 0.0;
			const Double_t initialParameters[4] = {timeshift, scalingFactor, verticalShift, sigma};

			
			// DEAL WITH CASE OF ZERO PEDESTAL
			//const bool nonZeroPedestal = (qAdc[0] + qAdc[1])*0.5 > 4.0 / sqrt(2.0) * 3.0;
	
			Double_t finalParameters[4];

			adcWaveform2TGraphErrors(adcData, _fitData);
			fitModel2Waveform(_fitModel, _fitData, initialParameters, finalParameters);
			fitParams2ResultantData(finalParameters, result);
		}
	private:
		void fitParams2ResultantData(Double_t *fitParameters, resultantHitData &result)
		{
			resultantPeakData peakData;

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			peakData._peakTime = fitParameters[0];
			peakData._peakHeight = fitParameters[1];

			result[0] = peakData;
		}

		TGraphErrors _fitData;
		TF1 _fitModel;
};

class FindDoublePeak : public FindPeakBase{
	public:
		FindDoublePeak(const configStruct &initParams) : FindPeakBase(initParams)
		{
			_fitModel = TF1("fitModel",doublePeak,0.0,_hitPeriod,5);
		}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result)
		{
			// Set initial fit parameters
			const double timeShift0 = 30.0;
			const double scalingFactor0 = TMath::Max(result[0]._peakHeight / 0.015, 1000.0);
			double verticalShift = (adcData[0] + adcData[1]) * 0.5; // This has implicit casting
			const double timeshift1 = result[1]._peakTime - result[0]._peakTime;
			const double scalingFactor1 = TMath::Max(result[1]._peakHeight / 0.015, 1000.0);
			const Double_t initialParameters[5] = {timeShift0, scalingFactor0, verticalShift, timeshift1, scalingFactor1};

			Double_t finalParameters[5];

			adcWaveform2TGraphErrors(adcData, _fitData);
			fitModel2Waveform(_fitModel, _fitData, initialParameters, finalParameters);
			fitParams2ResultantData(finalParameters, result);
		}

	protected:

		void fitParams2ResultantData(Double_t *fitParameters, resultantHitData &result)
		{
			resultantPeakData firstPeakData;

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			firstPeakData._peakTime = fitParameters[0];
			firstPeakData._peakHeight = fitParameters[1];

			resultantPeakData secondPeakData;
			secondPeakData._peakTime = fitParameters[3];
			secondPeakData._peakHeight = fitParameters[4];

			result[0] = firstPeakData;
			result[1] = secondPeakData;
		}

	private:
		TGraphErrors _fitData;
		TF1 _fitModel;
};


class FindDoublePeakWithDynamicPedestal : public FindDoublePeak{
	public:
		// FindDoublePeakWithDynamicPedestal normal constructor with configStruct initilization parameters
		FindDoublePeakWithDynamicPedestal(const configStruct &initParams) : FindDoublePeak(initParams)
		{
			_fitModel = TF1("fitModel",doublePeakWithDynamicPedestal,0.0,_hitPeriod,5);
		}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result)
		{
			const double timeShift0 = 30.0;
			const double scalingFactor0 = TMath::Max(result[1]._peakHeight * _bits2scalingFactor, 1000.0);
			const double Q = result[0]._peakHeight;
			const double timeshift1 = result[2]._peakTime - result[1]._peakTime;
			const double scalingFactor1 = TMath::Max(result[2]._peakHeight * _bits2scalingFactor, 1000.0);

			Double_t initialParameters[5] = {timeShift0, scalingFactor0, Q, timeshift1, scalingFactor1};
			Double_t finalParameters[5];

			adcWaveform2TGraphErrors(adcData, _fitData);
			fitModel2Waveform(_fitModel, _fitData, initialParameters, finalParameters);
			fitParams2ResultantData(finalParameters, result);
		}

	private:
		TGraphErrors _fitData;
		TF1 _fitModel;


		void fitParams2ResultantData(Double_t *fitParameters, resultantHitData &result)
		{
			resultantPeakData firstPeakData;

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			firstPeakData._peakTime = fitParameters[0];
			firstPeakData._peakHeight = fitParameters[1];

			resultantPeakData secondPeakData;
			secondPeakData._peakTime = fitParameters[3];
			secondPeakData._peakHeight = fitParameters[4];

			result[0] = firstPeakData;
			result[1] = secondPeakData;
			result.pop_back();
		}
};




class FindMultiplePeaks : public FindPeakBase{
	public:

		// FindMultiplePeaks normal constructor with configStruct initilization parameters
		FindMultiplePeaks(const configStruct &initParams) : FindPeakBase(initParams){}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result)
		{
			adcWaveform2TGraphErrors(adcData, _fitData);
			findPeaks(_fitData, initParams, result);

			const int numPeaks = result.size();

			if (result[0]._peakTime == 0.0) // If there is a dynamic pedestal
			{

				// If the only peak is a dynamic pedestal 
				// search for another peak
				if (numPeaks == 1) 
				{
					dynamicPedestalAddPeak(_fitData, result);
					FindSinglePeakWithDynamicPedestal singlePeak(initParams);
					singlePeak.process(adcData, result);
				}
				else if (numPeaks == 2)
				{
					FindSinglePeakWithDynamicPedestal singlePeak(initParams);
					singlePeak.process(adcData, result);
				}
				else if (numPeaks == 3)
				{
					FindDoublePeakWithDynamicPedestal doublePeak(initParams);
					doublePeak.process(adcData, result);
				}
			}
			// If there is no dynamic pedestal
			else
			{
				if (numPeaks == 1)
				{
					FindSinglePeak singlePeak(initParams);
					singlePeak.process(adcData, result);
				}
				if (numPeaks == 2)
				{
					FindDoublePeak doublePeak(initParams);
					doublePeak.process(adcData, result);
				}

			}	
		}

	private:
		// Performs explicit peak search on adc waveform data
		void findPeaks(TGraphErrors &gr, const configStruct &initParams, resultantHitData &result, double sigma = 3.0)
		{
  			int ientry = 0; // Start time at 0
  			const double *measurementTimes = gr.GetX();
  			const double *adcValues = gr.GetY();

  			while(ientry < FindPeakBase::_initParams._numSamplesPerHit)
  			{
    			double adcValue = adcValues[ientry];
				double tMax = measurementTimes[ientry];
    			double adcMax = adcValue;
    			double adcPrev = adcValue;

    			int jentry = ientry + 1;
    			bool descending = false;
    			while (jentry < FindPeakBase::_initParams._numSamplesPerHit)
				{
					adcValue = adcValues[jentry];
      				descending |= ((adcPrev-adcValue) > (TMath::Sqrt2()*FindPeakBase::_initParams._adcError*sigma));

      				if (descending && (adcValue-adcPrev > (TMath::Sqrt2()*FindPeakBase::_initParams._adcError*sigma)))
      				{
        			break;
      				}
      				else
      				{
        			if (adcValue > adcMax)
        			{
          			adcMax  = adcValue;
          			tMax = measurementTimes[jentry];
      				}
        			adcPrev = adcValue;
        			ientry = jentry;
        			++jentry;
        			}
  				}
  				resultantPeakData peakData(tMax, adcMax);
      			result.push_back(peakData);
      			++ientry;
  			}
		}

		// This function searches for another peak in the waveform data by subtracting out a dynamic pedestal 
		// from the adc waveform and finding the maximum adc value in the "subtracted data".
		// This function is applied when no peak is found in the explicit peak search (findPeaks).
		void dynamicPedestalAddPeak(TGraphErrors &gr, resultantHitData &result)
		{	
			// This maybe could be done using linear algebra vectors
			// instead of arrays
			const Double_t *adcValues = gr.GetY();
			const Double_t *measurementTimes = gr.GetX();
			Double_t subtractedValues[FindPeakBase::_initParams._numSamplesPerHit];

			Double_t dynamicPedstalParam[1] = {adcValues[0]};
			Double_t dynamicPedestalX[1];

			for (int i = 0; i < FindPeakBase::_initParams._numSamplesPerHit; ++i)
			{
				dynamicPedestalX[0] = measurementTimes[i];
				subtractedValues[i] = adcValues[i] - dynamicPedestal(dynamicPedestalX, dynamicPedstalParam);
			}

			// New peak is max value of difference between of adc values and dynamic pedestal
			const Float_t newAdcPeak = TMath::MaxElement(FindPeakBase::_initParams._numSamplesPerHit, subtractedValues);
			const Float_t newTPeak = TMath::LocMax(FindPeakBase::_initParams._numSamplesPerHit, subtractedValues);

			resultantPeakData newPeakData(newTPeak, newAdcPeak);
			result.push_back(newPeakData);
		}

		TGraphErrors _fitData;
		TF1 _fitModel;


};

