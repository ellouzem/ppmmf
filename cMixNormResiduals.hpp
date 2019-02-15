#pragma once
#ifndef _CNORMRESIDUALS_H_
#define _CNORMRESIDUALS_H_
#include "cAbstResiduals.h"
#include "cRegArchValue.h"
#include "cRegArchGradient.h"
/*!
 \file cNormResiduals.h
 \brief Definition of the class for N(0, 1) conditional distribution.

 \author Jean-Baptiste DURAND, Ollivier TARAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/

namespace RegArchLib {

	/*!
	 * \class cNormResiduals
	 * \brief  Class to implement the N(0, 1) residuals
	 */

	class _DLLEXPORT_ cNormResiduals: public cAbstResiduals
	{
	private:
		cDVector mvMixNorm ; // vetcor containing the parameters of the distribution p, sigma1, sigma2,
	public :

		cNormResiduals(const cDVector* theDistrParameter=NULL, bool theSimulFlag=true) ; ///< a simple constructor
		virtual ~cNormResiduals() ; ///< a simple destructor
		virtual cAbstResiduals* PtrCopy() const ; /// < Return a copy of *this
		void Print(ostream& theOut) const ; ///< print the distribution type
		void Generate(uint theNSample, cDVector& theYt) const ; ///< Draw a sample from residual distribution
		double LogDensity(double theX) const ;
		uint GetNParam(void) const ;
		void Set(double theValue, uint theIndex=0);
		double Get(uint theIndex=0);
		/** Compute the derivatives of the log-density with respect to the variable \e and the parameters */
		void ComputeGrad(uint theDate, const cRegArchValue& theData, cRegArchGradient& theGradData) ;
		void RegArchParamToVector(cDVector& theDestVect, uint theIndex) const ;
		void VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex = 0) ;
	} ;

}

#endif //_CNORMRESIDUALS_H_
