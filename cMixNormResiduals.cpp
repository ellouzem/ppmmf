#include "StdAfxRegArchLib.h"

/*!
 \file cMixNormResiduals.cpp
 \brief implementation of the class for N(0, 1) conditional distribution.

 \author Jean-Baptiste DURAND, Ollivier TARAMASCO
 \date dec-18-2006 - Last change feb-18-2011
*/
namespace RegArchLib {
	/*!
	 * \fn cMixNormResiduals::cMixNormResiduals(const cDVector* theDistrParameter, bool theSimulFlag):cAbstResiduals(eNormal, NULL, theSimulFlag)
	 * \param bool theSimulFlag: true if created for simulation
	 * \details: mvBool is initialised by cAbstResiduals constructor and theDistrParameter is never used.
	 */
	cMixNormResiduals::cMixNormResiduals(const cDVector* theDistrParameter, bool theSimulFlag):cAbstResiduals(eNormal, NULL, theSimulFlag)
	{
		int myNParam = theDistrParameter.GetSize();
		if ( myNParam != 3)
			throw cError("Bad number of distribution parameters") ;
		mvMixNorm.ReAlloc(myNParam) ;

		for (register uint t = 0 ; t < myNParam ; t++){
				mvMixNorm[i] = mDistrParameter[i] ;
		}
		MESS_CREAT("cMixNormResiduals")
	}

	/*!
	 * \fn cMixNormResiduals::~cMixNormResiduals
	 * \details: nothing to do here
	 */
	cMixNormResiduals::~cMixNormResiduals()
	{
		MESS_DESTR("cMixNormResiduals")
	}

	/*!
	 * \fn cAbstCondVar* cMixNormResiduals::::PtrCopy()
	 */

	cAbstResiduals* cMixNormResiduals::PtrCopy() const
	{
		cStudentResiduals *mycStudentResiduals = NULL ;
		cDVector* myDistrParameter = new cDVector(mDistrParameter) ;

		bool mySimulFlag = false ;

		if (mtR != NULL)
			mySimulFlag = true ;

		mycStudentResiduals = new cStudentResiduals(myDistrParameter, mySimulFlag);

		delete myDistrParameter ;

		return mycStudentResiduals ;
	}

	/*!
	 * \fn void cMixNormResiduals::Print(ostream& theOut) const
	 * \param ostream& theOut: the output stream, default cout.
	 */
	void cMixNormResiduals::Print(ostream& theOut) const
	{
		theOut << "Conditional standardized normal distribution" << endl ;

	}
	/*!
	 * \fn void cMixNormResiduals::Generate(uint theNSample, cDVector& theYt) const
	 * \brief Draw a sample of N(0, 1) residuals.
	 * \param uint theNSample: the sample size.
	 * \param cDVector& theYt: output parameter
	 * \details: theYt is reallocated to size theNSample here.
	 */


	 void cMixNormResiduals::Set(double theValue, uint theIndex)
 	{
 		if (theIndex >= mvMixNorm.GetSize())
 			throw cError("Bad index") ;
 		else
 			mvMixNorm[theIndex]=theValue ;
 	}

	double cMixNormResiduals::Get(uint theIndex)
	{
		return mvMixNorm[theIndex] ;
	}


	void cMixNormResiduals::Generate(uint theNSample, cDVector& theYt) const
	{
		theYt.ReAlloc(theNSample) ;
		double e;

		for (register uint t = 0 ; t < theNSample ; t++)
			e = gsl_ran_bernoulli(mtR, mvMixNorm[0]);
			theYt[t] = e*gsl_ran_gaussian(mtR, mvMixNorm[1]) + (1-e)*gsl_ran_gaussian(mtR, mvMixNorm[2]);
	}

	/*!
	 * \fn double cMixNormResiduals::LogDensity(double theX) const
	 * \param double theX: the point where density is computed.
	 * \brief Compute the log density of N(0, 1)
	 */
	double cMixNormResiduals::LogDensity(double theX) const
	{
		double logn1 = (-LOG_SQRT_2_PI - log(mvMixNorm[1])- theX*theX/2.0);
		double logn2 = (-LOG_SQRT_2_PI - log(mvMixNorm[2])- theX*theX/2.0);
		double p = mvMixNorm[0];
		return p * logn1 + (1-p) * logn2;
	}

	/*!
	 * \fn double cMixNormResiduals::GetNParam(void)
	 * \param void.
	 * \brief return 0: no parameter for N(0,1) residuals.
	 */
	uint cMixNormResiduals::GetNParam(void) const
	{
		return 3 ;
	}

	/*!
	 * \fn static void GradLogDensity(double theX, cDVector& theGrad)
	 * \brief Compute the derivative of log density of a Gaussian distribution with respect to the random variable (theGrad[0])
	 * \e and the gradient of log density with respect to the model parameters (other components in theGrad)
	 * \param theX double: value of the random variable
	 * \param theGrad cDVector&: concatenation of derivatives with respect to the random variable and the model parameters
	 */
	static void GradLogDensity(double theX, cDVector& theGrad)
	{
		double logn1 = (-LOG_SQRT_2_PI - mvMixNorm[1]- theX*theX/2.0);
		double logn2 = (-LOG_SQRT_2_PI - mvMixNorm[2]- theX*theX/2.0);
		theGrad[0] = - theX ;
		theGrad[1] = logn1 - logn2;
		theGrad[2] = - mvMixNorm[0]/mvMixNorm[1];
		theGrad[3] = - (1-mvMixNorm[0])/mvMixNorm[2];
	}

	/*!
	 * \fn void cMixNormResiduals::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData)
	 * \brief Compute the derivative of log density with respect to the random variable (theGradData[0]) \e and the gradient
	 * of log density with respect to the model parameters (other components in theGradData)
	 * \param theDate uint: time at which gradient is computed
	 * \param theValue const cRegArchValue&: value of the random variable
	 * \param theGradData cRegArchGradient&: concatenation of derivatives with respect to the random variable and the model parameters
	 */
	void cMixNormResiduals::ComputeGrad(uint theDate, const cRegArchValue& theValue, cRegArchGradient& theGradData)
	{
		GradLogDensity(theValue.mEpst[theDate], theGradData.mCurrentGradDens) ;
	}

	void cMixNormResiduals::RegArchParamToVector(cDVector& theDestVect, uint theIndex) const
	{
	}

	void cMixNormResiduals::VectorToRegArchParam(const cDVector& theSrcVect, uint theIndex)
	{
	}


}//namespace
