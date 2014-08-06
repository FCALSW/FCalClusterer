/*
 * showerTemplates.h
 *
 *  Created on: Jul 1, 2014
 *      Author: strahinja
 */

#ifndef SHOWERTEMPLATES_H_
#define SHOWERTEMPLATES_H_

#ifndef __CINT__
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#endif




/*********************************
 *
 * 1D templates
 *
 *********************************/


Double_t eShower1D_varx0 (Double_t *x, Double_t *par);



class ProfileTester : protected TF1
{
private:
	double _eCoeff, _p0a, _p1a, _p0b, _p1b;
	void InitParameters(double deposit);
public:
	ProfileTester(double eCoeff=56.72, double p0a=0.737, double p1a=1.75, double p0b=.02209, double p1b=1.8e6, double loEdge=0.5, double hiEdge=39.5);
	~ProfileTester();

	void Calibrate(double eCoeff, double p0a=0.737, double p1a=1.75, double p0b=.02209, double p1b=1.8e6);

	double Test( const TH1* profile, float &x0, float &maxCorrel) ;
	void Draw(double integral, double x0, Option_t *option="");

	double Get_k() const {return GetParameter(0);};
	double Get_norm() const {return Get_k();};
	double Get_a() const {return GetParameter(1);};
	double Get_b() const {return GetParameter(2);};
	double Get_x0() const {return GetParameter(3);};
};



#endif /* SHOWERTEMPLATES_H_ */
