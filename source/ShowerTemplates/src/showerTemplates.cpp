/*
 * showerTemplates.cxx
 *
 *  Created on: Jul 1, 2014
 *      Author: strahinja
 */

#ifndef __CINT__
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#endif


#include "showerTemplates.hh"


/*************************
 *
 * Implementation of 1D templates
 *
 *************************/


Double_t eShower1D_varx0 (Double_t *x, Double_t *par)
{
	Double_t k = par[0];
	Double_t a = par[1];
	Double_t b = par[2];
	Double_t x0 = par[3];

	Double_t s = x[0] - x0;

	static const double lowedge = 0.5;
	if(x[0]<lowedge) return 0;
//	static const double mipdepth = 0.00;
	double mipdepth = pow(.001/k, 1./(a-1));
	if(s<=mipdepth) s=mipdepth;//*mipdepth/(mipdepth+x0-lowedge);

	return k * (pow(s, a-1.) * exp(-b*s));
}



ProfileTester::ProfileTester(double eCoeff, double p0a, double p1a, double p0b, double p1b, double loEdge, double hiEdge) :
		TF1("ShowerTemplate", eShower1D_varx0, loEdge, hiEdge, 4),
		_eCoeff(eCoeff), _p0a(p0a), _p1a(p1a), _p0b(p0b), _p1b(p1b)
{
	InitParameters(100.);
}


ProfileTester::~ProfileTester() {}

void ProfileTester::Calibrate(double eCoeff, double p0a, double p1a, double p0b, double p1b)
{
	_eCoeff = eCoeff;
	_p0a = p0a;
	_p1a = p1a;
	_p0b = p0b;
	_p1b = p1b;
}

void ProfileTester::InitParameters(double deposit)
{
	if (!(deposit > 1.e-6)) deposit = 1.e-6;
	double energy = _eCoeff*deposit;
	double a = _p0a*log(_p1a*energy);
	double b = _p0b*log(_p1b*energy);

	FixParameter(1, a);
	FixParameter(2, b);
	SetParameter(3, 0);
	// Normalize to PDF
	FixParameter(0,  pow(b, a) / TMath::Gamma(a) );
}


double ProfileTester::Test(const TH1 * profile, float &xStart, float &maxCorrel)
{
	TH1* profileNormed = dynamic_cast< TH1D*> (profile->Clone("normed"));
	if(profileNormed == 0)
	{
		profileNormed = dynamic_cast< TH1F*> (profile->Clone("normed"));
		if(profileNormed == 0)
		{
			profileNormed = dynamic_cast< TH1I*> (profile->Clone("normed"));
			if(profileNormed == 0) return 0.;
		}
	}

	double deposit = profileNormed->Integral();
	if(!(deposit > 1.e-6)) return 0.;

	InitParameters(deposit); // Energy dependence of the ideal shape

	xStart = -5.; maxCorrel = 0.;
	// template function autocorrelation (scalar product with self)
	float autoFsqBest = 0;
	for(float x0 = -5.; x0<25.; x0 += .01)
	{
		float correl=0., autoFsq=0.;
		SetParameter(3,x0);

		for(int i=1; i<profileNormed->GetNbinsX(); i++)
		{
			float fi = Eval(float(i));
			float hi = profileNormed->GetBinContent(i);
			correl += fi * hi;
			autoFsq += fi*fi;
		}

		if (correl > maxCorrel)
		{ maxCorrel = correl; xStart = x0; autoFsqBest = autoFsq;}
	}

	// Histogram autocorrelation (scalar product with self)
	float  autoHsq=0.;
	for(int i=1; i<profileNormed->GetNbinsX(); i++)
	{
		float hi = profileNormed->GetBinContent(i);
		autoHsq += hi*hi;
	}


	delete profileNormed;
	// Return true correlation coefficient
	maxCorrel /= (sqrt(autoFsqBest) * sqrt(autoHsq));
	return maxCorrel;
}


void ProfileTester::Draw(double integral, double x0, Option_t *option)
{
	double a, b, k;
	a = GetParameter(1);
	b = GetParameter(2);
	k = pow(b, a) / TMath::Gamma(a);
	SetParameter(0, k*integral);

	SetParameter(3, x0);

	TF1::Draw(option);
}

