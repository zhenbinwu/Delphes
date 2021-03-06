/** \class ParticlePropagator
 *
 *  Propagates charged and neutral particles
 *  from a given vertex to a cylinder defined by its radius, 
 *  its half-length, centered at (0,0,0) and with its axis
 *  oriented along the z-axis.
 *
 *  $Date: 2013-02-12 15:13:59 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 907 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/ParticlePropagator.h"

//#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

static const double mm  = 1.;
static const double m = 1000.*mm;
static const double ns  = 1.;
static const double s = 1.e+9 *ns;
static const double c_light   = 2.99792458e+8 * m/s;


using namespace std;

//------------------------------------------------------------------------------

ParticlePropagator::ParticlePropagator() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

ParticlePropagator::~ParticlePropagator()
{
}

//------------------------------------------------------------------------------

void ParticlePropagator::Init()
{
  fRadius = GetDouble("Radius", 1.0);
  fRadius2 = fRadius*fRadius;
  fHalfLength = GetDouble("HalfLength", 3.0);
  fBz = GetDouble("Bz", 0.0);
  fKeepPileUp = GetInt("KeepPileUp",1);

  if(fRadius < 1.0E-2)
  { 
    cout << "ERROR: magnetic field radius is too low\n";
    return;
  }
  if(fHalfLength < 1.0E-2)
  {
    cout << "ERROR: magnetic field length is too low\n";
    return;
  }

  // import array with output from filter/classifier module

  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
  fChargedHadronOutputArray = ExportArray(GetString("ChargedHadronOutputArray", "chargedHadrons"));
  fElectronOutputArray = ExportArray(GetString("ElectronOutputArray", "electrons"));
  fMuonOutputArray = ExportArray(GetString("MuonOutputArray", "muons"));
}

//------------------------------------------------------------------------------

void ParticlePropagator::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ParticlePropagator::Process()
{
  Candidate *candidate, *mother;
  TLorentzVector candidatePosition, candidateMomentum;
  Double_t px, py, pz, pt, pt2, e, q;
  Double_t x, y, z, t, r, phi;
  Double_t x_c, y_c, r_c, phi_c, phi_0;
  Double_t x_t, y_t, z_t, r_t;
  Double_t t1, t2, t3, t4, t5, t6;
  Double_t t_z, t_r, t_ra, t_rb;
  Double_t tmp, discr, discr2;
  Double_t delta, gammam, omega, asinrho;

  float t_orig;
  
  //  const Double_t c_light = 2.99792458E8;
    
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {

    // For no pileup samples...
    if (fKeepPileUp == 0 && candidate->IsPU > 0) {
      continue;
    }

    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;
    x = candidatePosition.X()*1.0E-3;
    y = candidatePosition.Y()*1.0E-3;
    z = candidatePosition.Z()*1.0E-3;
    t_orig = candidatePosition.T(); // this should already be in ns; if not, fix it!
    q = candidate->Charge;

    // check that particle position is inside the cylinder
    if(TMath::Hypot(x, y) > fRadius || TMath::Abs(z) > fHalfLength)
    {
      continue;
    }

    px = candidateMomentum.Px();
    py = candidateMomentum.Py();
    pz = candidateMomentum.Pz();
    pt = candidateMomentum.Pt();
    pt2 = candidateMomentum.Perp2();
    e = candidateMomentum.E();

    if(pt2 < 1.0E-9)
    {
      continue;
    }

    if(TMath::Abs(q) < 1.0E-9 || TMath::Abs(fBz) < 1.0E-9)
    {
      // solve pt2*t^2 + 2*(px*x + py*y)*t + (fRadius2 - x*x - y*y) = 0
      tmp = px*y - py*x;
      discr2 = pt2*fRadius2 - tmp*tmp;
      
      if(discr2 < 0)
      {
        // no solutions
        continue;
      }

      tmp = px*x + py*y;
      discr = TMath::Sqrt(discr2);
      t1 = (-tmp + discr)/pt2;
      t2 = (-tmp - discr)/pt2;
      t = (t1 < 0) ? t2 : t1; 

      z_t = z + pz*t;
      if(TMath::Abs(z_t) > fHalfLength)
      {
        t3 = (+fHalfLength - z) / pz;
        t4 = (-fHalfLength - z) / pz;
        t = (t3 < 0) ? t4 : t3; 
      }

      x_t = x + px*t;
      y_t = y + py*t;
      z_t = z + pz*t;

      mother = candidate;
      candidate = static_cast<Candidate*>(candidate->Clone());

      // In the neutral part of the code, t has units of m / GeV
      // This is some parameterization of the path, but not the time in any way that's obvious to me (S. Zenz)
      // Probably easiest to recompute the time from the (straight-line) distance travelled...
      t = 1E9*(1/c_light)*TMath::Sqrt((x_t-x)*(x_t-x)+(y_t-y)*(y_t-y)+(z_t-z)*(z_t-z)); // 1E9/c_light converts meters --> ns
      //      if (pt > 1.) {
      //	cout << " SCZ NEUTRAL Debug t t_orig " << t << " " << t_orig << endl;
      //      }
      t += t_orig;

      candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, t);

      candidate->Momentum = candidateMomentum;
      candidate->AddCandidate(mother);
      
      /*
      cout << "In ParticlePropogator. Have just added neutral candidate with X Y Z T " << candidate->Position.X() << " " << candidate->Position.Y() << " "
	   << candidate->Position.Z() << " " << candidate->Position.T() << " " << endl;
      Candidate *prt = static_cast<Candidate*>(candidate->GetCandidates()->Last());
      const TLorentzVector &ini = prt->Position;
      cout << "                                                   Mother has X Y Z T " << ini.X() << " " << ini.Y() << " " << ini.Z() << " " << ini.T() << endl;
      */


      fOutputArray->Add(candidate);
      if(TMath::Abs(q) > 1.0E-9) 
      {
        switch(TMath::Abs(candidate->PID))
        {
          case 11:
            fElectronOutputArray->Add(candidate);
            break;
          case 13:
            fMuonOutputArray->Add(candidate);
            break;
          default:
            fChargedHadronOutputArray->Add(candidate);
        }
      }
    }
    else
    {

      // 1.  initial transverse momentum p_{T0} : Part->pt
      //     initial transverse momentum direction \phi_0 = -atan(p_X0/p_Y0) 
      //     relativistic gamma : gamma = E/mc� ; gammam = gamma \times m
      //     giration frequency \omega = q/(gamma m) fBz
      //     helix radius r = p_T0 / (omega gamma m)

      gammam = e*1.0E9 / (c_light*c_light);      // gammam in [eV/c�]
      omega = q * fBz / (gammam);                // omega is here in [ 89875518 / s]  
      r = pt / (q * fBz) * 1.0E9/c_light;        // in [m]  

      phi_0 = TMath::ATan2(py, px); // [rad] in [-pi; pi]

      // 2. helix axis coordinates
      x_c = x + r*TMath::Sin(phi_0);
      y_c = y - r*TMath::Cos(phi_0);
      r_c = TMath::Hypot(x_c, y_c);
      phi_c = TMath::ATan2(y_c, x_c);
      phi = phi_c;
      if(x_c < 0.0) phi += TMath::Pi();

      // 3. time evaluation t = TMath::Min(t_r, t_z)
      //    t_r : time to exit from the sides
      //    t_z : time to exit from the front or the back
      t_r = 0.0; // this is in s in this part of the code, I think...?
      int sign_pz = (pz > 0.0) ? 1 : -1;
      if(pz == 0.0) t_z = 1.0E99;
      else t_z = gammam / (pz*1.0E9/c_light) * (-z + fHalfLength*sign_pz);

      if(r_c + TMath::Abs(r)  < fRadius)
      {
        // helix does not cross the cylinder sides
        t = t_z;
      }
      else
      {
        asinrho = TMath::ASin( (fRadius*fRadius - r_c*r_c - r*r) / (2*TMath::Abs(r)*r_c)  );
        delta = phi_0 - phi;
        if(delta <-TMath::Pi()) delta += 2*TMath::Pi();
        if(delta > TMath::Pi()) delta -= 2*TMath::Pi();
        t1 = (delta + asinrho) / omega;
        t2 = (delta + TMath::Pi() - asinrho) / omega;
        t3 = (delta + TMath::Pi() + asinrho) / omega;
        t4 = (delta - asinrho) / omega;
        t5 = (delta - TMath::Pi() - asinrho) / omega;
        t6 = (delta - TMath::Pi() + asinrho) / omega;

        if(t1 < 0) t1 = 1.0E99;
        if(t2 < 0) t2 = 1.0E99;
        if(t3 < 0) t3 = 1.0E99;
        if(t4 < 0) t4 = 1.0E99;
        if(t5 < 0) t5 = 1.0E99;
        if(t6 < 0) t6 = 1.0E99;

        t_ra = TMath::Min(t1, TMath::Min(t2, t3));
        t_rb = TMath::Min(t4, TMath::Min(t5, t6));
        t_r = TMath::Min(t_ra, t_rb);
        t = TMath::Min(t_r, t_z); 
      }

      // 4. position in terms of x(t), y(t), z(t)
      x_t = x_c + r * TMath::Sin(omega * t - phi_0);
      y_t = y_c + r * TMath::Cos(omega * t - phi_0);
      z_t = z + pz*1.0E9 / c_light / gammam * t;
      r_t = TMath::Hypot(x_t, y_t);

      if(r_t > 0.0)
      {
        mother = candidate;
        candidate = static_cast<Candidate*>(candidate->Clone());

	//	if (pt > 1.) {
	  //			  cout << " SCZ CHARGED Debug: t_r t_z t 1E9*t t_orig " << t_r << " " << t_z << " " << t << " "  << 1E9*t << " " << t_orig << endl;
	  //	                  cout << "   SCZ Debug line2: pt pz x y z " << pt << " " << pz << " " << x << " " << y << " " << z << endl;
	  //	}
	t = 1E9*t + t_orig; // Delphes had omega in inverse s, so 1E9 converts to ns

	//	cout << " Original (z,t)=(" << z << ", " << t_orig << ") - New (z,t)=(" << z_t*1.0E3 << ", " << t << ")" << endl;

        candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, t);

        candidate->Momentum = candidateMomentum;
        candidate->AddCandidate(mother);

	/*
	cout << "In ParticlePropogator. Have just added charged candidate with X Y Z T " << candidate->Position.X() << " " << candidate->Position.Y() << " "
	     << candidate->Position.Z() << " " << candidate->Position.T() << " " << endl;
	Candidate *prt = static_cast<Candidate*>(candidate->GetCandidates()->Last());
	const TLorentzVector &ini = prt->Position;
	cout << "                                                   Mother has X Y Z T " << ini.X() << " " << ini.Y() << " " << ini.Y() << " " << ini.Z() << endl;
	*/

        fOutputArray->Add(candidate);
        switch(TMath::Abs(candidate->PID))
        {
          case 11:
            fElectronOutputArray->Add(candidate);
	    /*
	    if (pt > 0.5) {
  	      double t_straight = t_orig + 1E9*(1/c_light)*TMath::Sqrt((x_t-x)*(x_t-x)+(y_t-y)*(y_t-y)+(z_t-z)*(z_t-z)); // 1E9/c_light converts meters --> ns
              cout << " SCZ Debug ParticlePropogator pt t_orig t_straight t IsPU PID " << pt << " " << t_orig << " " << t_straight << " "
		   << t << " " << candidate->IsPU << " " << candidate->PID << endl;
	    }
	    */
            break;
          case 13:
            fMuonOutputArray->Add(candidate);
            break;
          default:
            fChargedHadronOutputArray->Add(candidate);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
